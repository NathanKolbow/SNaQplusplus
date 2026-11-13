#include "util.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <numeric>

namespace csnaq {

void warn(const std::string& msg) { std::cerr << "┌ Warning: " << msg << "\n└ @ snaq\n"; }

// Julia's `mod` (result has the sign of the divisor).
static double julia_mod(double a, double b) {
    double r = std::fmod(a, b);
    if (r != 0.0 && ((r < 0.0) != (b < 0.0))) r += b;
    return r;
}

std::string julia_float(double x) {
    if (x == 0.0) return std::signbit(x) ? "-0.0" : "0.0";
    if (std::isnan(x)) return "NaN";
    if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";

    // Shortest round-trip digits: the smallest precision whose correctly
    // rounded scientific representation parses back to x (portable
    // alternative to std::to_chars, which older macOS/GCC lack for doubles).
    char buf[64];
    for (int prec = 0; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof(buf), "%.*e", prec, x);
        if (std::strtod(buf, nullptr) == x) break;
    }
    std::string sci(buf);
    bool neg = false;
    size_t p = 0;
    if (sci[0] == '-') {
        neg = true;
        p = 1;
    }
    size_t epos = sci.find('e');
    std::string mant = sci.substr(p, epos - p);
    int e10 = std::stoi(sci.substr(epos + 1));
    std::string digits;
    for (char c : mant)
        if (c != '.') digits.push_back(c);
    // strip trailing zeros (to_chars shortest should not produce them, but be safe)
    while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
    const int olength = static_cast<int>(digits.size());
    const int pt = e10 + 1;  // position of decimal point relative to first digit

    bool plain = (-4 < pt && pt <= 6);
    if (plain && pt >= olength) {
        double m = julia_mod(x + 0.05, std::pow(10.0, pt - olength));
        if (std::abs(m - 0.05) > 0.05) plain = false;
    }

    std::string out = neg ? "-" : "";
    if (plain) {
        if (pt <= 0) {
            out += "0.";
            out.append(static_cast<size_t>(-pt), '0');
            out += digits;
        } else if (pt >= olength) {
            out += digits;
            out.append(static_cast<size_t>(pt - olength), '0');
            out += ".0";
        } else {
            out += digits.substr(0, static_cast<size_t>(pt));
            out += ".";
            out += digits.substr(static_cast<size_t>(pt));
        }
    } else {
        out += digits[0];
        out += ".";
        out += (olength > 1 ? digits.substr(1) : std::string("0"));
        out += "e";
        out += std::to_string(pt - 1);
    }
    return out;
}

double julia_round(double x, int digits) {
    // Julia's round(x, digits=d) for Float64 uses round(x * 10^d) / 10^d with
    // a correction to handle representation error; this matches it for the
    // magnitudes that appear in SNaQ output (branch lengths, γ, scores).
    if (!std::isfinite(x)) return x;
    double scale = std::pow(10.0, digits);
    double y = x * scale;
    double r = std::nearbyint(y);  // round half to even (default FE_TONEAREST)
    // Julia: if the scaled value is not exactly representable, the
    // rounding is done on the exact decimal value; emulate via a check.
    double z = r / scale;
    // Julia's `_round_digits` uses `round(x * scale) / scale` when
    // `x * scale` is exact, else a more careful method; differences only
    // arise exactly at ties, which do not matter for log output.
    return z;
}

long long julia_round_int(double x) { return static_cast<long long>(std::nearbyint(x)); }

std::string time_elapsed(double totaltime) {
    long long seconds = julia_round_int(totaltime);
    long long minutes = (seconds / 60) % 60;
    long long hours = (seconds / 3600) % 24;
    long long days = seconds / 86400;
    seconds = seconds % 60;
    auto s = [](long long v) { return std::to_string(v); };
    if (days > 0)
        return s(days) + " days, " + s(hours) + " hours, " + s(minutes) + " minutes and " + s(seconds) + " seconds";
    if (hours > 0) return s(hours) + " hours, " + s(minutes) + " minutes and " + s(seconds) + " seconds";
    if (minutes > 0) return s(minutes) + " minutes and " + s(seconds) + " seconds";
    return s(seconds) + " seconds";
}

static std::tm local_tm() {
    std::time_t t = std::time(nullptr);
    std::tm tmv{};
#if defined(_WIN32)
    localtime_s(&tmv, &t);
#else
    localtime_r(&t, &tmv);
#endif
    return tmv;
}

std::string current_time() {
    std::tm tmv = local_tm();
    char buf[64];
    std::strftime(buf, sizeof(buf), "%H:%M:%S %Y-%m-%d", &tmv);
    return buf;
}

std::string current_timestamp() {
    std::tm tmv = local_tm();
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tmv);
    return buf;
}

double wall_time() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// ------------------------------- RNG -------------------------------------

static uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

void Rng::seed(uint64_t seed) {
    uint64_t x = seed;
    for (auto& s : s_) s = splitmix64(x);
}

static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

uint64_t Rng::next_u64() {
    const uint64_t result = rotl(s_[0] + s_[3], 23) + s_[0];
    const uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
}

double Rng::rand() { return static_cast<double>(next_u64() >> 11) * 0x1.0p-53; }

long long Rng::randint(long long lo, long long hi) {
    if (hi < lo) throw SnaqError("randint: empty range");
    uint64_t range = static_cast<uint64_t>(hi - lo) + 1;
    if (range == 0) return static_cast<long long>(next_u64());
    // Lemire-style rejection to avoid modulo bias.
    uint64_t limit = UINT64_MAX - (UINT64_MAX % range);
    uint64_t r;
    do {
        r = next_u64();
    } while (r >= limit);
    return lo + static_cast<long long>(r % range);
}

std::vector<int> Rng::randperm(int n) {
    std::vector<int> p(static_cast<size_t>(n));
    std::iota(p.begin(), p.end(), 0);
    for (int i = n - 1; i > 0; --i) {
        int j = static_cast<int>(randint(0, i));
        std::swap(p[static_cast<size_t>(i)], p[static_cast<size_t>(j)]);
    }
    return p;
}

std::vector<long long> Rng::sample_without_replacement(long long n, long long k) {
    if (k > n) throw SnaqError("cannot sample more values than available without replacement");
    std::vector<long long> out;
    out.reserve(static_cast<size_t>(k));
    if (k * 4 > n || n < 64) {
        // partial Fisher-Yates
        std::vector<long long> pool(static_cast<size_t>(n));
        std::iota(pool.begin(), pool.end(), 0LL);
        for (long long i = 0; i < k; ++i) {
            long long j = randint(i, n - 1);
            std::swap(pool[static_cast<size_t>(i)], pool[static_cast<size_t>(j)]);
            out.push_back(pool[static_cast<size_t>(i)]);
        }
    } else {
        // Floyd's algorithm with a sorted vector check would be O(k^2); use hash set.
        std::vector<long long> chosen;
        chosen.reserve(static_cast<size_t>(k));
        // Simple rejection sampling with a sorted index for k << n.
        std::vector<uint64_t> seen;  // open addressing table
        size_t cap = 1;
        while (cap < static_cast<size_t>(k) * 2) cap <<= 1;
        seen.assign(cap, UINT64_MAX);
        auto insert = [&](uint64_t v) -> bool {
            size_t h = static_cast<size_t>((v * 0x9e3779b97f4a7c15ULL) >> 7) & (cap - 1);
            while (seen[h] != UINT64_MAX) {
                if (seen[h] == v) return false;
                h = (h + 1) & (cap - 1);
            }
            seen[h] = v;
            return true;
        };
        while (static_cast<long long>(out.size()) < k) {
            long long v = randint(0, n - 1);
            if (insert(static_cast<uint64_t>(v))) out.push_back(v);
        }
    }
    return out;
}

size_t Rng::weighted_index(const std::vector<double>& w) {
    double total = 0.0;
    for (double x : w) total += x;
    if (!(total > 0.0)) throw SnaqError("weighted sampling with non-positive total weight");
    double r = rand() * total;
    double acc = 0.0;
    for (size_t i = 0; i < w.size(); ++i) {
        acc += w[i];
        if (r < acc) return i;
    }
    return w.size() - 1;
}

long long binomial(long long n, long long k) {
    if (k < 0 || k > n) return 0;
    long long r = 1;
    for (long long i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

}  // namespace csnaq
