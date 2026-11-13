// Small shared utilities: errors, Julia-compatible number formatting, RNG,
// an inline small vector, and timing helpers.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace csnaq {

// Generic error (mirrors Julia's `error(...)`).
struct SnaqError : public std::runtime_error {
    explicit SnaqError(const std::string& msg) : std::runtime_error(msg) {}
};

// Mirrors PhyloNetworks.RootMismatch.
struct RootMismatch : public SnaqError {
    explicit RootMismatch(const std::string& msg) : SnaqError(msg) {}
};

// Mirrors Julia's ArgumentError (raised for invalid user inputs).
struct ArgumentError : public SnaqError {
    explicit ArgumentError(const std::string& msg) : SnaqError(msg) {}
};

// Print a warning to stderr, formatted like Julia's `@warn`.
void warn(const std::string& msg);

// Julia's `string(x::Float64)`: shortest round-trip representation, with
// Julia's choice between plain ("0.0001", "123456.0") and exponent
// ("1.0e-5", "1.234567e6") forms.
std::string julia_float(double x);

// Julia's `round(x, digits=d)` for Float64.
double julia_round(double x, int digits);

// Julia's `Int(round(x))` (round half to even).
long long julia_round_int(double x);

// Elapsed-time description used in SNaQ log files (`timeelapsed`).
std::string time_elapsed(double seconds);

// "HH:MM:SS yyyy-mm-dd" (`currenttime` in SNaQ).
std::string current_time();

// "yyyy-mm-dd HH:MM:SS" (timestamps of `logtext`).
std::string current_timestamp();

// Wall clock in seconds.
double wall_time();

// Peak resident set size of this process, in bytes (0 if unavailable).
size_t peak_rss_bytes();

// -------------------------------------------------------------------------
// RNG: xoshiro256++ seeded with splitmix64. SNaQ.jl results are not meant
// to be reproduced bit-for-bit, so this does not mimic Julia's streams.
// -------------------------------------------------------------------------
class Rng {
public:
    explicit Rng(uint64_t seed = 0) { this->seed(seed); }
    void seed(uint64_t seed);
    uint64_t next_u64();
    // Uniform in [0, 1).
    double rand();
    // Uniform integer in [lo, hi] (inclusive).
    long long randint(long long lo, long long hi);
    // Uniform index in [0, n).
    size_t index(size_t n) { return static_cast<size_t>(randint(0, static_cast<long long>(n) - 1)); }
    // Random Int64, like Julia's `rand(Int)`.
    long long rand_int64() { return static_cast<long long>(next_u64()); }
    // Random permutation of 0..n-1.
    std::vector<int> randperm(int n);
    // k distinct values sampled from 0..n-1, in sampling order.
    std::vector<long long> sample_without_replacement(long long n, long long k);
    // Index sampled with probability proportional to `weights`.
    size_t weighted_index(const std::vector<double>& weights);

    template <class T>
    const T& sample(const std::vector<T>& v) { return v[index(v.size())]; }

private:
    uint64_t s_[4];
};

// -------------------------------------------------------------------------
// SmallVec: vector with N inline slots and heap overflow (node edge lists).
// -------------------------------------------------------------------------
template <class T, int N>
class SmallVec {
public:
    SmallVec() = default;
    SmallVec(std::initializer_list<T> il) {
        for (const T& x : il) push_back(x);
    }
    SmallVec(const SmallVec& o) { copy_from(o); }
    SmallVec& operator=(const SmallVec& o) {
        if (this != &o) {
            clear();
            copy_from(o);
        }
        return *this;
    }
    SmallVec(SmallVec&& o) noexcept { move_from(std::move(o)); }
    SmallVec& operator=(SmallVec&& o) noexcept {
        if (this != &o) {
            free_heap();
            move_from(std::move(o));
        }
        return *this;
    }
    ~SmallVec() { free_heap(); }

    size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }
    T* data() { return heap_ ? heap_ : inline_; }
    const T* data() const { return heap_ ? heap_ : inline_; }
    T& operator[](size_t i) { return data()[i]; }
    const T& operator[](size_t i) const { return data()[i]; }
    T* begin() { return data(); }
    T* end() { return data() + size_; }
    const T* begin() const { return data(); }
    const T* end() const { return data() + size_; }
    T& back() { return data()[size_ - 1]; }

    void push_back(const T& x) {
        if (size_ == cap_) grow();
        data()[size_++] = x;
    }
    void clear() { size_ = 0; }
    // Remove element at position i, shifting the rest (Julia deleteat!).
    void erase_at(size_t i) {
        T* d = data();
        for (size_t k = i; k + 1 < size_; ++k) d[k] = d[k + 1];
        --size_;
    }
    // Index of first element equal to x, or -1.
    int find(const T& x) const {
        const T* d = data();
        for (size_t k = 0; k < size_; ++k)
            if (d[k] == x) return static_cast<int>(k);
        return -1;
    }
    bool contains(const T& x) const { return find(x) >= 0; }
    void replace(const T& from, const T& to) {
        T* d = data();
        for (size_t k = 0; k < size_; ++k)
            if (d[k] == from) d[k] = to;
    }

private:
    T inline_[N];
    T* heap_ = nullptr;
    size_t size_ = 0;
    size_t cap_ = N;

    void grow() {
        size_t nc = cap_ * 2;
        T* nh = new T[nc];
        for (size_t k = 0; k < size_; ++k) nh[k] = data()[k];
        delete[] heap_;
        heap_ = nh;
        cap_ = nc;
    }
    void free_heap() {
        delete[] heap_;
        heap_ = nullptr;
        cap_ = N;
        size_ = 0;
    }
    void copy_from(const SmallVec& o) {
        for (size_t k = 0; k < o.size_; ++k) push_back(o.data()[k]);
    }
    void move_from(SmallVec&& o) {
        if (o.heap_) {
            heap_ = o.heap_;
            cap_ = o.cap_;
            size_ = o.size_;
            o.heap_ = nullptr;
            o.cap_ = N;
            o.size_ = 0;
        } else {
            heap_ = nullptr;
            cap_ = N;
            size_ = 0;
            for (size_t k = 0; k < o.size_; ++k) push_back(o.inline_[k]);
            o.size_ = 0;
        }
    }
};

// Index of the lowest set bit of x (x != 0).
inline int ctz64(uint64_t x) {
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned long i;
    _BitScanForward64(&i, x);
    return static_cast<int>(i);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int i = 0;
    while (!(x & 1)) {
        x >>= 1;
        ++i;
    }
    return i;
#endif
}

// Number of 4-taxon subsets of n taxa.
inline long long nchoose4(long long n) { return n < 4 ? 0 : n * (n - 1) * (n - 2) * (n - 3) / 24; }

// Julia's `binomial(n, k)` for small k.
long long binomial(long long n, long long k);

}  // namespace csnaq
