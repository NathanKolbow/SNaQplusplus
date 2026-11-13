// Minimal fork-join parallel loop over index ranges using std::thread.
#pragma once

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace csnaq {

// Whether a loop over n items is worth splitting across nthreads
// (SNaQ's `useparallelloop`).
inline bool useparallelloop(size_t n, int nthreads) { return nthreads > 1 && n >= 4 * static_cast<size_t>(nthreads); }

// Calls f(lo, hi, tid) on contiguous chunks of [0, n). Exceptions thrown by
// workers are rethrown in the calling thread.
template <class F>
void parallel_for(size_t n, int nthreads, F&& f) {
    if (!useparallelloop(n, nthreads)) {
        f(size_t{0}, n, 0);
        return;
    }
    const size_t nt = static_cast<size_t>(nthreads);
    const size_t chunk = (n + nt - 1) / nt;
    std::vector<std::thread> th;
    std::exception_ptr err = nullptr;
    std::mutex m;
    for (size_t k = 1; k < nt; ++k) {
        size_t lo = k * chunk, hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        th.emplace_back([&, lo, hi, k] {
            try {
                f(lo, hi, static_cast<int>(k));
            } catch (...) {
                std::lock_guard<std::mutex> g(m);
                if (!err) err = std::current_exception();
            }
        });
    }
    try {
        f(size_t{0}, std::min(n, chunk), 0);
    } catch (...) {
        std::lock_guard<std::mutex> g(m);
        if (!err) err = std::current_exception();
    }
    for (auto& t : th) t.join();
    if (err) std::rethrow_exception(err);
}

}  // namespace csnaq
