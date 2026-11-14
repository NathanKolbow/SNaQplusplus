# One-command build for macOS/Linux: `make` produces ./bin/snaq.
# (This wraps CMake; on Windows, or without `make`, run: cmake -P build.cmake)

BUILD_DIR ?= build
PREFIX    ?= $(HOME)/.local
JOBS      ?= 4

.PHONY: all snaq install test test-julia clean tools bench bench-quick bench-fit bench-julia bench-report

all: snaq

snaq:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release
	cmake --build $(BUILD_DIR) --config Release --parallel $(JOBS)
	@mkdir -p bin
	@cp $(BUILD_DIR)/snaq bin/snaq
	@echo ""
	@echo "Done: the snaq program is at bin/snaq"
	@echo "Try: ./bin/snaq --help"

install: snaq
	cmake --install $(BUILD_DIR) --prefix $(PREFIX)
	@echo "Installed snaq to $(PREFIX)/bin"

test:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DCSNAQ_BUILD_TESTS=ON
	cmake --build $(BUILD_DIR) --config Release --parallel $(JOBS)
	ctest --test-dir $(BUILD_DIR) --output-on-failure

# The unit tests, plus the tests that check snaq's output files with SNaQ.jl
# (needs julia and the environment of tools/julia: see tools/julia/setup.jl).
test-julia:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DCSNAQ_BUILD_TESTS=ON -DCSNAQ_JULIA_TESTS=ON
	cmake --build $(BUILD_DIR) --config Release --parallel $(JOBS)
	ctest --test-dir $(BUILD_DIR) --output-on-failure

clean:
	rm -rf $(BUILD_DIR) bin

# ---------------------------------------------------------------- benchmarks
# See benchmarks/README.md. LABEL names the version being benchmarked,
# THREADS sets the number of threads, ARGS passes extra flags to snaq.
LABEL   ?=
THREADS ?= 1
ARGS    ?=
SUITE   ?= standard
BENCH_OPTS = --threads $(THREADS) $(if $(LABEL),--label "$(LABEL)") $(if $(ARGS),--extra-args "$(ARGS)")

tools:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DCSNAQ_BUILD_TOOLS=ON
	cmake --build $(BUILD_DIR) --config Release --parallel $(JOBS)

bench: tools
	python3 benchmarks/bench.py run --suite $(SUITE) $(BENCH_OPTS)

bench-quick: tools
	python3 benchmarks/bench.py run --suite quick $(BENCH_OPTS)

bench-fit: tools
	python3 benchmarks/bench.py run --suite fit $(BENCH_OPTS)

bench-julia: tools
	python3 benchmarks/bench.py run-julia --suite $(SUITE) --threads $(THREADS)

bench-report:
	python3 benchmarks/bench.py report
