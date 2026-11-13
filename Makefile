# One-command build for macOS/Linux: `make` produces ./bin/snaq.
# (This wraps CMake; on Windows, or without `make`, run: cmake -P build.cmake)

BUILD_DIR ?= build
PREFIX    ?= $(HOME)/.local
JOBS      ?= 4

.PHONY: all snaq install test clean

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

clean:
	rm -rf $(BUILD_DIR) bin
