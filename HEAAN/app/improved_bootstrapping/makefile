# This sample should be run after make libHEAAN.a file in lib folder

# All Target
all: clean build

# Tool invocations
build:
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ test_improved_bootstrapping.cpp -o test_improved_bootstrapping -L/usr/local/lib -I ../../../eigen ../../lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O3 -I/../../src -pthread
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	rm -rf test_improved_bootstrapping

run:
	./test_improved_bootstrapping
