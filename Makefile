CXX = clang++
CXXFLAGS = -g -O0 -std=c++20
INCLUDES = -I/opt/homebrew/include/eigen3 -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lfmt
TEST_LDFLAGS = $(LDFLAGS) -lgtest -lgtest_main -pthread

bravyi_kitaev : bravyi_kitaev.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

test_bravyi_kitaev: test_bravyi_kitaev.cpp bravyi_kitaev.o spin_op.o xtensor_impl.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ $(TEST_LDFLAGS) -o $@

xtensor_impl.o: xtensor_impl.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

spin_op.o: spin_op.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -L/opt/homebrew/lib/fmt -lfmt -o $@

.PHONY: clean test

test: test_bravyi_kitaev
	./test_bravyi_kitaev

clean:
	rm -f test_bravyi_kitaev *.o
