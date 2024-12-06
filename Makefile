CXX = clang++
CXXFLAGS = -g -O0 -std=c++20
INCLUDES = -I/opt/homebrew/include/eigen3 -I/opt/homebrew/include 
LDFLAGS = -L/opt/homebrew/lib -lfmt

bk: bk.cc spin_op.o xtensor_impl.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ $(LDFLAGS) -o $@

xtensor_impl.o: xtensor_impl.cpp
	$(CXX) -c -std=c++20 $(INCLUDES) $< -o $@

spin_op.o: spin_op.cpp
	$(CXX) -c -std=c++20 $(INCLUDES) $< -L/opt/homebrew/lib/fmt -lfmt -o $@

.PHONY: clean
clean:
	rm -f bk *.o
