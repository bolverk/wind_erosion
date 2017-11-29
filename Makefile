LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wno-long-long -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations -Wconversion -O3 -march=native
CXX = g++
build/rich.exe: build/rich.o
	$(CXX)  $< $(RICH_ROOT)/library_production/librich.a -L $(HDF5_LIB_PATH) -lhdf5 -lhdf5_cpp -lz -o $@
build/rich.o: source/rich.cpp build
	$(CXX) -c $(LINT_FLAGS) $< -o $@ -I $(RICH_ROOT)
build:
	mkdir build