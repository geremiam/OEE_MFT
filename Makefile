# Makefile
# Automates the compilation and building of the MFT solver

# Compiler name
CXX=g++-7
# Flags for including in path search
INC_FLAGS=-I${LAPACKE_INC}# -I${NETCDF_INC}
# Compiler flags (add -fopenmp to compilation and linking for OpenMP)
CXXFLAGS=-std=c++11 -O2
# Linker flags (add -fopenmp to compilation and linking for OpenMP)
LD_FLAGS=-L${LAPACKE_LIB}# -L${NETCDF_LIB}
# Flags for linking with libraries
LD_LIBS=-llapacke# -lnetcdf
HAM=ham1
# List of object files to be linked together
OBJECTS= $(HAM).o math_routines.o alloc_dealloc.o init_routines.o diag_routines.o driver.o
HEADERS= $(HAM).h math_routines.h alloc_dealloc.h init_routines.h diag_routines.h


## all: Default target; builds the final executable
.PHONY: all
all: driver

# Linking of the object files into the final executable. Depends on all .o files.
driver: $(OBJECTS)
	$(CXX) $(LD_FLAGS) $(LD_LIBS) $(OBJECTS) -o driver

# Compilation of the driver.cc source code into the driver.o object file.
driver.o: driver.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) -c driver.cc -o driver.o

# #######################################################################################
# MODULE HAM1

# Creation of the ham1.o object file, which depends on the module's header file and its 
# source file
ham1.o: ham1.cc ham1.h
	$(CXX) $(CXXFLAGS) -c ham1.cc -o ham1.o

## ut_ham1: Runs the testing suite for the module ham1
.PHONY: ut_ham1
ut_ham1: ham1_test # Runs the testing suite's executable
	./ham1_test

# Linking of the ham1_test.o object and ham1.o object files to create 
# the executable file.
ham1_test: ham1_test.o ham1.o
	$(CXX) ham1_test.o ham1.o -o ham1_test

# Creation of the ham1_test.o object file, which depends on its source file and 
# on the ham1.h header file.
ham1_test.o: ham1_test.cc ham1.h
	$(CXX) $(CXXFLAGS) -c ham1_test.cc -o ham1_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_ham1_clean
ut_ham1_clean:
	rm -f ham1_test.o ham1.o ham1_test

# #######################################################################################
# MODULE ALLOC_DEALLOC

# Creation of the alloc_dealloc.o object file, which depends on the module's header file 
# and its source file
alloc_dealloc.o: alloc_dealloc.cc alloc_dealloc.h
	$(CXX) $(CXXFLAGS) -c alloc_dealloc.cc -o alloc_dealloc.o

## ut_alloc_dealloc: Runs the testing suite for the module alloc_dealloc
.PHONY: ut_alloc_dealloc
ut_alloc_dealloc: alloc_dealloc_test # Runs the testing suite's executable
	./alloc_dealloc_test

# Linking of the alloc_dealloc_test.o object and alloc_dealloc.o object files to create 
# the executable file.
alloc_dealloc_test: alloc_dealloc_test.o alloc_dealloc.o
	$(CXX) alloc_dealloc_test.o alloc_dealloc.o -o alloc_dealloc_test

# Creation of the alloc_dealloc_test.o object file, which depends on its source file and 
# on the alloc_dealloc.h header file.
alloc_dealloc_test.o: alloc_dealloc_test.cc alloc_dealloc.h
	$(CXX) $(CXXFLAGS) -c alloc_dealloc_test.cc -o alloc_dealloc_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_alloc_dealloc_clean
ut_alloc_dealloc_clean:
	rm -f alloc_dealloc_test.o alloc_dealloc.o alloc_dealloc_test

# #######################################################################################
# MODULE INIT_ROUTINES

# Creation of the init_routines.o object file, which depends on the module's header file 
# and its source file
init_routines.o: init_routines.cc init_routines.h
	$(CXX) $(CXXFLAGS) -c init_routines.cc -o init_routines.o

## ut_init_routines: Runs the testing suite for the module init_routines
.PHONY: ut_init_routines
ut_init_routines: init_routines_test # Runs the testing suite's executable
	./init_routines_test

# Linking of the init_routines_test.o object and init_routines.o object files to create 
# the executable file.
init_routines_test: init_routines_test.o init_routines.o
	$(CXX) init_routines_test.o init_routines.o -o init_routines_test

# Creation of the init_routines_test.o object file, which depends on its source file and 
# on the init_routines.h header file.
init_routines_test.o: init_routines_test.cc init_routines.h
	$(CXX) $(CXXFLAGS) -c init_routines_test.cc -o init_routines_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_init_routines_clean
ut_init_routines_clean:
	rm -f init_routines_test.o init_routines.o init_routines_test


# #######################################################################################
# MODULE MATH_ROUTINES

# Creation of the math_routines.o object file, which depends on the module's header file 
# and its source file
math_routines.o: math_routines.cc math_routines.h
	$(CXX) $(CXXFLAGS) -c math_routines.cc -o math_routines.o

## ut_math_routines: Runs the testing suite for the module math_routines
.PHONY: ut_math_routines
ut_math_routines: math_routines_test # Runs the testing suite's executable
	./math_routines_test

# Linking of the math_routines_test.o object and math_routines.o object files to create 
# the executable file.
math_routines_test: math_routines_test.o math_routines.o
	$(CXX) math_routines_test.o math_routines.o -o math_routines_test

# Creation of the math_routines_test.o object file, which depends on its source file and 
# on the math_routines.h header file.
math_routines_test.o: math_routines_test.cc math_routines.h
	$(CXX) $(CXXFLAGS) -c math_routines_test.cc -o math_routines_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_math_routines_clean
ut_math_routines_clean:
	rm -f math_routines_test.o math_routines.o math_routines_test

# #######################################################################################
# MODULE DIAG_ROUTINES

# diag_routines.o object file depends on header file and source file
# Must include LAPACKE_INC in the path search
diag_routines.o: diag_routines.cc diag_routines.h
	$(CXX) $(CXXFLAGS) -I${LAPACKE_INC} -c diag_routines.cc -o diag_routines.o

## ut_diag_routines: Runs the testing suite for the module diag_routines
.PHONY: ut_diag_routines
ut_diag_routines: diag_routines_test # Runs the testing suite's executable
	./diag_routines_test

# Testing suite executable depends on diag_routines_test.o diag_routines.o
# Must include LAPACKE_LIB in path search and lapacke library
diag_routines_test: diag_routines_test.o diag_routines.o
	$(CXX) -L${LAPACKE_LIB} -llapacke diag_routines_test.o diag_routines.o \
	-o diag_routines_test

# diag_routines_test.o object file depends on source file and diag_routines.h header file
diag_routines_test.o: diag_routines_test.cc diag_routines.h
	$(CXX) $(CXXFLAGS) -c diag_routines_test.cc -o diag_routines_test.o

# Clean target for this unit test
.PHONY: ut_diag_routines_clean
ut_diag_routines_clean:
	rm -f diag_routines_test.o diag_routines.o diag_routines_test

# #######################################################################################

## clean: remove object files and executable files (other than unit test files)
.PHONY: clean
clean: 
	rm -f $(OBJECTS) driver

## ut_clean: Run clean rules for all unit tests
.PHONY: ut_clean
ut_clean: ut_ham1_clean \
          ut_math_routines_clean \
          ut_alloc_dealloc_clean \
          ut_init_routines_clean \
          ut_diag_routines_clean

## help: Show targets and their descriptions
.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<