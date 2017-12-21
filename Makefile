# Makefile
# Automates the compilation and building of the MFT solver

# Compiler name
CXX=clang++#g++-7
# Flags for including in path search
INC_FLAGS=-I${LAPACKE_INC} -I${NETCDF_INC}
# Compiler flags (add -fopenmp to compilation and linking for OpenMP)
CXXFLAGS=-std=c++11 -O2
# Linker flags (add -fopenmp to compilation and linking for OpenMP)
LD_FLAGS=-L${LAPACKE_LIB} -L${NETCDF_LIB}
# Flags for linking with libraries
LD_LIBS=-llapacke -lnetcdf
# List of object files and header files belonging to modules
OBJECTS=alloc_dealloc.o init_routines.o math_routines.o diag_routines.o kspace.o nc_IO.o
HEADERS=alloc_dealloc.h init_routines.h math_routines.h diag_routines.h kspace.h nc_IO.h


## all: Default target; empty
.PHONY: all
all: help

# #######################################################################################
# DRIVER_HAM1

## driver_ham1: Builds the final executable for driver_ham1
# Linking of the object files into the final executable. Depends on all .o files.
driver_ham1: driver_ham1.o $(OBJECTS)
	$(CXX) $(LD_FLAGS) $(LD_LIBS) driver_ham1.o $(OBJECTS) -o driver_ham1

# Creation of the driver_ham1.o object file, which depends on the module's header file 
# and the other headers
driver_ham1.o: driver_ham1.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) -c driver_ham1.cc -o driver_ham1.o

# Deletion of the object file and executable file for this driver
.PHONY: driver_ham1_clean
driver_ham1_clean:
	rm -f driver_ham1.o driver_ham1

# #######################################################################################
# DRIVER_HAM2

## driver_ham2: Builds the final executable for driver_ham2
# Linking of the object files into the final executable. Depends on all .o files.
driver_ham2: driver_ham2.o $(OBJECTS)
	$(CXX) $(LD_FLAGS) $(LD_LIBS) driver_ham2.o $(OBJECTS) -o driver_ham2

# Creation of the driver_ham2.o object file, which depends on the header files
driver_ham2.o: driver_ham2.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) -c driver_ham2.cc -o driver_ham2.o

# Deletion of the object file and executable file for this driver
.PHONY: driver_ham2_clean
driver_ham2_clean:
	rm -f driver_ham2.o driver_ham2

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
# MODULE KSPACE

# kspace.o object file depends on header file, source file, and all included header files
kspace.o: kspace.cc kspace.h alloc_dealloc.h init_routines.h
	$(CXX) $(CXXFLAGS) -c kspace.cc -o kspace.o

## ut_kspace: Runs the testing suite for the module kspace
.PHONY: ut_kspace
ut_kspace: kspace_test # Runs the testing suite's executable
	./kspace_test

# Testing suite executable depends on kspace_test.o, kspace.o, and the other modules used 
# in the source code.
kspace_test: kspace_test.o kspace.o alloc_dealloc.o init_routines.o
	$(CXX) kspace_test.o kspace.o alloc_dealloc.o init_routines.o -o kspace_test

# kspace_test.o object file depends on source file and kspace.h header
kspace_test.o: kspace_test.cc kspace.h
	$(CXX) $(CXXFLAGS) -c kspace_test.cc -o kspace_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_kspace_clean
ut_kspace_clean:
	rm -f alloc_dealloc.o init_routines.o kspace_test.o kspace.o kspace_test
# #######################################################################################
# MODULE NC_IO

# nc_IO.o object file depends on header file, source file, and all included header files
nc_IO.o: nc_IO.cc nc_IO.h
	$(CXX) $(CXXFLAGS) -I${NETCDF_INC} -c nc_IO.cc -o nc_IO.o

## ut_nc_IO: Runs the testing suite for the module nc_IO
.PHONY: ut_nc_IO
ut_nc_IO: nc_IO_test # Runs the testing suite's executable
	./nc_IO_test

# Testing suite executable depends on nc_IO_test.o, Z.o, and the other modules used 
# in the source code.
nc_IO_test: nc_IO_test.o nc_IO.o
	$(CXX) -L${NETCDF_LIB} -lnetcdf nc_IO_test.o nc_IO.o -o nc_IO_test

# nc_IO_test.o object file depends on source file and nc_IO.h header
nc_IO_test.o: nc_IO_test.cc nc_IO.h
	$(CXX) $(CXXFLAGS) -c nc_IO_test.cc -o nc_IO_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_nc_IO_clean
ut_nc_IO_clean:
	rm -f nc_IO_test.o nc_IO.o nc_IO_test
# #######################################################################################

## clean: Removes module object files as well as driver object files and executables
.PHONY: clean
clean: driver_ham1_clean driver_ham2_clean
	rm -f $(OBJECTS)

## ut_clean: Runs clean rules for all unit tests
.PHONY: ut_clean
ut_clean: ut_alloc_dealloc_clean \
          ut_init_routines_clean \
          ut_math_routines_clean \
          ut_diag_routines_clean \
          ut_kspace_clean \
          ut_nc_IO_clean

## help: Shows targets and their descriptions
.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<