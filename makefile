# Define the target to build:
TARGET = PALET

# Define paths *PLEASE UPDATE WITH YOUR DIRECTORIES*
HOME = /home/grads/k/k.kitzmiller

GMP = $(HOME)/spack/opt/spack/linux-ubuntu18.04-x86_64_v4/gcc-8.4.0/gmp-6.2.1-vsw45xa4h6x5usr5r7r4uicuvbqbv7vi
MPFR = $(HOME)/spack/opt/spack/linux-ubuntu18.04-x86_64_v4/gcc-8.4.0/mpfr-4.1.0-gi2lorkelalblnnn3dqfexhhnk7zvyax
GIVARO = $(HOME)/givaro-4.2.0
OBLAS = $(HOME)/OpenBLAS-0.3.21
FFLAS = $(HOME)/fflas-ffpack-2.5.0

# Define path to header files
HEADERS = -I$(GMP)/include -I$(MPFR)/include -I$(GIVARO)/lib/include -I$(OBLAS)/lib/include -I$(FFLAS)/lib/include

# Define compile time library paths
CMP_LIBPATHS = -lm -L$(GMP)/lib -lgmp -lgmpxx -L$(MPFR)/lib -lmpfr -L$(GIVARO)/lib/lib -lgivaro -L$(OBLAS)/lib/lib -lopenblas -L$(FFLAS)/lib/lib -lgomp

# Define run time library paths
RPATH = -Wl,-rpath -Wl,
RUN_LIBPATHS = $(RPATH)$(GMP)/lib $(RPATH)$(MPFR)/lib $(RPATH)$(GIVARO)/lib/lib $(RPATH)$(OBLAS)/lib/lib $(RPATH)$(FFLAS)/lib/lib

# Define the compiler
CMP = g++
EXT = .cpp

# Define compiler flags:
FLAGS = -Wall -std=c++11 -fopenmp
FLAGS_DBG = $(FLAGS) -g3 -Og
FLAGS_RUN = $(FLAGS) -O2 -DNDEBUG

all: run

run:	$(TARGET)$(EXT)
	$(CMP) $(FLAGS_RUN) $(HEADERS) -o $(TARGET) $(TARGET)$(EXT) $(CMP_LIBPATHS) $(RUN_LIBPATHS)

clean:
	rm $(TARGET)
