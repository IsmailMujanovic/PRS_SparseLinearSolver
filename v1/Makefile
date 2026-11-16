CXX ?= c++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -Iinclude

SUPERLU_PREFIX := $(shell brew --prefix superlu 2>/dev/null)
OPENBLAS_PREFIX := $(shell brew --prefix openblas 2>/dev/null)

ifneq ($(strip $(SUPERLU_PREFIX)),)
    ifeq ($(wildcard $(SUPERLU_PREFIX)/include/slu_zdefs.h),)
        SUPERLU_PREFIX :=
    endif
endif

ifneq ($(strip $(OPENBLAS_PREFIX)),)
    ifeq ($(wildcard $(OPENBLAS_PREFIX)/lib/libopenblas.*),)
        OPENBLAS_PREFIX :=
    endif
endif

SUPERLU_CFLAGS ?=
SUPERLU_LIBS ?=
BLAS_LIBS ?= -lopenblas

ifeq ($(strip $(SUPERLU_CFLAGS)),)
    SUPERLU_CFLAGS := $(shell pkg-config --cflags superlu 2>/dev/null)
endif

ifeq ($(strip $(SUPERLU_LIBS)),)
    SUPERLU_LIBS := $(shell pkg-config --libs superlu 2>/dev/null)
endif

ifeq ($(strip $(SUPERLU_CFLAGS)),)
    ifneq ($(strip $(SUPERLU_PREFIX)),)
        ifneq ($(wildcard $(SUPERLU_PREFIX)/include/slu_zdefs.h),)
            SUPERLU_CFLAGS := -I$(SUPERLU_PREFIX)/include
        endif
    endif
endif

ifeq ($(strip $(SUPERLU_CFLAGS)),)
    ifneq ($(wildcard /usr/include/superlu/slu_zdefs.h),)
        SUPERLU_CFLAGS := -I/usr/include/superlu
    else ifneq ($(wildcard /usr/local/include/superlu/slu_zdefs.h),)
        SUPERLU_CFLAGS := -I/usr/local/include/superlu
    else ifneq ($(wildcard /usr/include/slu_zdefs.h),)
        SUPERLU_CFLAGS := -I/usr/include
    else ifneq ($(wildcard /usr/local/include/slu_zdefs.h),)
        SUPERLU_CFLAGS := -I/usr/local/include
    endif
endif

ifeq ($(strip $(SUPERLU_LIBS)),)
    ifneq ($(strip $(SUPERLU_PREFIX)),)
        SUPERLU_LIBS := -L$(SUPERLU_PREFIX)/lib -lsuperlu
    else
        SUPERLU_LIBS := -lsuperlu
    endif
endif

ifeq ($(strip $(OPENBLAS_PREFIX)),)
    OPENBLAS_PREFIX := /usr
endif

BENCH_CPPFLAGS := $(SUPERLU_CFLAGS)
BENCH_LDFLAGS := $(SUPERLU_LIBS) -L$(OPENBLAS_PREFIX)/lib
BENCH_LDLIBS := $(BLAS_LIBS)

BENCH_RPATH :=
ifneq ($(strip $(SUPERLU_PREFIX)),)
    ifneq ($(SUPERLU_PREFIX),/usr)
        BENCH_RPATH += -Wl,-rpath,$(SUPERLU_PREFIX)/lib
    endif
endif
ifneq ($(OPENBLAS_PREFIX),/usr)
    BENCH_RPATH += -Wl,-rpath,$(OPENBLAS_PREFIX)/lib
endif

SRC := src/matrix.cpp src/md.cpp src/permute.cpp src/symbolic.cpp src/numeric.cpp
PROGRAM_SRC := src/main.cpp src/validate.cpp $(SRC)
BENCH_SRC := benchmark/bench_superlu.cpp benchmark/superlu_adapter.cpp $(SRC)
HEADERS := $(wildcard include/*.hpp) benchmark/superlu_adapter.hpp

.PHONY: all program bench run clean

all: program

program: build/program

bench: build/bench

build/program: $(PROGRAM_SRC) $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(PROGRAM_SRC) -o $@

build/bench: $(BENCH_SRC) $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(BENCH_CPPFLAGS) $(BENCH_SRC) $(BENCH_LDFLAGS) $(BENCH_RPATH) $(BENCH_LDLIBS) -o $@

run: build/bench
	./benchmark/run.sh

clean:
	rm -rf build
