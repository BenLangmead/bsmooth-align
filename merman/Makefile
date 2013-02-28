#!gmake

#
# Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
#
# This file is part of Merman.
#
# Merman is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Merman is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Merman.  If not, see <http://www.gnu.org/licenses/>.
#

USE_PTHREADS = 1

# Detect Cygwin or MinGW
WINDOWS = 0
ifneq (,$(findstring CYGWIN,$(shell uname)))
WINDOWS = 1
else
ifneq (,$(findstring MINGW,$(shell uname)))
WINDOWS = 1
endif
endif

SRC=mer_index.cpp \
    alphabet.cpp \
    check.cpp \
    merman.cpp \
    qual.cpp \
    read.cpp \
    seeds.cpp \
    edit.cpp \
    hit_set.cpp \
    refmap.cpp \
    annot.cpp \
    threading.cpp \
    color_dec.cpp \
    color.cpp \
    output.cpp \
    align_sw.cpp \
    align_naive.cpp \
    ref.cpp \
	orient.cpp \
	merman_main.cpp

HEAD=mer_ent.h \
    mer_index.h \
    alphabet.h \
    check.h \
    tokenize.h \
    filebuf.h \
    timer.h \
    qual.h \
    read.h \
    seeds.h \
    edit.h \
    hit_set.h \
    refmap.h \
    annot.h \
    threading.h \
    parsort.h \
    color_dec.h \
    color.h \
    output.h \
    ds.h \
    toa.h \
    random_source.h \
    ref.h \
    align_naive.h \
    align_sw.h \
    sstring.h \
	merman.h \
	orient.h

BITS_FLAG =
ifeq (32,$(BITS))
BITS_FLAG = -m32
endif
ifeq (64,$(BITS))
BITS_FLAG = -m64
endif

PTHREAD_PKG =
PTHREAD_LIB =
PTHREAD_DEF =
ifeq (1,$(USE_PTHREADS))
PTHREAD_DEF = -DUSE_PTHREADS
ifeq (1,$(WINDOWS))
# pthreads for windows forces us to be specific about the library
PTHREAD_LIB = -lpthreadGC2
PTHREAD_PKG = pthreadGC2.dll
else
# There's also -pthread, but that only seems to work on Linux
PTHREAD_LIB = -lpthread
endif
endif

VERSION = $(shell cat VERSION)
DEFS = $(PTHREAD_DEF) -DMERMAN_VERSION="\"$(VERSION)\""
LIBS = $(PTHREAD_LIB)

all: merman merman-debug

merman: $(SRC) $(HEAD)
	g++ $(BITS_FLAG) $(DEFS) $(LIBS) -Wall -g -DNDEBUG -O3 -o $@ $(SRC)

merman-debug: $(SRC) $(HEAD)
	g++ $(BITS_FLAG) $(DEFS) $(LIBS) -Wall -g -O0 -o $@ $(SRC)

merman-prof: $(SRC) $(HEAD)
	g++ $(BITS_FLAG) $(DEFS) $(LIBS) -g -pg -Wall -DNDEBUG -O3 -o $@ $(SRC)

COLOR_DEC_SRC=color_dec.cpp \
              alphabet.cpp \
              color.cpp

color-dec: $(COLOR_DEC_SRC)
	g++ $(BITS_FLAG) $(DEFS) $(LIBS) -DMAIN_COLOR_DEC -g -pg -Wall -DNDEBUG -O0 $(COLOR_DEC_SRC) -o $@

.PHONY: clean
clean:
	rm -f merman merman-debug merman-prof
