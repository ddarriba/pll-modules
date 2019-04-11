# Copyright (C) 2019 Alexey Kozlov
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
# Heidelberg Institute for Theoretical Studies,
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

INSTALLPATH=../../install

# Profiling
OS = $(shell uname -s)
PROFILING=-g
CC ?= gcc

# Compiler warnings
ifeq ($(CC),clang)
	ADD_WARN=-Weverything -Wno-padded -Wno-missing-noreturn \
	         -Wno-format-nonliteral -Wno-unused-parameter
endif
WARN=-Wall -Wsign-compare $(ADD_WARN)

CFLAGS = -g -O3 -static -Wall -Wsign-compare $(PROFILING) $(WARN)
LDFLAGS = -lpll_algorithm -lpll_optimize -lpll_tree -lpll_util -lpll -lm

# temp workaround
CFLAGS += -I$(INSTALLPATH)/include -L$(INSTALLPATH)/lib

OBJS = genotype.o
PROG = genotype

default: all
all : $(PROG)

genotype : $(OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(PROG)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) $(PROG)
