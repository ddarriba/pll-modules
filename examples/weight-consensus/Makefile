# Copyright (C) 2016 Diego Darriba
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
# Contact: Tomas Flouri <Diego.Darriba@h-its.org>,
# Heidelberg Institute for Theoretical Studies,
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

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

CFLAGS = -g -O3 -Wall -Wsign-compare $(PROFILING) $(WARN)
LDFLAGS = -lpll_tree -lpll -lm

PROG = weight-consensus

default: all
all :
	cc -g -O3 -Wall -Wsign-compare -g -Wall -Wsign-compare  -c -o $(PROG).o $(PROG).c
	cc  -g -O3 -Wall -Wsign-compare -g -Wall -Wsign-compare  $(PROG).o -lpll_tree -lpll -lm -o $(PROG)


rf-distance : $(OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(PROG)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) $(PROG)
