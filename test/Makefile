#    Copyright (C) 2015 Diego Darriba
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
#    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
#    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

CC = gcc
CFLAGS = -g -O3 -Wall -std=c99
CLIBS = -lpll -lm -lpll_optimize -lpll_tree -lpll_binary

ifdef LIBPLL_INC
  CFLAGS += -I$(LIBPLL_INC)
  ifdef LIBPLL_LIB
    CFLAGS += -L$(LIBPLL_LIB)
  else
    CFLAGS += -L$(LIBPLL_INC)
  endif
else
  CFLAGS += -I../install/include/libpll -L../install/lib
endif

MODULES = binary optimize tree

CFILES = src/binary/binary-sequential.c \
         src/binary/binary-random.c \
         src/binary/binary-skeleton.c \
         src/optimize/blopt-minimal.c \
         src/optimize/blopt-5states.c \
         src/tree/random-tree.c \
         src/tree/parsimony-tree.c \
         src/tree/treemove-nni.c \
         src/tree/treemove-spr.c \
         src/tree/rtreemove-spr.c \
         src/tree/treemove-tbr.c \
         src/tree/serialize.c \
		 src/tree/split-reconstruct.c

OBJFILES = $(patsubst src/%.c, obj/%, $(CFILES))

OBJCOMMON = src/common.o

DATADIR = testdata
RESULTDIR = result

DEPFILES = worms16s.fas \
           small.fas small.tree small.rooted.tree \
           246x4465.fas 246x4465.tree \
           ribosomal_l5_pf00673.fas \
           medium.fas medium.tree medium.rooted.tree
REQFILES = $(patsubst %, $(DATADIR)/%, $(DEPFILES))
ASSETS = https://raw.githubusercontent.com/xflouris/assets/master/libpll

all: $(OBJCOMMON) $(OBJFILES) $(REQFILES)
	@mkdir -p $(RESULTDIR) $(foreach m,$(MODULES), $(RESULTDIR)/$(m))

$(DATADIR)/%:
	@mkdir -p "$(@D)"
	wget -O $@ $(ASSETS)/$@

obj/%: src/%.c $(DEPS)
	@mkdir -p "$(@D)"
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(OBJCOMMON) $(CLIBS) $(LDFLAGS) -DDATADIR=\"$(DATADIR)\"

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $< $(CLIBS) $(LDFLAGS) -DDATADIR=\"$(DATADIR)\"
	
clean:
	rm -rf obj result $(OBJCOMMON)
