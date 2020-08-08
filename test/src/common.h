/*
    Copyright (C) 2015 Diego Darriba, Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
#ifndef COMMON_H_
#define COMMON_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

/* parse attributes from the arguments */
unsigned int get_attributes(int argc, char **argv);
/* skip current test */
void skip_test();

/* callback function for traverse the utree */
int cb_full_traversal (pll_unode_t * node);
int cb_rfull_traversal (pll_rnode_t * node);

/* displays a tree */
void show_tree (pll_unode_t * tree, int SHOW_ASCII_TREE);
void show_rtree (pll_rnode_t * tree, int SHOW_ASCII_TREE);
/* print error and exit */
void fatal(const char * format, ...) __attribute__ ((noreturn));

#endif /* COMMON_H_ */
