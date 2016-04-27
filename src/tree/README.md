# libpll/tree

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

Module for tree management. Exportable functions in this module start either with the prefix `pll_tree_` for general tree management functions, `pll_utree_` for unrooted tree operations, or `pll_rtree_` for rooted tree operations.

## Compilation instructions

The module can be compiled using the included Makefile:

`make`

Microsoft Windows compatibility was tested with a cross-compiler and seems to
work out-of-the-box using [MingW](http://www.mingw.org/).

## Code

The code is written in C.

    File               | Description
-----------------------|----------------
**pll_tree.c**         | Functions for performing complete operations.
**utree_operations.c** | Operations on unrooted trees.
