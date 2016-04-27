# libpll/optimize

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

Module for model parameters optimization. Exportable functions in this module start either with the prefix `pll_optimize_` for high level optimization algorithms, or `pll_minimize_` for low level optimization functions.

## Compilation instructions

The module can be compiled using the included Makefile:

`make`

Microsoft Windows compatibility was tested with a cross-compiler and seems to
work out-of-the-box using [MingW](http://www.mingw.org/).

## Code

The code is written in C.

    File              | Description
----------------------|----------------
**pll_optimize.c**    | High level optimization algorithms.
**opt_algorithms.c**  | Low level optimization algorithms.
**lbfgsb/**           | Core implementation of L-BFGS-B algorithm.
