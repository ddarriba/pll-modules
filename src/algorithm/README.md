# libpll/algorithm

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

Module for miscellaneous algorithms involving several modules together
Exportable functions in this module start either with the prefix `pllmod_algo_`.

## Compilation instructions

The module can be compiled using the included Makefile:

`make`

Microsoft Windows compatibility was tested with a cross-compiler and seems to
work out-of-the-box using [MingW](http://www.mingw.org/).

## Code

The code is written in C.

    File              | Description
----------------------|----------------
**pll_modalgorithm.c** | High level algorithms.
**algo_callback.c**    | Internal callback functions
