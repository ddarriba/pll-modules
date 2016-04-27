# libpll/msa

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

Module for MSA management. Exportable functions in this module start with the prefix `pll_msa_`.

## Compilation instructions

The module can be compiled using the included Makefile:

`make`

Microsoft Windows compatibility was tested with a cross-compiler and seems to
work out-of-the-box using [MingW](http://www.mingw.org/).

## Code

The code is written in C.

    File              | Description
----------------------|----------------
**pll_msa.c**         | Functions for analyzing MSA embedded in `pll_partition_t` structure.
