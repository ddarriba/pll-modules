# libpll/binary

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

Module containing functions for binary I/O. Exportable functions in this module start with the prefix `pllmod_binary_`.

## Compilation instructions

Read the compilation instructions at the main README file

## Code

The code is written in C.

|    File                   | Description                     |
|---------------------------|---------------------------------|
|**pll_binary.c**           | Interface functions.            |
|**binary_io_operations.c** | Operations with binary files.   |

## Type definitions

* struct `pllmod_subst_model_t`
* struct `pllmod_mixture_model_t`

## Flags

* `PLLMOD_BIN_BLOCK_PARTITION`
* `PLLMOD_BIN_BLOCK_CLV`
* `PLLMOD_BIN_BLOCK_TREE`
* `PLLMOD_BIN_BLOCK_CUSTOM`

* `PLLMOD_BIN_ACCESS_SEQUENTIAL`
* `PLLMOD_BIN_ACCESS_RANDOM`
* `PLLMOD_BIN_ACCESS_SEEK`

* `PLLMOD_BIN_INVALID_OFFSET`

* `PLLMOD_BIN_ATTRIB_UPDATE_MAP`
* `PLLMOD_BIN_ATTRIB_PARTITION_DUMP_CLV`
* `PLLMOD_BIN_ATTRIB_PARTITION_DUMP_WGT`
* `PLLMOD_BIN_ATTRIB_ALIGNED`

## Functions

* `FILE * pllmod_binary_create`
* `FILE * pllmod_binary_open`
* `FILE * pllmod_binary_append_open`
* `int pllmod_binary_close`
* `pll_block_map_t * pllmod_binary_get_map`
* `int pllmod_binary_partition_dump`
* `pll_partition_t * pllmod_binary_partition_load`
* `int pllmod_binary_clv_dump`
* `int pllmod_binary_clv_load`
* `int pllmod_binary_utree_dump`
* `pll_utree_t * pllmod_binary_utree_load`
* `int pllmod_binary_custom_dump`
* `void * pllmod_binary_custom_load`

## Error codes

* 4001: `PLLMOD_BIN_ERROR_BLOCK_MISMATCH`
* 4002: `PLLMOD_BIN_ERROR_BLOCK_LENGTH`
* 4003: `PLLMOD_BIN_ERROR_BINARY_IO`
* 4010: `PLLMOD_BIN_ERROR_INVALID_INDEX`
* 4011: `PLLMOD_BIN_ERROR_INVALID_SIZE`
* 4012: `PLLMOD_BIN_ERROR_LOADSTORE`
* 4020: `PLLMOD_BIN_ERROR_MISSING_BLOCK`
