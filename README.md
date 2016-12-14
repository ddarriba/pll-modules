# pll-modules
High Level modules for the Low Level Phylogenetic Likelihood Library

## Introduction

## Clone & compile

PLL-Modules depends on the PLL library submodule (libs/libpll), so you have to clone the repository as follows:

```bash
git clone --recursive https://github.com/ddarriba/pll-modules
```

or

```
git clone https://github.com/ddarriba/pll-modules
git submodule update --init --recursive
```

To compile and install module with the bundled version `libpll`, run:

```bash
./install-with-libpll.sh [TARGET_DIR]
```

Alternatively, if you want to use a different version of `libpll`, you can run:

```bash
./autogen.sh
./configure CPPFLAGS="-Ipath/to/libpll" LDFLAGS="-Lpath/to/libpll"
make
make install    # as root, otherwise run: sudo make install
```

The library will be installed on the operating system's standard paths.  For
some GNU/Linux distributions it might be necessary to add that standard path
(typically `/usr/local/lib`) to `/etc/ld.so.conf` and run `ldconfig`.

Microsoft Windows compatibility was tested with a cross-compiler and seems to
work out-of-the-box using [MingW](http://www.mingw.org/).

## Documentation

PLL-Modules is documented with Doxygen. Run the following command to generate
the API description in HTML and LaTeX formats:

```bash
doxygen doxy.cfg
```

Afterwards, the html API documentation is available at `doc/html/index.html`.
Additionally, you can generate pdf documentation by compiling the LaTeX output:

```bash
cd doc/latex
make
```

## Usage examples 

Please refer to the [wiki page](https://github.com/ddarriba/pll-modules/wiki) and/or the [examples directory](https://github.com/ddarriba/pll-modules/tree/master/examples).

## Documentation

Please refer to the [wiki page](https://github.com/ddarriba/pll-modules/wiki).

## Available functionality

Below is a list of available modules in the current version.
Check each individual README.md file in the correspondent subdirectory for more information.


| Module    | Prefix       | Description             |
|-----------|--------------|-------------------------|
| binary    | pllmod_bin_  | Binary I/O              |
| msa       | pllmod_msa_  | MSA management          |
| optimize  | pllmod_opt_  | Optimization algorithms |
| tree      | pllmod_tree_ | Tree management         |
| util      | pllmod_util_ | Convenience functions   |
| algorithm | pllmod_algo_ | High level algorithms   |

