/*
 Copyright (C) 2016 Diego Darriba, Alexey Kozlov

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

 /**
  * @file rtree_operations.c
  *
  * @brief Common functions for PLL modules
  *
  * @author Diego Darriba
  * @author Alexey Kozlov
  */
#include <stdarg.h>

#include "pll.h"
#include "pllmod_common.h"

/**
 * @brief Set pll error (pll_errno and pll_errmsg)
 *
 * @param[in] errno the error code
 * @param[in] errmsg_fmt formatted error message
 */
__attribute__((format(printf, 2, 3)))
void pllmod_set_error(int errno, const char* errmsg_fmt, ...)
{
  pll_errno = errno;

  va_list args;
  va_start(args, errmsg_fmt);
  vsnprintf(pll_errmsg, PLLMOD_ERRMSG_LEN, errmsg_fmt, args);
  va_end(args);
}

/**
 * Reset pll error and messages.
 *
 * Call this function whithin operations whose error status depends on `pll_errno`
 * such that no error leaks in from previous operations.
 */
void pllmod_reset_error()
{
  pll_errno = 0;
  strcpy(pll_errmsg, "");
}
