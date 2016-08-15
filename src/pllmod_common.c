#include <stdarg.h>

#include "pll.h"
#include "pllmod_common.h"

void pllmod_set_error(int errno, const char* errmsg_fmt, ...)
{
  pll_errno = errno;

  va_list args;
  va_start(args, errmsg_fmt);
  vsnprintf(pll_errmsg, PLLMOD_ERRMSG_LEN, errmsg_fmt, args);
  va_end(args);
}
