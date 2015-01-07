// Author: Mingcheng Chen (linyufly@gmail.com)

#include "util.h"

#include <cstdio>
#include <cstdarg>
#include <cstdlib>

namespace {

void report_va_error(const char *format, va_list args) {
  vfprintf(stderr, format, args);

  exit(EXIT_FAILURE);
}

}

void report_error(const char *format, ...) {
  va_list args;
  va_start(args, format);
  report_va_error(format, args);
  va_end(args);
}
