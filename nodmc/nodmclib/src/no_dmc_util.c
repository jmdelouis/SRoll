/******************************************************************************
 * "no_dmc_util.c"
 * 
 * author:  Christan Madsen
 * date:    2015-01-06 (initial)
 * version: STABLE
 *****************************************************************************/

#define _POSIX_SOURCE


#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>

#include "no_dmc_util.h"


//-----------------------------------------------------------------------------
// String manipulation
//-----------------------------------------------------------------------------

/*
 * Local implementation of strdup(). This definition prevent to define
 * specific compilation tag...
 */
char *strdupnodmc(const char *s) {
  int n = strlen(s);
  char *out = malloc((n+1) * sizeof(char));

  if (out == NULL) {
    return NULL;
  }

  strcpy(out, s);

  return out;
}

/*
 * Allow to remove leading and trailing space.
 * Note that this funtion does not allocate a new string!
 */
char *trim(char *s) {
  char *p = s;
  int l = strlen(p);

  while(isspace(p[l - 1])) p[--l] = 0;
  while(* p && isspace(* p)) ++p, --l;

  memmove(s, p, l + 1);

  return s;
}

/*
 * Returns 0 if 'str' starts with the prefix 'pre', another non zero value otherwise.
 */
int startsWith(const char *str, const char *pre) {
  size_t lenpre = strlen(pre),
         lenstr = strlen(str);
  return lenstr < lenpre ? -1 : strncmp(pre, str, lenpre);
}

/*
 * Replace an 'old' char by a 'new' one in a specified string 'str'.
 * No further allocation required!
 * Input string is modified!
 */
char *replace(char *s, char old, char new)
{
  char *p = s;

  while(*p)
  {
    if(*p == old)
      *p = new;

    ++p;
  }

  return s;
}


//-----------------------------------------------------------------------------
// File manipulation
//-----------------------------------------------------------------------------

/*
 * Returns the size (in bytes) of the file stream pass in arg.
 */
off_t getSizeOfFile(FILE *f) {
  int fd = fileno(f); // to retrieve file descriptor from stream.
  struct stat st;
  fstat(fd, &st);

  return st.st_size;
}
