
#ifndef _NO_DMC_UTIL_H_
#define _NO_DMC_UTIL_H_

#include <sys/types.h>
#include <stdio.h>


static const int MAX_STR_LEN = 1024; // Max size for string var


char *strdupnodmc(const char *s);
char *trim(char *s);
int startsWith(const char *str, const char *pre);
char *replace(char *s, char old, char new);

off_t getSizeOfFile(FILE *f); // Size in bytes


#endif
