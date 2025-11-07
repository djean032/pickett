#ifndef SLIB_H
#define SLIB_H

#ifndef FILE
#include <stdio.h>
#endif
/************** SLIB interface ***************************************** */
#ifdef __cplusplus
extern "C"
{
#endif

  int chtime (char str[], const int n);
  int caldelay (int delay);
  /*@out@*/ /*@only@*/ /*@notnull@*/ void *mallocq (size_t nl);
  int maxmem (/*@null@*/ /*@out@*/ size_t *nl);
  int filget (const int argc, char *argv[], const int nfile,
              /*@out@*/ char *cfil[], const char *cext[]);
  int fgetstr (/*@out@*/ char buffer[], const int n, FILE *stream);
  int rqexit (const int ival);
  void brkqr (int);

  /*@dependent@*/ FILE *fopenq (const char *fname, const char *opt)
      /*@modifies fileSystem@*/;

#ifdef __cplusplus
}
#endif

#endif
