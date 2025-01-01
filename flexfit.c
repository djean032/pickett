/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   25 March 1999: read option cards with fgetstr */
/*   30 Dec.  1999: include changes for dlsq */
/*   10 Oct.  2001: change fit diverging code */
/*   21 Sept. 2002: fix NRJ criterion */
/*   18 Aug.  2003: code cleanup, @comment@ is for splint */ 

/**************************************************************************/
/*                                                                        */
/*   THIS IS A GENERALIZED LINE FITTING PROGRAM                           */
/*   IT FITS LINES TO PARAMETERS IN A MODEL HAMILTONIAN                   */
/*   BLENDED LINES ARE TREATED SPECIALLY:                                 */
/*     IF THE EXPTL.FREQ. IS THE SAME TO 1 HZ THEN ONLY THE INVERSE ERROR */
/*             AVERAGED FREQUENCIES ARE USED IN THE FIT FOR THE BLEND     */
/*                                                                        */
/**************************************************************************/
/* "trust region" Marquardt fitting is described in John. E. Dennis and   */
/* Robert B. Schnabel, Numerical Methods for Unsconstrained Optimization  */
/* and Non-linear Equations, Prentice-Hall, 1983.                         */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "calpgm.h"
#include "lsqfit.h"
#define NDCARD 130
#define PR_DELAY 6    /* seconds delay between informational messages */

/************** CALFIT interfaces ***********************************/
int qnfmt2(int nqn, short *qnum, /*@out@*/ char *aqnum);
int parer(double par, double errx, double dif, 
                 /*@out@*/ char *ptmp);
int linein(FILE * luin, int *nline, int iqnfmt);
int lineix(FILE * lu, int flg, int nline, int nblkpf, int iqnfmt);
/********************************************************************/
static char card[NDCARD];

int main(int argc, char *argv[])
{
  #define NFILE 5
  #define LBLEN 10
  #define FMT_xbgnMW "%5d: %s%14.5f%14.5f %10.5f %10.5f %9.5f"
  #define FMT_xblnMW "%14.5f %10.5f %6.4f\n"
  #define FMT_xbgnIR "%5d: %s%14.7f%14.7f %10.7f %10.7f %9.7f"
  #define FMT_xblnIR "%14.7f %10.7f %6.4f\n"
  static double zero = 0.;
  static double tiny = 1.5e-38;
  static const char *ext[NFILE] = { "par", "lin", "fit", "bak", "var" };
  enum efile {epar, elin, efit, ebak, evar};
  SXLINE *linebuf;
  bcd_t *idpar;
  double *par, *erp, *oldpar, *erpar, *dpar, *delbgn;
  double *fitbgn, *var, *oldfit, *fit, *teig, *pmix;
  int *iperm;
  char *parlbl;
  FILE *lupar, *lulin, *lufit, *lubak, *luvar;
  double *egy, *egyder;
  double *pvar, *pfit, *pfitb, *pfitd;
  char *fname[NFILE+1], *tlbl, *tlblnxt;
  double dvec[8], fqfac[4], varv[5], adif, cerr, afrq, rerr, xsqbest, cwid;
  double parfac, ex, xsqir, xsqmw, xerr, xfrq, xsqt, scale, avgir, xwid;
  double avgmw, fqfacq, xerrmx, dif, marqp[3], frq, val, xwt, marqlast, parfac0;
  double bigd, bigf, bigdIR, bigfIR, bigdMW, bigfMW;
  size_t nl, nlsq;
  int ifac, iblk, ndfit, ndiag, lblk, iflg, line, icnt, nfit, indx, nfir;
  int lstf, nitr, i, k, iblnd, lblnd, nline, initl, inpcor, n, marqflg;
  int nsize, nxpar, nxfit, lnext, noptn, ibase, nf, nblkpf, itd, ibcd, nfmt;
  int limlin, iqnfmt, maxdm, maxf, nrj, nqn, itr, npar, nsize_p, ndbcd;
  int supblnd, catqn, ndfree, ndfree0;
  short qnum[2*MAXQN];
  char ch, pare[64], aqnum[6*MAXQN+2], namfil[NDCARD];

  bigdIR = 9.9999999;
  bigdMW = 100. * bigdIR;
  bigfIR = 99999.9999999;
  bigfMW = 100. * bigfIR;
  fqfac[0] = 1;
  fqfac[1] = -1;
  marqflg = 0;
  nsize_p = maxmem(&nl);
/**
 * @note The following gets the filenames from the command line arguments.
 * The .par file is backed up to a .bak file.
 * Then the .fit and .bak files are opened for writing and reading, respectively.
 */
  filget(argc, argv, NFILE, fname, ext);
  filbak(fname[epar], fname[ebak]);
  lubak = fopenq(fname[ebak], "r");
  lufit = fopenq(fname[efit], "w");

/**
 * @note
 * The first line of the .par file is attempted to be read.
 */

  if (fgetstr(card, NDCARD, lubak) <= 0) {
    puts(" Unable to read title of .par file");
    exit(EXIT_FAILURE);
  }
/**
 * @note
 * The dvec array is initialized and then the second line of the .par file is attempted to be read.
 * If the second line read is successful then the values are read into dvec.
 */
  chtime(card, 82);
  fputs(card, lufit);
  puts(card);
  catqn = MAXCAT;
  dvec[0] = 100.;
  dvec[1] = 32767.;
  dvec[2] = 1.;
  dvec[3] = 0.;
  dvec[4] = 0.;
  dvec[5] = 1e6;
  dvec[6] = 1.;
  dvec[7] = 1.;
  n = fgetstr(card, NDCARD, lubak);
  if (n != 0)
    n = pcard(card, dvec, 8, NULL);
  if (n == 0) {
    puts(" Unable to read second line of .par file");
    exit(EXIT_FAILURE);
  }
  /**
   * @note
   * The values of the dvec array are assigned to the corresponding variables.
   * number of params dvec[0]
   * number of lines dvec[1]
   * number of iterations dvec[2]
   * number of excluded params dvec[3]
   * marquardt parameter dvec[4]
   * error max dvec[5]
   * scaling factor for paramter errors dvec[6]
   * scaling factor for frequency errors dvec[7]
   */
    npar = (int) dvec[0];
    if (dvec[1] < 0.) {
      catqn = MAXQN; dvec[1] = -dvec[1];
    }
    nl = (size_t) dvec[1];
    nitr = (int) dvec[2];
    nxpar = (int) dvec[3];
    marqp[0] = dvec[4];
    xerrmx = dvec[5];
    parfac = dvec[6];
    parfac0 = parfac;
    fqfacq = dvec[7];
    fqfac[2] = fqfacq / 29979.2458;
    fqfac[3] = -fqfac[2];
    marqp[1] = -1;
    if (marqp[0] < 0.)
      marqp[0] = 0.;
    limlin = nline = (int) nl;
    if (fabs(parfac - 1.) > 1e-10) {
      fprintf(lufit, "PARAMETER ERRORS SCALED BY %15.6f", fabs(parfac));
      if (parfac < 0.) fputs(" times the standard error",lufit);
      fputc('\n', lufit);
    }
    return 0;
}

int qnfmt2(int nqn, short *qnum, char *aqnum)
{
  /* Local variables */
  int i;

  /*  formats quantum numbers for output */
  /*     NQN   = number of quanta */
  /*     QNUM  = vector of quanta */
  /*     AQNUM = string of quanta in character form */
  for (i = 0; i < nqn; ++i) {
    sprintf(aqnum, "%3d", (int) qnum[i]);
    aqnum += 3;
  }
  for (i = nqn; i < 12; ++i) {
    aqnum[2] = aqnum[1] = aqnum[0] = ' ';
    aqnum += 3;
  }
  aqnum[0] = '\0';
  return 0;
}                               /* qnfmt2 */

int parer(par, errx, dif, ptmp)
double par, errx, dif;
char *ptmp;
{
  static int czero = (int) '0';
  char *pfmt;
  double adif, apar, aten, aerr;
  char chexp[6], fmt[34];
  int msd, id, ie, efield, ip, lsd, k;


  /*      sets up special format for parameters and errors */
  /*     PAR  = parameter value */
  /*     ERRX = parameter error */
  /*     DIF  = parameter change */
  /*     PTMP = output string for printing */

  apar = par;
  aerr = errx;
  adif = dif;
  efield = 0;
  aten = 1.;
  /*     compute exponent fields */
  ie = (int) (log10(fabs(aerr) + 1.e-37) - 102.5) + 100;
  id = (int) (log10(fabs(adif) + 1.e-37) - 100.0) + 100;
  ip = (int) (log10(fabs(apar) + 1.e-37) - 100.0) + 100;
  lsd = -ie;
  if (lsd < 0)
    lsd = 0;
  msd = (ip > id) ? ip : id;
  /*  check for too many digits */
  k = 14 - ip;
  if (k < lsd)
    lsd = k;
  k = 10 - id;
  if (k < lsd)
    lsd = k;
  if (msd <= -2) {              /* number too small without exponent */
    k = (1 - msd) / 3;
    efield = -3 * k;
    while ((--k) >= 0)
      aten *= 1000;
  } else if (lsd < 0) {         /* number too big without exponent */
    k = (1 + msd) / 3;
    if (k > 0)
      efield = 3 * k;
    while ((--k) >= 0)
      aten *= 0.001;
  }
  if (efield != 0) {            /* E format */
    lsd += efield;
    memcpy(chexp, "0fE+00", 6);
    if (efield < 0) {
      chexp[3] = '-';
      efield = -efield;
    }
    msd = efield / 10;
    if (msd > 0) {
      efield -= msd * 10;
      chexp[4] = (char) (msd + czero);
    }
    chexp[5] = (char) (efield + czero);
    apar *= aten;
    aerr *= aten;
    adif *= aten;
  } else {                      /* F format */
    memcpy(chexp, "0f    ", 6);
  }
  if (lsd > 9)
    lsd = 9;
  if (lsd > 0)
    chexp[0] = (char) (lsd + czero);
  while ((lsd--) > 0)
    aerr *= 10.;
  ie = (int) (aerr + 0.5);
  pfmt = fmt;
  memcpy(pfmt, "%16.", 4);
  pfmt += 4;
  memcpy(pfmt, chexp, 2);
  pfmt += 2;
  memcpy(pfmt, "(%3d)", 5);
  pfmt += 5;
  memcpy(pfmt, chexp + 2, 4);
  pfmt += 4;
  memcpy(pfmt, " %12.", 5);
  pfmt += 5;
  memcpy(pfmt, chexp, 6);
  pfmt += 6;
  *pfmt = '\0';
  sprintf(ptmp, fmt, apar, ie, adif);
  return 0;
}                               /* parer */


int linein(luin, nline, iqnfmt)
FILE *luin;
int *nline;
int iqnfmt;
{
  /* Local variables */
  SXLINE *xline;
  double xfrqn, xerrn, xwtn, xfrqx, xerrx;
  int nqn, nqnu, nqnl, kqnu, kqnl, i, iqf, ipace, mxline, mxqn, isblnd, icmp;
  short nbln, nqnt[20], *iqnum;

  /*   get lines from input  and stores them */

  /*     LUIN= unit for finding lines */
  /*     NLINE = number of lines */
  /*     IQNFMT= qunatum number format for line input */
  /*     RETURN: largest quantum number */
  /*******************************************************************/

  mxline = *nline;
  mxqn = 1;
  nbln = 1;
  nqn = deflin(iqnfmt, nqnt);
  nqnu = nqn - 1;
  if (nqnt[nqnu] < 0)
    nqnu = 0;
  kqnu = nqnt[nqnu];
  nqnl = nqnu + nqn;
  kqnl = nqnt[nqnl];
  ipace = 100;
  xfrqx = xerrx = 0.;
  icmp = 0;
  for (i = 1; i <= mxline; ++i) {       /*  loop for reading lines */
    xline = lbufof(1, i);
    iqnum = xline->qn;
    if (getlin(luin, nqn, nqnt, iqnum, &xfrqn, &xerrn, &xwtn, 
               card, NDCARD) < 0) {
      *nline = i - 1;
      return mxqn;
    }
    iqf = iqnum[nqnu];
    if (iqf == -1) {
      if (kqnu >= 0) {
        iqf = -iqnum[kqnu];
        if (iqf >= 0)
          iqf = -1;
      }
      iqnum[nqnu] = (short) iqf;
    }
    if (iqf < 0)
      iqf = -iqf;
    if (iqf > mxqn)
      mxqn = iqf;
    iqf = iqnum[nqnl];
    if (iqf == -1) {
      if (kqnl > 0) {
        iqf = -iqnum[kqnl];
        if (iqf >= 0)
          iqf = -1;
      }
      iqnum[nqnl] = (short) iqf;
    }
    if (iqf < 0)
      iqf = -iqf;
    if (mxqn < iqf)
      mxqn = iqf;
    xline->xfrq = xfrqn;
    xline->xerr = (float) xerrn;
    xline->xwt = (float) fabs(xwtn);
    xline->linku = 0;
    xline->linkl = 0;
    isblnd = 0;
    if (icmp != 0 && fabs(xfrqn - xfrqx) < fabs(xfrqn) * 1.e-14 + 1.e-8) { 
      /* frq match */
      if (fabs(xerrn - xerrx) < 1e-7) {
        isblnd = 1;
      } else if ((xerrn / xerrx) > 2.0 && nbln > 2) {
        isblnd = 1; ++nbln; icmp = 0;
        xline->xwt = (float)0.;
        iqnum[0] = (short)-1;
        iqnum[nqn] = iqnum[0];
      }
    }
    if (isblnd != 0) {
      xline->bln = nbln;
      xline = lbufof(1, i - 1);
      xline->bln = -2;
      nbln += 2;
    } else {
      xline->bln = 0;
      nbln = 2; icmp = 1;
    }
    if (ipace <= i) {
      ipace += 100;
      printf("Reading Line %d\n", i);
      fflush(stdout);
    }
    xerrx = xerrn;
    xfrqx = xfrqn;
  }
  return mxqn;
}                               /* linein */

int lineix(lu, flg, nline, nblkpf, iqnfmt)
FILE *lu;
int flg, nline, nblkpf, iqnfmt;
{  /*   get lines from input and store them */
  /*     LU = unit for printout of lines ( if > 0 ) */
  /*     NLINE = number of lines */
  /*     NBLKPF= number of blocks per F */
  /*     IQNFMT= qunatum number format for line input */
  /******************************************************************/
  static int nsort = 2048;
  SXLINE *xline;
  double xfrqn, xerrn, xwtn, xnorm;
  int nblk, ipos, i, j, ipace, nread, iblkl, iblku, ncat;
  int linkx, indxl, linky, indxu, orgblk, nqn, nqn2, nbad;
  /*@owned@*/ int *prvblk;
  short *iqnum;
  char aqnum[6*MAXQN+2];

  nbad = 0;
  nblk = 0;
  prvblk = (int *) mallocq((size_t) (nsort + 1) * sizeof(int));
  prvblk[0] = 0;
  for (i = 1; i <= nsort; ++i) {
    prvblk[i] = 0;
  }
  nqn = iqnfmt % 10;
  if (nqn == 0) nqn = 10;
  nqn2 = nqn + nqn; ncat = nqn2;
  if (ncat < 12) ncat = 12;
  i = (iqnfmt / 100) % 5;
  if (i >= nqn) {
    ipos = 1;
  } else {
    ipos = nqn;
  }
  if (flg < 0) {
    fputs(" LINE,BLKU,INDXU,BLKL,INDXL,QUANTUM NUMBERS", lu);
    for (i = 0; i < 19; ++i)
      fputc(' ', lu);
    fputs("ENERGY    EXP. ERROR    WEIGHTS\n", lu);
  }
  xnorm = 0.;
  ipace = 50;
  /*       loop for converting lines */
  for (nread = 1; nread <= nline; ++nread) {
    xline = lbufof(1, nread);
    xfrqn = xline->xfrq;
    xerrn = xline->xerr;
    xwtn = xline->xwt;
    /* find blocks and index for upper and lower states */
    iqnum = xline->qn;
    getblk(&iblku, &indxu, iqnum, nblkpf, ipos, nqn);
    getblk(&iblkl, &indxl, &iqnum[nqn], nblkpf, ipos, nqn);
    if (iblkl == 0 && iqnum[nqn] >= 0)
      iblku = 0;
    xline->ibu = iblku;
    xline->inu = (short) indxu;
    xline->ibl = iblkl;
    xline->inl = (short) indxl;
    if (iblku == 0 && (xline->bln & 1) == 0) {
      /*  print out bad line and try for next */
      ++nbad;
      xline->xwt = 0.;
      xwtn = 0.;
      qnfmt2(nqn2, iqnum, aqnum);
      printf(    "Bad Line(%3d): %s %14.5f %8.5f\n",
                 nread, aqnum, xfrqn, xerrn);
      fprintf(lu,"Bad Line(%3d): %s %14.5f %8.5f\n",
                 nread, aqnum, xfrqn, xerrn);
    } else {
      /*  set up links for calculating in order of block */
      if (iblku <= iblkl) {
        lnlink(prvblk, nsort, iblku, nread);
        lnlink(prvblk, nsort, iblkl, -nread);
        if (nblk < iblkl)
          nblk = iblkl;
      } else {
        lnlink(prvblk, nsort, iblkl, -nread);
        lnlink(prvblk, nsort, iblku, nread);
        if (nblk < iblku)
          nblk = iblku;
      }
      if (flg < 0) {
        iqnum = xline->qn;
        fprintf(lu," %4d%4d%4d%4d%4d:", nread, iblku, indxu, iblkl, indxl);
        for (i = 0; i < ncat; ++i) {
          j = iqnum[i];
          fprintf(lu, "%3d", j);
        }
        fprintf(lu, " %14.4f %9.4f %9.4f", xfrqn, xerrn, xwtn);
        j = xline->bln;
        if (j != 0) {
          fprintf(lu, "   Line Blended with %3d\n", nread - (j >> 1));
        } else {
          fputc('\n', lu);
        }
      }
    }
    /* let the user know something is happening */
    if (ipace <= nread || nread == nline) {
      ipace += 50;
      printf("Converting Line %d\n", nread);
      fflush(stdout);
    }
    j = xline->bln;
    if (j != 0) {
      xnorm += xwtn;      
      if (j > 0) {     /* normalize weights */
        xnorm = 1. / xnorm; 
        for (j = nread - (j >> 1); j <= nread; ++j) {
          xline = lbufof(1, j);
          xline->xwt = (float) (xline->xwt * xnorm);
        }
        xnorm = 0.;
      }
    } else {
      xline->xwt = 1.;
    }
  }                             /* end loop for converting lines */
  orgblk = nsort;
  while (nblk > orgblk) {       /* finish up links */
    --orgblk;
    linkx = 0;
    for (j = 0; j < nsort; ++j) {
      if (prvblk[j] != 0) {
        linkx = prvblk[j];
        prvblk[j] = 0;
      }
    }
    prvblk[0] = linkx;
    prvblk[nsort] = 0;
    linkx = lnlink(prvblk, nsort, 1, 0);
    while (linkx != 0) {
      linky = linkx;
      getdbk(&linkx, &iblkl, &j, &j, &j);
      iblkl -= orgblk;
      lnlink(prvblk, nsort, iblkl, linky);
    }
    orgblk += nsort;
  }
  free(prvblk);
  return nbad;
}                               /* lineix */
