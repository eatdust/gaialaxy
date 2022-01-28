//   This file is part of gaialaxy
//
//   Copyright (C) 2022 C. Ringeval
//   
//   gaialaxy is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   gaialaxy is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with gaialaxy.  If not, see <https://www.gnu.org/licenses/>.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <wcs.h>
#include <wcserr.h>
#include <wcsfix.h>
#include <wcsprintf.h>
#include <wcsunits.h>
#include <wcsutil.h>
#include <wcshdr.h>


//let's hide the structure from fortran and makes it global
struct wcsprm wcs;


void ini_wcsprm(int NAXIS, double CRPIX[], double CDELT[], double CRVAL[], \
		double *PC, char (*CUNIT)[72], char (*CTYPE)[72], char DATEREF[72])
{
  
  int i,j;
  int status;

  // to be compared with fortran input
  /* printf("NAXIS= %i\n",NAXIS);
     printf("CRPIX[0] %e\n",CRPIX[0]);
     printf("CRPIX[1] %e\n",CRPIX[1]);
     printf("CDELT[0] %e\n",CDELT[0]);
     printf("CDELT[1] %e\n",CDELT[1]);
     printf("CRVAL[0] %e\n",CRVAL[0]);
     printf("CRVAL[1] %e\n",CRVAL[1]);
  */

  wcs.flag = -1;
  
  wcsini(1,NAXIS,&wcs);

  for (i = 0; i < NAXIS; i++) {
    wcs.crpix[i] = CRPIX[i];
    wcs.cdelt[i] = CDELT[i];    
    wcs.crval[i] = CRVAL[i];
  }

  wcs.pc = PC;

  // to be compared with fortran input
  /* for (i = 0; i < NAXIS; i++) {
    for (j = 0; j < NAXIS; j++) {
      printf("wcs.pc= %e\n",*(wcs.pc++));
    }
  }
  */
  
  for (i = 0; i < NAXIS; i++) {
    strcpy(wcs.cunit[i], &CUNIT[i][0]);
  }

  for (i = 0; i < NAXIS; i++) {
    strcpy(wcs.ctype[i], &CTYPE[i][0]);
  }

  strcpy(wcs.dateref, DATEREF);
  
  
  status = wcsset(&wcs);
  if (status) {
    wcsprintf("\n");
    wcsperr(&wcs, 0x0);
  }

  return;
}


void destroy_wcsprm()
{
  wcsfree(&wcs);
}


void size_wcsprm(int sizes[2])
{
  wcssize(&wcs, sizes);
}


void fix_wcsprm()
{
  int stat[NWCSFIX];
  struct wcserr info[NWCSFIX];

  wcserr_enable(1);
  int status = wcsfixi(7, 0, &wcs, stat, info);

  
  // some logs
  /* wcsprintf("wcsfix status returns: (");
  for (int i = 0; i < NWCSFIX; i++) {
    wcsprintf(i ? ", %d" : "%d", stat[i]);
  }
  wcsprintf(")\n");
  */

  for (int i = 0; i < NWCSFIX; i++) {
    if (info[i].status < -1 || 0 < info[i].status) {
      wcsprintf("\n");
      wcserr_prt(info+i, 0x0);

      // Free memory used to store the message.
      if (info[i].msg) wcsdealloc(info[i].msg);
    }
  }

  if (status) {
    wcsprintf("\nfix_wcsprm error %d", status);
  }

  // Print error messages from the wcsprm structure
  if (wcsset(&wcs)) {
    wcsprintf("\n");
    wcsperr(&wcs, 0x0);
  }

  
  // Print the structure content
  /*
    wcsprintf("\n");
    wcsprt(&wcs);
    wcsprintf("\n------------------------------------"
    "------------------------------------\n");

  */

  return;
  
}

void hdo_alloc_wcsprm(char **header, int *nkeyrec)
{

  int status;
  char *hptr;
  int i;
  
  status = wcshdo(WCSHDO_none, &wcs, nkeyrec, header);

  if (status) {
    wcsprintf("\nhdo_wcsprm error %d", status);
  }


  // to compare with fortran output
  /*
  printf("nkeys= %i\n",*nkeyrec);    
  hptr = *header;
  printf("\n");
  for (i = 0; i < *nkeyrec; i++, hptr += 80) {
    printf("%.80s\n", hptr);
  }
  */

  return;
}

void hdo_free_wcsprm(char *header)
{
  wcsdealloc(header);  
  return;
}

void p2s_wcsprm(int ncoord, int nelem, const double pixcrd[ncoord][nelem], \
		double imgcrd[ncoord][nelem], double phi[ncoord], double theta[ncoord], \
		double world[ncoord][nelem])
{

  int status;
  int stat[ncoord];

  status = wcsp2s(&wcs,ncoord,nelem,pixcrd[0],imgcrd[0],phi,theta,world[0],stat);

  if (status) {
    wcsprintf("\np2s_wcsprm error %d", status);
  }
    

}


void s2p_wcsprm(int ncoord, int nelem, const double world[ncoord][nelem], \
		double phi[ncoord], double theta[ncoord], double imgcrd[ncoord][nelem], \
		double pixcrd[ncoord][nelem])
{
  int status;
  int stat[ncoord];

  status = wcss2p(&wcs,ncoord,nelem,world[0],phi,theta,imgcrd[0],pixcrd[0],stat);

  if (status) {
    wcsprintf("\ns2p_wcsprm error %d", status);
  }

}


void ccs_wcsprm(double lng2p1, double lat2p1, double lng1p2, \
		const char *clng, const char *clat, const char *radesys, \
		double equinox, const char *alt)  
{

  int status;
  
  // In case we need a copy of the previous wcsprm struct.
  /* struct wcsprm wcsold;
  
  wcsold.flag = -1;
  if (wcssub(1, &wcs, 0x0, 0x0, &wcsold)) {
    wcsperr(&wcsold, 0x0);
    return;
  }
  */

  status = wcsccs(&wcs, lng2p1, lat2p1, lng1p2, clng, clat, radesys, equinox, alt);

  if (status) {
    wcsprintf("\ns2p_wcsprm error %d", status);
  }
    
}
