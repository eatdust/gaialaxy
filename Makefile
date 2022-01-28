#   This file is part of gaialaxy
#
#   Copyright (C) 2022 C. Ringeval
#   
#   gaialaxy is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   gaialaxy is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with gaialaxy.  If not, see <https://www.gnu.org/licenses/>.


ext=$(shell uname | cut -c1-3)

FC=gfortran
CC=gcc
FFLAGS = -O2 -fopenmp
CFLAGS = -O2
INCLUDE= -I/usr/include/wcslib
LFLAGS= -L/usr/lib64 -lcfitsio -lwcs


OBJS= precision.o iofits.o gaiaconst.o iogaia.o wcswrap.o fwcs.o 

gaialaxy.$(ext) : $(OBJS) gaialaxy.o
	$(FC) $(FFLAGS) $(OBJS) gaialaxy.o $(LFLAGS) -o $@

%.o: %.F08
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.F03
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.F90
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.f03
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

clean:
	rm *.$(ext) *.o *.mod


