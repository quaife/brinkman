MEX=/Applications/MATLAB_R2017a.app/bin/mex
FORTRAN = ifort
# Fortran compiler
FFLAGS = -c

OBJS=\
bhcommon_oldtree.o\
bhfmm2dpart.o\
bhfmm2drouts.o\
bhrouts.o\
d2mtreeplot.o\
d2tstrcr_omp.o\
dlaran.o\
hkrand.o\
l2dterms.o\
laprouts2d.o\
prini.o

all: bhcommon_oldtree bhfmm2dpart bhfmm2drouts bhrouts d2mtreeplot d2tstrcr_omp dlaran hkrand l2dterms laprouts2d prini mex

mex: stokesDLPfmm.F $(OBJS)
	$(MEX) -O stokesDLPfmm.F $(OBJS) $(MEXLIB)

bhcommon_oldtree: bhcommon_oldtree.f
	$(FORTRAN) $(FFLAGS) bhcommon_oldtree.f

bhfmm2dpart: bhfmm2dpart.f
	$(FORTRAN) $(FFLAGS) bhfmm2dpart.f

bhfmm2dpart_dr: bhfmm2dpart_dr.f
	$(FORTRAN) $(FFLAGS) bhfmm2dpart_dr.f

bhfmm2drouts: bhfmm2drouts.f
	$(FORTRAN) $(FFLAGS) bhfmm2drouts.f

bhrouts: bhrouts.f
	$(FORTRAN) $(FFLAGS) bhrouts.f

d2mtreeplot: d2mtreeplot.f
	$(FORTRAN) $(FFLAGS) d2mtreeplot.f

d2tstrcr_omp: d2tstrcr_omp.f
	$(FORTRAN) $(FFLAGS) d2tstrcr_omp.f

dlaran: dlaran.f
	$(FORTRAN) $(FFLAGS) dlaran.f

hkrand: hkrand.f
	$(FORTRAN) $(FFLAGS) hkrand.f

l2dterms: l2dterms.f
	$(FORTRAN) $(FFLAGS) l2dterms.f

laprouts2d: laprouts2d.f
	$(FORTRAN) $(FFLAGS) laprouts2d.f

prini: prini.f
	$(FORTRAN) $(FFLAGS) prini.f

clean:
	rm *.o
	cp *.mex* ../src
	rm *.mex*
