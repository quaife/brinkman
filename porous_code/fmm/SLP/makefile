MEX=/Applications/MATLAB_R2017a.app/bin/mex
#location of mex compiler on ronaldo

FORTRAN = gfortran
FFLAGS = -c 
#MEX=/Applications/MATLAB_R2015a.app/bin/mex
#location of mex compiler on ronaldo

OBJS=\
cfmm2dpart.o\
d2mtreeplot.o\
d2tstrcr_omp.o\
l2dterms.o\
laprouts2d.o\
lfmm2drouts.o\
prini.o\
rfmm2dpart.o\
stokesSLP.o

all: cfmm2dpart d2mtreeplot d2tstrcr_omp l2dterms laprouts2d lfmm2drouts prini rfmm2dpart stokesSLP mex

mex: stokesSLPfmm.F $(OBJS) 
	$(MEX) -compatibleArrayDims stokesSLPfmm.F $(OBJS) 

cfmm2dpart: cfmm2dpart.f
	$(FORTRAN) $(FFLAGS) cfmm2dpart.f

d2mtreeplot: d2mtreeplot.f
	$(FORTRAN) $(FFLAGS) d2mtreeplot.f

d2tstrcr_omp: d2tstrcr_omp.f
	$(FORTRAN) $(FFLAGS) d2tstrcr_omp.f

l2dterms: l2dterms.f
	$(FORTRAN) $(FFLAGS) l2dterms.f

laprouts2d: laprouts2d.f
	$(FORTRAN) $(FFLAGS) laprouts2d.f

lfmm2drouts: lfmm2drouts.f
	$(FORTRAN) $(FFLAGS) lfmm2drouts.f

prini: prini.f
	$(FORTRAN) $(FFLAGS) prini.f

rfmm2dpart: rfmm2dpart.f
	$(FORTRAN) $(FFLAGS) rfmm2dpart.f

stokesSLP: stokesSLP.f
	$(FORTRAN) $(FFLAGS) stokesSLP.f

clean: 
	rm *.o
