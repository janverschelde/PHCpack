ALL: default
##### User configurable options #####

SHELL       = /bin/sh
ARCH        = LINUX
MPIR_HOME   = /huis/cluster/mpich-1.2.4/mpich-1.2.4
CC          = /huis/cluster/mpich-1.2.4/mpich-1.2.4/bin/mpicc
CLINKER     = $(CC)
F77         = /huis/cluster/mpich-1.2.4/mpich-1.2.4/bin/mpif77
FLINKER     = $(F77)
F90         = /huis/cluster/mpich-1.2.4/mpich-1.2.4/bin/mpif90
F90LINKER   = $(F90)      
CCC         = /huis/cluster/mpich-1.2.4/mpich-1.2.4/bin/mpiCC
CCLINKER    = $(CCC)
AR          = ar crl
RANLIB      = ranlib
LOG_LIB     = -mpilog -lm
PROF_LIB    = -lmpe -lm
OPTFLAGS    = 
MPE_DIR     = /huis/cluster/mpich-1.2.4/mpich-1.2.4/mpe
MAKE        = make --no-print-directory
#
srcdir=/huis/cluster/mpich-1.2.4/mpich-1.2.4/examples/basic

VPATH=.:$(srcdir)
### End User configurable options ###
.SUFFIXES:	.cc .C .f90

CFLAGS	  = $(OPTFLAGS) 
CFLAGSMPE = $(CFLAGS) -I$(MPE_DIR)/include 
CCFLAGS	  = $(CFLAGS)
#FFLAGS	  = '-qdpc=e' 
FFLAGS	  = $(OPTFLAGS)
EXECS	  = adatest dbltest ts_phc_sys_rw phctest broadcast
ALL_EXECS = ${EXECS}

# This is added stuff

PHCLIB=../Ada
BIN=/tmp
INCLULIBS=-I$(PHCLIB)/System \
 -I$(PHCLIB)/Math_Lib/Numbers -I$(PHCLIB)/Math_Lib/Matrices \
 -I$(PHCLIB)/Math_Lib/Polynomials -I$(PHCLIB)/Math_Lib/Supports \
 -I$(PHCLIB)/Homotopy -I$(PHCLIB)/Continuation \
 -I$(PHCLIB)/Root_Counts/Product -I$(PHCLIB)/Root_Counts/Implift \
 -I$(PHCLIB)/Root_Counts/Stalift -I$(PHCLIB)/Root_Counts/Dynlift \
 -I$(PHCLIB)/Root_Counts/Symmetry -I$(PHCLIB)/Root_Counts/ExtSoft \
 -I$(PHCLIB)/Schubert/SAGBI -I$(PHCLIB)/Schubert/Pieri \
 -I$(PHCLIB)/Schubert/Induction \
 -I$(PHCLIB)/Components/Samplers -I$(PHCLIB)/Components/Interpolators \
 -I$(PHCLIB)/Components/Factorization -I$(PHCLIB)/Components/Decomposition \
 -I$(PHCLIB)/Differentials -I$(PHCLIB)/Main -I$(PHCLIB)/CtoPHC
GNATFLAGS=-gnatv -O3 -gnatp

all: adatest dbltest ts_phc_sys_rw

# how to make the executables :

adatest: adatest.o
	gcc -c call_hello.c
	gnatmake -c hello.adb
	gnatbind -n -C hello
	gnatlink -C hello call_hello.o adatest.o -o adatest \
                 --LINK=$(CLINKER)

broadcast: broadcast.o
	$(CLINKER) $(OPTFLAGS) -o broadcast broadcast.o -lm

dbltest: dbltest.o
	$(CLINKER) $(OPTFLAGS) -o dbltest dbltest.o -lm

ts_phc_sys_rw:
	gcc -c $(PHCLIB)/CtoPHC/ts_phc_sys_rw.c
	gnatmake -c $(INCLULIBS) $(GNATFLAGS) phc_sys_rw.adb
	gnatbind -n -C phc_sys_rw
	gnatlink -C phc_sys_rw ts_phc_sys_rw.o -o ts_phc_sys_rw

phctest: phctest.o
	gcc -c call_phc_rw.c
	gnatmake -c $(INCLULIBS) $(GNATFLAGS) phc_sys_rw.adb
	gnatbind -n -C phc_sys_rw
	gnatlink -C phc_sys_rw call_phc_rw.o phctest.o -o phctest \
                 --LINK=$(CLINKER)

# cleaning up everything :

clean:
	rm -f *.o *ali *~ PI* $(ALL_EXECS) core
