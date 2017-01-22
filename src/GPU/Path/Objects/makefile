# This makefile works on RHEL 6.6 with GNATGPL 2014, NVIDIA CUDA 6.5.

# Location of the QD library and its GPU version:

QD_ROOT=/usr/local/qd-2.3.17
QD_LIB=/usr/local/lib
GQD_HOME=/usr/local/gqd_1_2
SDK_HOME=/usr/local/cuda/sdk

# QD_ROOT=/usr/local/qd-2.3.13
# QD_LIB=/usr/local/lib
# GQD_HOME=/usr/local/gqd_1_2
# SDK_HOME=/usr/local/cuda/sdk

gpp=/usr/bin/g++

QDFLAG=-I$(QD_ROOT)/include
GQDFLAG=-I$(QD_ROOT)/include -I/usr/local/cuda/include -I$(GQD_HOME)/inc

CFLAG   = ../Complex
CDFLAG  = ../DefineTypesD
CDDFLAG = ../DefineTypesDD
CQDFLAG = ../DefineTypesQD
PATHCPU = ../Path_CPU
PATHGPU = ../Path_GPU
PATHTEST= ../Main
POLY    = ../Poly
BIN     = ../bin
TEST    = ../Test

IFLAG=-I./$(POLY) -I./$(PATHTEST) -I./$(PATHCPU) -I./$(PATHGPU) -I./$(CFLAG)

PHCv2 = ./../../..
PHCLib=$(PHCv2)/Lib
PHCAda=$(PHCv2)/Ada
MPILib=$(PHCv2)/MPI
PHCgpu=$(PHCv2)/GPU
PHCBIN=$(PHCv2)/bin
INCLULIBS=-I$(PHCAda)/System -I$(PHCAda)/System/Unix_Timer \
 -I$(PHCAda)/Math_Lib/Numbers -I$(PHCAda)/Math_Lib/QD \
 -I$(PHCAda)/Math_Lib/Vectors -I$(PHCAda)/Math_Lib/Matrices \
 -I$(PHCAda)/Math_Lib/Divisors -I$(PHCAda)/Math_Lib/Reduction \
 -I$(PHCAda)/Math_Lib/Polynomials -I$(PHCAda)/Math_Lib/Functions \
 -I$(PHCAda)/Math_Lib/Circuits -I$(PHCAda)/Math_Lib/Supports \
 -I$(PHCAda)/Deformations/Solutions -I$(PHCAda)/Deformations/Homotopy \
 -I$(PHCAda)/Deformations/Newton -I$(PHCAda)/Deformations/Curves \
 -I$(PHCAda)/Deformations/End_Games -I$(PHCAda)/Deformations/Sweep \
 -I$(PHCAda)/Deformations/Trackers -I$(PHCAda)/Deformations/Continuation \
 -I$(PHCAda)/Root_Counts/Product -I$(PHCAda)/Root_Counts/Binomials \
 -I$(PHCAda)/Root_Counts/Implift -I$(PHCAda)/Root_Counts/Stalift \
 -I$(PHCAda)/Root_Counts/Dynlift -I$(PHCAda)/Root_Counts/Symmetry \
 -I$(PHCAda)/Root_Counts/MixedVol -I$(PHCAda)/Root_Counts/Puiseux \
 -I$(PHCAda)/Schubert/SAGBI -I$(PHCAda)/Schubert/Pieri \
 -I$(PHCAda)/Schubert/Induction \
 -I$(PHCAda)/Components/Samplers -I$(PHCAda)/Components/Interpolators \
 -I$(PHCAda)/Components/Factorization -I$(PHCAda)/Components/Decomposition \
 -I$(PHCAda)/Components/Solver -I$(PHCAda)/Components/Tropical \
 -I$(PHCAda)/CtoPHC/Funky -I$(PHCAda)/CtoPHC/State -I$(MPILib) \
 -I$(PHCAda)/PHCtoC -I$(PHCAda)/Tasking -I$(PHCAda)/Main
# flags with full optimization
GNATFLAGS=-gnatv -O3 -gnatp -gnatf
GNAT95FLAGS=-gnat95 -gnatv -O3 -gnatp -gnatf

all: path_d path_dd path_qd

test:
	$(gpp) -O2 $(TEST)/test.cpp -o test

gqd_qd_util.o: $(CFLAG)/gqd_qd_util.h $(CFLAG)/gqd_qd_util.cpp
	@-echo ">>> compiling utilities ..."
	$(gpp) -O2 $(GQDFLAG) -c $(CFLAG)/gqd_qd_util.cpp 
            
# D

workspace_host_d.o: $(PATHCPU)/workspace_host.h $(PATHCPU)/workspace_host.cpp
	$(gpp) -O2 -I $(CDFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/workspace_host.cpp -o workspace_host_d.o

eval_host_d.o: $(PATHCPU)/eval_host.h $(PATHCPU)/eval_host.cpp $(CDFLAG)/*.h  \
        path_data_d.o poly_d.o workspace_host_d.o
	$(gpp) -O2 -I $(CDFLAG) $(IFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/eval_host.cpp -o eval_host_d.o
	         
path_data_d.o: $(PATHTEST)/path_data.cpp $(PATHTEST)/path_data.h $(CDFLAG)/*.h
	$(gpp) -O2 -I $(CDFLAG)  $(IFLAG) $(QDFLAG) \
               -c $(PATHTEST)/path_data.cpp -o path_data_d.o

path_host_d.o: $(PATHCPU)/path_host.h $(PATHCPU)/path_host.cpp \
               predictor_host_d.o newton_host_d.o
	$(gpp) -O2  -I $(CDFLAG) $(IFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/path_host.cpp -o path_host_d.o

newton_host_d.o: $(PATHCPU)/newton_host.h $(PATHCPU)/newton_host.cpp \
        eval_host_d.o mgs_host_d.o 
	$(gpp) -O2 -I $(CDFLAG) $(IFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/newton_host.cpp -o newton_host_d.o

predictor_host_d.o: $(PATHCPU)/predictor_host.h $(PATHCPU)/predictor_host.cpp \
        workspace_host_d.o
	$(gpp) -O2 -I $(CDFLAG) $(QDFLAG) \
               -c $(PATHCPU)/predictor_host.cpp -o predictor_host_d.o

mgs_host_d.o: $(PATHCPU)/mgs_host.h $(PATHCPU)/mgs_host.cpp 
	$(gpp) -O2 -I $(CDFLAG) $(GQDFLAG) \
	       -c $(PATHCPU)/mgs_host.cpp -o mgs_host_d.o

families_d.o: $(POLY)/families.h $(POLY)/families.cpp
	$(gpp) -O2  -I $(CDFLAG) $(QDFLAG) \
	       -c $(POLY)/families.cpp -o families_d.o

poly_d.o: $(POLY)/poly.h $(POLY)/poly.cpp utilities_d.o $(CDFLAG)/*.h
	$(gpp) -O2  -I $(CDFLAG) $(QDFLAG) \
	       -c $(POLY)/poly.cpp -o poly_d.o

polysol_d.o: $(POLY)/poly.h $(POLY)/polysol.cpp utilities_d.o $(CDFLAG)/*.h
	$(gpp) -O2  -I $(CDFLAG) $(QDFLAG) \
	       -c $(POLY)/polysol.cpp -o polysol_d.o

path_gpu_d.o: $(PATHGPU)/*.cu $(PATHGPU)/*.h $(CDFLAG)/*.h
	@-echo ">>> compiling kernels ..."
	nvcc -O2 --define-macro path_precision=cd \
	         -I $(CDFLAG) $(IFLAG)\
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
	         -o path_gpu_d.o \
             -c $(PATHGPU)/path_gpu.cu --ptxas-options=-v

utilities_d.o: $(POLY)/utilities.h $(POLY)/utilities.cpp
	$(gpp) -O2  -I $(CDFLAG)  $(GQDFLAG)\
	       -c $(POLY)/utilities.cpp -o utilities_d.o 
	         
err_check_d.o: $(PATHTEST)/err_check.cpp $(PATHTEST)/err_check.h $(CDFLAG)/*.h
	$(gpp) -O2  -I $(CDFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/err_check.cpp -o err_check_d.o
	         
eval_test_d.o: $(PATHTEST)/eval_test.cpp $(PATHTEST)/err_check.h $(CDFLAG)/*.h
	$(gpp) -O2  -I $(CDFLAG)  $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/eval_test.cpp -o eval_test_d.o
	         
mgs_test_d.o: $(PATHTEST)/mgs_test.cpp $(PATHTEST)/err_check.h \
              $(CDFLAG)/*.h $(PATHGPU)/mgs.cu
	$(gpp) -O2  -I $(CDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/mgs_test.cpp -o mgs_test_d.o
	         
newton_test_d.o: $(PATHTEST)/newton_test.cpp $(PATHTEST)/err_check.h \
                 $(CDFLAG)/*.h 
	$(gpp) -O2  -I $(CDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/newton_test.cpp -o newton_test_d.o
	         
path_test_d.o: $(PATHTEST)/path_test.cpp $(PATHTEST)/path_test.h \
               $(CDFLAG)/*.h path_host_d.o
	$(gpp) -O2  -I $(CDFLAG)  $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/path_test.cpp -o path_test_d.o
	         
path_multi_test_d.o: $(PATHTEST)/path_multi_test.cpp $(PATHTEST)/err_check.h $(CDFLAG)/*.h path_test_d.o
	$(gpp) -O2  -I $(CDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/path_multi_test.cpp -o path_multi_test_d.o
	         
predict_test_d.o: $(PATHTEST)/predict_test.cpp $(PATHTEST)/err_check.h $(CDFLAG)/*.h
	$(gpp) -O2  -I $(CDFLAG) $(IFLAG) $(GQDFLAG)\
	         -c $(PATHTEST)/predict_test.cpp -o predict_test_d.o
	         
pieri_test_d.o: $(PATHTEST)/pieri_test.cpp $(PATHTEST)/err_check.h $(CDFLAG)/*.h
	$(gpp) -O2  -I $(CDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/pieri_test.cpp -o pieri_test_d.o
	         
# (PATHTEST)/init_test.cpp $(PATHTEST)/err_check.h $(CDFLAG)/*.h

init_test_d.o: poly_d.o
	$(gpp) -O2  -I $(CDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/init_test.cpp -o init_test_d.o
	         
witness_set_test_d.o: $(PATHTEST)/witness_set_test.cpp $(PATHTEST)/err_check.h $(CDFLAG)/*.h
	$(gpp) -O2  -I $(CDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/witness_set_test.cpp -o witness_set_test_d.o
	         
main_test_d.o: $(PATHTEST)/main_test.cpp $(PATHCPU)/parameter.h families_d.o \
               path_host_d.o path_gpu_d.o newton_host_d.o err_check_d.o \
               eval_test_d.o mgs_test_d.o newton_test_d.o predict_test_d.o \
               path_test_d.o path_multi_test_d.o pieri_test_d.o init_test_d.o \
               witness_set_test_d.o
	$(gpp) -Dcompilewgpp=1 -O2 -I$(CDFLAG) -I$(PHCLib) $(IFLAG) $(GQDFLAG) \
	       -c $(PATHTEST)/main_test.cpp -o main_test_d.o
	         
path_d: families_d.o main_test_d.o poly_d.o polysol_d.o utilities_d.o \
        path_gpu_d.o gqd_qd_util.o path_host_d.o eval_host_d.o mgs_test_d.o
	$(gpp) -O2 $(GQDFLAG) \
		main_test_d.o families_d.o poly_d.o polysol_d.o utilities_d.o \
                path_gpu_d.o gqd_qd_util.o predictor_host_d.o mgs_host_d.o \
                newton_host_d.o path_host_d.o path_data_d.o eval_host_d.o \
		workspace_host_d.o err_check_d.o mgs_test_d.o \
                predict_test_d.o eval_test_d.o path_test_d.o \
                path_multi_test_d.o newton_test_d.o pieri_test_d.o \
                init_test_d.o witness_set_test_d.o \
		$(QD_LIB)/libqd.a -o $(BIN)/path_d \
		-lcutil_x86_64 -lcudart \
		-L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib
            
# DD

workspace_host_dd.o: $(PATHCPU)/workspace_host.h $(PATHCPU)/workspace_host.cpp
	$(gpp) -O2 -I $(CDDFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/workspace_host.cpp -o workspace_host_dd.o

eval_host_dd.o: $(PATHCPU)/eval_host.h $(PATHCPU)/eval_host.cpp \
                $(CDDFLAG)/*.h  path_data_dd.o poly_dd.o workspace_host_dd.o
	$(gpp) -O2 -I $(CDDFLAG) $(IFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/eval_host.cpp -o eval_host_dd.o
	         
path_data_dd.o: $(PATHTEST)/path_data.cpp $(PATHTEST)/path_data.h \
                $(CDDFLAG)/*.h
	$(gpp) -O2 -I $(CDDFLAG)  $(IFLAG) $(QDFLAG) \
               -c $(PATHTEST)/path_data.cpp -o path_data_dd.o

path_host_dd.o: $(PATHCPU)/path_host.h $(PATHCPU)/path_host.cpp \
                predictor_host_dd.o newton_host_dd.o
	$(gpp) -O2  -I $(CDDFLAG) $(IFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/path_host.cpp -o path_host_dd.o

newton_host_dd.o: $(PATHCPU)/newton_host.h $(PATHCPU)/newton_host.cpp \
                  eval_host_dd.o mgs_host_dd.o 
	$(gpp) -O2 -I $(CDDFLAG) $(IFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/newton_host.cpp -o newton_host_dd.o

predictor_host_dd.o: $(PATHCPU)/predictor_host.h \
                     $(PATHCPU)/predictor_host.cpp workspace_host_dd.o
	$(gpp) -O2 -I $(CDDFLAG) $(QDFLAG) \
               -c $(PATHCPU)/predictor_host.cpp -o predictor_host_dd.o

mgs_host_dd.o: $(PATHCPU)/mgs_host.h $(PATHCPU)/mgs_host.cpp 
	$(gpp) -O2 -I $(CDDFLAG) $(GQDFLAG) \
	       -c $(PATHCPU)/mgs_host.cpp -o mgs_host_dd.o

families_dd.o: $(POLY)/families.h $(POLY)/families.cpp
	$(gpp) -O2  -I $(CDDFLAG) $(QDFLAG) \
	       -c $(POLY)/families.cpp -o families_dd.o

poly_dd.o: $(POLY)/poly.h $(POLY)/poly.cpp utilities_dd.o $(CDDFLAG)/*.h
	$(gpp) -O2  -I $(CDDFLAG) $(QDFLAG) \
	       -c $(POLY)/poly.cpp -o poly_dd.o

polysol_dd.o: $(POLY)/poly.h $(POLY)/polysol.cpp utilities_dd.o $(CDFLAG)/*.h
	$(gpp) -O2  -I $(CDFLAG) $(QDFLAG) \
	       -c $(POLY)/polysol.cpp -o polysol_dd.o

path_gpu_dd.o: $(PATHGPU)/*.cu $(PATHGPU)/*.h $(CDDFLAG)/*.h
	@-echo ">>> compiling kernels ..."
	nvcc -O2 --define-macro path_precision=cdd \
	         -I $(CDDFLAG) $(IFLAG)\
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
	         -o path_gpu_dd.o \
             -c $(PATHGPU)/path_gpu.cu --ptxas-options=-v

utilities_dd.o: $(POLY)/utilities.h $(POLY)/utilities.cpp
	$(gpp) -O2  -I $(CDDFLAG)  $(GQDFLAG)\
	       -c $(POLY)/utilities.cpp -o utilities_dd.o 
	         
err_check_dd.o: $(PATHTEST)/err_check.cpp $(PATHTEST)/err_check.h \
                $(CDDFLAG)/*.h
	$(gpp) -O2  -I $(CDDFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/err_check.cpp -o err_check_dd.o
	         
eval_test_dd.o: $(PATHTEST)/eval_test.cpp $(PATHTEST)/err_check.h \
                $(CDDFLAG)/*.h
	$(gpp) -O2  -I $(CDDFLAG)  $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/eval_test.cpp -o eval_test_dd.o
	         
mgs_test_dd.o: $(PATHTEST)/mgs_test.cpp $(PATHTEST)/err_check.h \
               $(CDDFLAG)/*.h  $(PATHGPU)/mgs.cu
	$(gpp) -O2  -I $(CDDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/mgs_test.cpp -o mgs_test_dd.o
	         
newton_test_dd.o: $(PATHTEST)/newton_test.cpp $(PATHTEST)/err_check.h \
                  $(CDDFLAG)/*.h 
	$(gpp) -O2  -I $(CDDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/newton_test.cpp -o newton_test_dd.o
	         
path_test_dd.o: $(PATHTEST)/path_test.cpp $(PATHTEST)/path_test.h \
                $(CDDFLAG)/*.h path_host_dd.o
	$(gpp) -O2  -I $(CDDFLAG)  $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/path_test.cpp -o path_test_dd.o
	         
path_multi_test_dd.o: $(PATHTEST)/path_multi_test.cpp \
                      $(PATHTEST)/err_check.h $(CDDFLAG)/*.h path_test_dd.o
	$(gpp) -O2  -I $(CDDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/path_multi_test.cpp -o path_multi_test_dd.o
	         
predict_test_dd.o: $(PATHTEST)/predict_test.cpp $(PATHTEST)/err_check.h \
                   $(CDDFLAG)/*.h
	$(gpp) -O2  -I $(CDDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/predict_test.cpp -o predict_test_dd.o
	         
pieri_test_dd.o: $(PATHTEST)/pieri_test.cpp $(PATHTEST)/err_check.h \
                 $(CDDFLAG)/*.h
	$(gpp) -O2  -I $(CDDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/pieri_test.cpp -o pieri_test_dd.o
	         
init_test_dd.o: poly_dd.o $(PATHTEST)/init_test.cpp $(PATHTEST)/err_check.h \
                          $(CDDFLAG)/*.h
	$(gpp) -O2  -I $(CDDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/init_test.cpp -o init_test_dd.o
	         
witness_set_test_dd.o: $(PATHTEST)/witness_set_test.cpp \
                       $(PATHTEST)/err_check.h $(CDDFLAG)/*.h
	$(gpp) -O2  -I $(CDDFLAG) $(IFLAG) $(GQDFLAG) \
	       -c $(PATHTEST)/witness_set_test.cpp -o witness_set_test_dd.o
	         
main_test_dd.o: $(PATHTEST)/main_test.cpp $(PATHCPU)/parameter.h \
                families_dd.o path_host_dd.o path_gpu_dd.o \
		newton_host_dd.o err_check_dd.o eval_test_dd.o mgs_test_dd.o \
                newton_test_dd.o predict_test_dd.o \
		path_test_dd.o path_multi_test_dd.o pieri_test_dd.o \
                init_test_dd.o witness_set_test_dd.o
	$(gpp) -Dcompilewgpp=1 -O2 \
               -I$(CDDFLAG) -I$(PHCLib) $(IFLAG) $(GQDFLAG) \
	       -c $(PATHTEST)/main_test.cpp -o main_test_dd.o
	         
path_dd: families_dd.o main_test_dd.o poly_dd.o polysol_dd.o utilities_dd.o \
         path_gpu_dd.o \
         gqd_qd_util.o path_host_dd.o eval_host_dd.o mgs_test_dd.o
	$(gpp) -O2 $(GQDFLAG) main_test_dd.o \
                families_dd.o poly_dd.o polysol_dd.o utilities_dd.o \
                path_gpu_dd.o gqd_qd_util.o predictor_host_dd.o mgs_host_dd.o \
                newton_host_dd.o path_host_dd.o path_data_dd.o eval_host_dd.o \
		workspace_host_dd.o err_check_dd.o mgs_test_dd.o \
                predict_test_dd.o  eval_test_dd.o path_test_dd.o \
                path_multi_test_dd.o newton_test_dd.o pieri_test_dd.o \
                init_test_dd.o witness_set_test_dd.o \
		$(QD_LIB)/libqd.a -o $(BIN)/path_dd \
		-lcutil_x86_64 -lcudart \
		-L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

# QD

workspace_host_qd.o: $(PATHCPU)/workspace_host.h $(PATHCPU)/workspace_host.cpp
	$(gpp) -O2 -I./$(CQDFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/workspace_host.cpp -o workspace_host_qd.o

eval_host_qd.o: $(PATHCPU)/eval_host.h $(PATHCPU)/eval_host.cpp \
                $(CQDFLAG)/*.h path_data_qd.o poly_qd.o workspace_host_qd.o
	$(gpp) -O2 -I./$(CQDFLAG) $(IFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/eval_host.cpp -o eval_host_qd.o
	         
path_data_qd.o: $(PATHTEST)/path_data.cpp $(PATHTEST)/path_data.h \
                $(CQDFLAG)/*.h
	$(gpp) -O2 -I./$(CQDFLAG)  $(IFLAG) $(QDFLAG) \
               -c $(PATHTEST)/path_data.cpp -o path_data_qd.o

path_host_qd.o: $(PATHCPU)/path_host.h $(PATHCPU)/path_host.cpp \
                predictor_host_qd.o newton_host_qd.o
	$(gpp) -O2  -I./$(CQDFLAG) $(IFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/path_host.cpp -o path_host_qd.o

newton_host_qd.o: $(PATHCPU)/newton_host.h $(PATHCPU)/newton_host.cpp \
                  eval_host_qd.o mgs_host_qd.o 
	$(gpp) -O2 -I./$(CQDFLAG) $(IFLAG) $(QDFLAG) \
	       -c $(PATHCPU)/newton_host.cpp -o newton_host_qd.o

predictor_host_qd.o: $(PATHCPU)/predictor_host.h \
                     $(PATHCPU)/predictor_host.cpp workspace_host_qd.o
	$(gpp) -O2 -I./$(CQDFLAG) $(QDFLAG) \
               -c $(PATHCPU)/predictor_host.cpp -o predictor_host_qd.o

mgs_host_qd.o: $(PATHCPU)/mgs_host.h $(PATHCPU)/mgs_host.cpp 
	$(gpp) -O2 -I./$(CQDFLAG) $(GQDFLAG) \
	       -c $(PATHCPU)/mgs_host.cpp -o mgs_host_qd.o

families_qd.o: $(POLY)/families.h $(POLY)/families.cpp
	$(gpp) -O2  -I./$(CQDFLAG) $(QDFLAG) \
	       -c $(POLY)/families.cpp -o families_qd.o

poly_qd.o: $(POLY)/poly.h $(POLY)/poly.cpp utilities_qd.o $(CQDFLAG)/*.h
	$(gpp) -O2  -I./$(CQDFLAG) $(QDFLAG) \
	       -c $(POLY)/poly.cpp -o poly_qd.o

polysol_qd.o: $(POLY)/poly.h $(POLY)/polysol.cpp utilities_qd.o $(CDFLAG)/*.h
	$(gpp) -O2  -I $(CDFLAG) $(QDFLAG) \
	       -c $(POLY)/polysol.cpp -o polysol_qd.o

path_gpu_qd.o: $(PATHGPU)/*.cu $(PATHGPU)/*.h $(CQDFLAG)/*.h
	@-echo ">>> compiling kernels ..."
	nvcc -O2 --define-macro path_precision=cqd \
	         -I./$(CQDFLAG) $(IFLAG)\
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
	         -o path_gpu_qd.o \
             -c $(PATHGPU)/path_gpu.cu --ptxas-options=-v

utilities_qd.o: $(POLY)/utilities.h $(POLY)/utilities.cpp
	$(gpp) -O2  -I./$(CQDFLAG)  $(GQDFLAG)\
	       -c $(POLY)/utilities.cpp -o utilities_qd.o 
	         
err_check_qd.o: $(PATHTEST)/err_check.cpp $(PATHTEST)/err_check.h \
                $(CQDFLAG)/*.h
	$(gpp) -O2  -I./$(CQDFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/err_check.cpp -o err_check_qd.o
	         
eval_test_qd.o: $(PATHTEST)/eval_test.cpp $(PATHTEST)/err_check.h \
                $(CQDFLAG)/*.h
	$(gpp) -O2  -I./$(CQDFLAG)  $(IFLAG) $(GQDFLAG) \
	       -c $(PATHTEST)/eval_test.cpp -o eval_test_qd.o
	         
mgs_test_qd.o: $(PATHTEST)/mgs_test.cpp $(PATHTEST)/err_check.h \
               $(CQDFLAG)/*.h $(PATHGPU)/mgs.cu
	$(gpp) -O2  -I./$(CQDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/mgs_test.cpp -o mgs_test_qd.o
	         
newton_test_qd.o: $(PATHTEST)/newton_test.cpp $(PATHTEST)/err_check.h \
                  $(CQDFLAG)/*.h 
	$(gpp) -O2  -I./$(CQDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/newton_test.cpp -o newton_test_qd.o
	         
path_test_qd.o: $(PATHTEST)/path_test.cpp $(PATHTEST)/path_test.h \
                $(CQDFLAG)/*.h path_host_qd.o
	$(gpp) -O2  -I./$(CQDFLAG)  $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/path_test.cpp -o path_test_qd.o
	         
path_multi_test_qd.o: $(PATHTEST)/path_multi_test.cpp \
                      $(PATHTEST)/err_check.h $(CQDFLAG)/*.h path_test_qd.o
	$(gpp) -O2  -I./$(CQDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/path_multi_test.cpp -o path_multi_test_qd.o
	         
predict_test_qd.o: $(PATHTEST)/predict_test.cpp $(PATHTEST)/err_check.h \
                   $(CQDFLAG)/*.h
	$(gpp) -O2  -I./$(CQDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/predict_test.cpp -o predict_test_qd.o
	         
pieri_test_qd.o: $(PATHTEST)/pieri_test.cpp $(PATHTEST)/err_check.h \
                 $(CQDFLAG)/*.h
	$(gpp) -O2  -I./$(CQDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/pieri_test.cpp -o pieri_test_qd.o
	         
init_test_qd.o: poly_qd.o $(PATHTEST)/init_test.cpp $(PATHTEST)/err_check.h \
                          $(CQDFLAG)/*.h
	$(gpp) -O2  -I./$(CQDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/init_test.cpp -o init_test_qd.o
	         
witness_set_test_qd.o: $(PATHTEST)/witness_set_test.cpp \
                       $(PATHTEST)/err_check.h $(CQDFLAG)/*.h
	$(gpp) -O2  -I./$(CQDFLAG) $(IFLAG) $(GQDFLAG)\
	       -c $(PATHTEST)/witness_set_test.cpp -o witness_set_test_qd.o
	         
main_test_qd.o: $(PATHTEST)/main_test.cpp $(PATHCPU)/parameter.h \
                families_qd.o path_host_qd.o path_gpu_qd.o \
		newton_host_qd.o err_check_qd.o eval_test_qd.o \
                mgs_test_qd.o newton_test_qd.o predict_test_qd.o \
		path_test_qd.o path_multi_test_qd.o pieri_test_qd.o \
                init_test_qd.o witness_set_test_qd.o
	$(gpp) -Dcompilewgpp=1 -O2 \
               -I./$(CQDFLAG) -I$(PHCLib) $(IFLAG) $(GQDFLAG) \
	       -c $(PATHTEST)/main_test.cpp -o main_test_qd.o
	         
path_qd: families_qd.o main_test_qd.o poly_qd.o polysol_qd.o utilities_qd.o \
         path_gpu_qd.o \
         gqd_qd_util.o path_host_qd.o eval_host_qd.o mgs_test_qd.o
	$(gpp) -O2 $(GQDFLAG) main_test_qd.o \
                families_qd.o poly_qd.o polysol_qd.o utilities_qd.o \
                path_gpu_qd.o gqd_qd_util.o predictor_host_qd.o mgs_host_qd.o \
                newton_host_qd.o path_host_qd.o path_data_qd.o eval_host_qd.o \
		workspace_host_qd.o err_check_qd.o mgs_test_qd.o \
                predict_test_qd.o  eval_test_qd.o path_test_qd.o \
                path_multi_test_qd.o newton_test_qd.o pieri_test_qd.o \
                init_test_qd.o witness_set_test_qd.o \
		$(QD_LIB)/libqd.a -o $(BIN)/path_qd \
		-lcutil_x86_64 -lcudart \
		-L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib
		
# ada interface (stil under construction)

lib_syscon: poly_d.o
	gnatmake -O2 -c $(INCLULIBS) $(GNATFLAGS) use_c2phc.adb
	g++ -Dcompilewgpp=1 -O2 -c $(PHCLib)/syscon.c
	g++ -Dcompilewgpp=1 -O2 -c $(PHCLib)/solcon.c
	g++ -Dcompilewgpp=1 -O2 -c $(PHCLib)/phcpack.c
	g++ -Dcompilewgpp=1 -O2 -I$(PHCLib) -I$(POLY) -I$(CDFLAG) $(QDFLAG) \
            -c $(TEST)/lib_syscon.cpp
	gnatbind -n use_c2phc
	gnatlink -O2 use_c2phc phcpack.o solcon.o syscon.o \
                     lib_syscon.o poly_d.o utilities_d.o \
                 -o $(BIN)/test --LINK=g++
        
main_test_ada.o: $(PATHTEST)/main_test.cpp $(PATHCPU)/parameter.h \
                 families_d.o path_host_d.o path_gpu_d.o newton_host_d.o \
                 err_check_d.o eval_test_d.o mgs_test_d.o newton_test_d.o \
                 predict_test_d.o path_test_d.o path_multi_test_d.o \
                 pieri_test_d.o init_test_d.o witness_set_test_d.o
	g++ -O2 -I$(CDFLAG) $(IFLAG) $(GQDFLAG) -I$(PHCLib) \
                -Dcompilewgpp=1 \
	        -c $(PATHTEST)/main_test.cpp -o main_test_ada.o
	         
path_ada: main_test_ada.o families_d.o poly_d.o \
          utilities_d.o path_gpu_d.o gqd_qd_util.o path_host_d.o \
          eval_host_d.o mgs_test_d.o	
	gnatmake -O2 -c $(INCLULIBS) $(GNATFLAGS) use_c2phc.adb
	g++ -Dcompilewgpp=1 -O2 -c $(PHCLib)/syscon.c
	g++ -Dcompilewgpp=1 -O2 -c $(PHCLib)/solcon.c
	g++ -Dcompilewgpp=1 -O2 -c $(PHCLib)/phcpack.c
	g++ -Dcompilewgpp=1 -O2 -c $(PHCLib)/jump_track.c
	# gcc -c $(PHCLib)/lib_syscon.c
	g++ -Dcompilewgpp=1 -O2 -I$(PHCLib) -I$(POLY) -I$(CDFLAG) $(QDFLAG) \
                               -c $(PATHTEST)/ada_test.cpp
	gnatbind -n use_c2phc
	gnatlink -O2 $(GQDFLAG) -I$(PHCLib) \
		 -C use_c2phc phcpack.o jump_track.o solcon.o syscon.o \
                    ada_test.o main_test_ada.o families_d.o poly_d.o \
                    utilities_d.o path_gpu_d.o gqd_qd_util.o \
                    predictor_host_d.o mgs_host_d.o newton_host_d.o \
                    path_host_d.o path_data_d.o eval_host_d.o \
                    workspace_host_d.o err_check_d.o mgs_test_d.o \
                    predict_test_d.o  eval_test_d.o path_test_d.o \
                    path_multi_test_d.o newton_test_d.o pieri_test_d.o \
                    init_test_d.o witness_set_test_d.o \
		$(QD_LIB)/libqd.a -o $(BIN)/path_ada \
		-lcutil_x86_64 -lcudart \
		-L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib \
		--LINK=g++
             
clean_path:
	rm -rf *.o
	rm -rf $(BIN)/path_d $(BIN)/path_dd $(BIN)/path_qd
             
clean:
	rm -rf GNAT-* *.o *.ali *.adb *.ads test 
	rm -rf $(BIN)/*
