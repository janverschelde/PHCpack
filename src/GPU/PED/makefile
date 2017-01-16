# This makefile works for the K20C and C2050.
# For the P100, remove the -arch=sm_13 option.

QD_ROOT=/usr/local/qd-2.3.17
QD_LIB=/usr/local/lib
GQD_HOME=/usr/local/gqd_1_2
SDK_HOME=/usr/local/cuda/sdk

all: run_ped_d run_ped_dd run_ped_qd test_ped

use_complex:
	g++ -o /tmp/use_complex -I$(QD_ROOT)/include use_complex.cpp \
            /usr/local/lib/libqd.a

test_ped:
	@-echo ">>> compiling the kernels ..."
	nvcc -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c all_ped_kernels.cu -arch=sm_13 -o ped_kernels.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling ped_host for double precision ..."
	g++ --define-macro precision=d -c ped_host.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -o ped_host_d.o
	@-echo ">>> compiling ped_host for double double precision ..."
	g++ --define-macro precision=dd -c ped_host.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -o ped_host_dd.o
	@-echo ">>> compiling ped_host for quad double precision ..."
	g++ --define-macro precision=qd -c ped_host.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -o ped_host_qd.o
	@-echo ">>> compiling test_ped ..."
	g++ -c test_ped.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include test_ped.o \
		gqd_qd_util.o ped_host_d.o ped_host_dd.o ped_host_qd.o \
		ped_kernels.o $(QD_LIB)/libqd.a -o test_ped \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_ped_d:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=d -I./DefineTypesD \
	     -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c ped_kernelsT.cu -arch=sm_13
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_ped ..."
	g++ --define-macro precision=d -I./DefineTypesD -c run_ped.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		ped_kernelsT.o run_ped.o gqd_qd_util.o \
		$(QD_LIB)/libqd.a -o run_ped_d \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_ped_dd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=dd -I./DefineTypesDD \
	     -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c ped_kernelsT.cu -arch=sm_13 --ptxas-options=-v
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_ped ..."
	g++ --define-macro precision=dd -I./DefineTypesDD -c run_ped.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		ped_kernelsT.o run_ped.o gqd_qd_util.o \
		$(QD_LIB)/libqd.a -I./DefineTypesDD -o run_ped_dd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_ped_qd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=qd -I./DefineTypesQD \
	     -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c ped_kernelsT_qd.cu -arch=sm_13 --ptxas-options=-v
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_ped ..."
	g++ --define-macro precision=qd -I./DefineTypesQD -c run_ped.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I./DefineTypesQD -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		ped_kernelsT_qd.o run_ped.o gqd_qd_util.o \
		$(QD_LIB)/libqd.a -o run_ped_qd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

clean:
	test -z "run_ped_d" || rm -f run_ped_d
	test -z "run_ped_dd" || rm -f run_ped_dd
	test -z "run_ped_qd" || rm -f run_ped_qd
	test -z "test_ped" || rm -f test_ped
	rm -f *.o
