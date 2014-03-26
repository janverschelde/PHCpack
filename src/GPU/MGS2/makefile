QD_ROOT=/usr/local/qd-2.3.13
QD_LIB=/usr/local/lib
GQD_HOME=/usr/local/gqd_1_2
SDK_HOME=/usr/local/cuda/sdk

norm: run_norm_d run_norm_dd run_norm_qd
gram: run_gram_d run_gram_dd run_gram_qd
mgs2: run_mgs2_d run_mgs2_dd run_mgs2_qd
mgs3: run_mgs3_d run_mgs3_dd run_mgs3_qd
mgs: run_mgs_d run_mgs_dd run_mgs_qd
chan: chandra_d chandra_dd chandra_qd
all: norm gram mgs3 mgs2 mgs chan

chandra_d:
	@-echo ">>> compiling evaluation and differentiation routines ..."
	g++ -I./DefineTypesD -c chandra.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o /tmp/chandra_d.o
	@-echo ">>> compiling the main test program ..."
	g++ -I./DefineTypesD -c chandra_test.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o /tmp/chandra_test_d.o
	@-echo ">>> linking ..."
	g++ /tmp/chandra_d.o /tmp/chandra_test_d.o -o /tmp/chandra_d

chandra_dd:
	@-echo ">>> compiling evaluation and differentiation routines ..."
	g++ -I./DefineTypesDD -c chandra.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o /tmp/chandra_dd.o
	@-echo ">>> compiling the main test program ..."
	g++ -I./DefineTypesDD -c chandra_test.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o /tmp/chandra_test_dd.o
	@-echo ">>> linking ..."
	g++ /tmp/chandra_dd.o /tmp/chandra_test_dd.o \
            $(QD_LIB)/libqd.a -o /tmp/chandra_dd

chandra_qd:
	@-echo ">>> compiling evaluation and differentiation routines ..."
	g++ -I./DefineTypesQD -c chandra.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o /tmp/chandra_qd.o
	@-echo ">>> compiling the main test program ..."
	g++ -I./DefineTypesQD -c chandra_test.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o /tmp/chandra_test_qd.o
	@-echo ">>> linking ..."
	g++ /tmp/chandra_qd.o /tmp/chandra_test_qd.o \
            $(QD_LIB)/libqd.a -o /tmp/chandra_qd

run_mgs2_d:
	@-echo ">>> compiling kernels ..."
	nvcc -O2 --define-macro precision=d -I./DefineTypesD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include -o mgs2_kernels_d.o \
             -c mgs2_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -O2 -I./DefineTypesD \
            -Wall -c mgs2_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o mgs2_host_d.o
	@-echo ">>> compiling utilities ..."
	g++ -O2 -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs2 ..."
	g++ -O2 -I$(GQD_HOME)/inc -I$(QD_ROOT)/include -I./DefineTypesD \
            -I/usr/local/cuda/include -c run_mgs2.cpp -o run_mgs2_d.o
	@-echo ">>> compiling chandra evaluation and differentiation ..."
	g++ -O2 -I./DefineTypesD -c chandra.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o chandra_d.o
	@-echo ">>> linking ..."
	g++ -O2 -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs2_kernels_d.o mgs2_host_d.o gqd_qd_util.o chandra_d.o \
                run_mgs2_d.o $(QD_LIB)/libqd.a -o run_mgs2_d \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs3_d:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=d -I./DefineTypesD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include -o mgs3_kernels_d.o \
             -c mgs3_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesD \
            -Wall -c mgs2_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o mgs2_host_d.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs3 ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include -I./DefineTypesD \
            -I/usr/local/cuda/include -c run_mgs3.cpp -o run_mgs3_d.o
	@-echo ">>> compiling chandra evaluation and differentiation ..."
	g++ -I./DefineTypesD -c chandra.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o chandra_d.o
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs3_kernels_d.o mgs2_host_d.o gqd_qd_util.o chandra_d.o \
                run_mgs3_d.o $(QD_LIB)/libqd.a -o run_mgs3_d \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs2_dd:
	@-echo ">>> compiling kernels ..."
	nvcc -O2 --define-macro precision=dd -I./DefineTypesDD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include -o mgs2_kernels_dd.o \
             -c mgs2_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -O2 -I./DefineTypesDD \
            -Wall -c mgs2_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o mgs2_host_dd.o
	@-echo ">>> compiling utilities ..."
	g++ -O2 -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs2 ..."
	g++ -O2 -I$(GQD_HOME)/inc -I$(QD_ROOT)/include -I./DefineTypesDD \
            -I/usr/local/cuda/include -c run_mgs2.cpp -o run_mgs2_dd.o
	@-echo ">>> compiling chandra evaluation and differentiation ..."
	g++ -O2 -I./DefineTypesDD -c chandra.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o chandra_dd.o
	@-echo ">>> linking ..."
	g++ -O2 -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs2_kernels_dd.o mgs2_host_dd.o gqd_qd_util.o chandra_dd.o \
                run_mgs2_dd.o $(QD_LIB)/libqd.a -o run_mgs2_dd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs3_dd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=dd -I./DefineTypesDD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include -o mgs3_kernels_dd.o \
             -c mgs3_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesDD \
            -Wall -c mgs2_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o mgs2_host_dd.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs2 ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include -I./DefineTypesDD \
            -I/usr/local/cuda/include -c run_mgs3.cpp -o run_mgs3_dd.o
	@-echo ">>> compiling chandra evaluation and differentiation ..."
	g++ -I./DefineTypesDD -c chandra.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o chandra_dd.o
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs3_kernels_dd.o mgs2_host_dd.o gqd_qd_util.o chandra_dd.o \
                run_mgs3_dd.o $(QD_LIB)/libqd.a -o run_mgs3_dd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs2_qd:
	@-echo ">>> compiling kernels ..."
	nvcc -O2 --define-macro precision=qd -I./DefineTypesQD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include -o mgs2_kernels_qd.o \
             -c mgs2_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -O2 -I./DefineTypesQD \
            -Wall -c mgs2_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o mgs2_host_qd.o
	@-echo ">>> compiling utilities ..."
	g++ -O2 -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs2 ..."
	g++ -O2 -I$(GQD_HOME)/inc -I$(QD_ROOT)/include -I./DefineTypesQD \
            -I/usr/local/cuda/include -c run_mgs2.cpp -o run_mgs2_qd.o
	@-echo ">>> compiling chandra evaluation and differentiation ..."
	g++ -O2 -I./DefineTypesQD -c chandra.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o chandra_qd.o
	@-echo ">>> linking ..."
	g++ -O2 -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs2_kernels_qd.o mgs2_host_qd.o gqd_qd_util.o chandra_qd.o \
                run_mgs2_qd.o $(QD_LIB)/libqd.a -o run_mgs2_qd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs3_qd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=qd -I./DefineTypesQD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include -o mgs3_kernels_qd.o \
             -c mgs3_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesQD \
            -Wall -c mgs2_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o mgs2_host_qd.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs2 ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include -I./DefineTypesQD \
            -I/usr/local/cuda/include -c run_mgs3.cpp -o run_mgs3_qd.o
	@-echo ">>> compiling chandra evaluation and differentiation ..."
	g++ -I./DefineTypesQD -c chandra.cpp -I/usr/local/cuda/include \
            -I$(QD_ROOT)/include -I$(GQD_HOME)/inc -o chandra_qd.o
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs3_kernels_qd.o mgs2_host_qd.o gqd_qd_util.o chandra_qd.o \
                run_mgs3_qd.o $(QD_LIB)/libqd.a -o run_mgs3_qd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_gram_d:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=d -I./DefineTypesD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c gram_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesD \
            -Wall -c gram_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o gram_host.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_gram ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include -I./DefineTypesD \
            -I/usr/local/cuda/include -c run_gram.cpp
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		gram_kernels.o gram_host.o gqd_qd_util.o run_gram.o \
		$(QD_LIB)/libqd.a -o run_gram_d \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_gram_dd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=dd -I./DefineTypesDD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c gram_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesDD \
            -Wall -c gram_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o gram_host.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_gram ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include -I./DefineTypesDD \
            -I/usr/local/cuda/include -c run_gram.cpp
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		gram_kernels.o gram_host.o gqd_qd_util.o run_gram.o \
		$(QD_LIB)/libqd.a -o run_gram_dd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_gram_qd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=qd -I./DefineTypesQD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c gram_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesQD \
            -Wall -c gram_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o gram_host.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_gram ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include -I./DefineTypesQD \
            -I/usr/local/cuda/include -c run_gram.cpp
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		gram_kernels.o gram_host.o gqd_qd_util.o run_gram.o \
		$(QD_LIB)/libqd.a -o run_gram_qd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_norm_d:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=d -I./DefineTypesD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c norm_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesD \
            -Wall -c norm_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o norm_host.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_norm ..."
	g++ -I./DefineTypesD -c run_norm.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		norm_kernels.o run_norm.o gqd_qd_util.o norm_host.o \
		$(QD_LIB)/libqd.a -o run_norm_d \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_norm_dd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=dd -I./DefineTypesDD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c norm_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesDD \
            -Wall -c norm_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o norm_host.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_norm ..."
	g++ -I./DefineTypesDD -c run_norm.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		norm_kernels.o run_norm.o gqd_qd_util.o norm_host.o \
		$(QD_LIB)/libqd.a -o run_norm_dd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_norm_qd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=qd -I./DefineTypesQD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c norm_kernels.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesQD \
            -Wall -c norm_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o norm_host.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_norm ..."
	g++ -I./DefineTypesQD -c run_norm.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		norm_kernels.o run_norm.o gqd_qd_util.o norm_host.o \
		$(QD_LIB)/libqd.a -o run_norm_qd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs_d:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=d -I./DefineTypesD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c mgs_kernelsT.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesD \
            -Wall -c mgs_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o mgs_host.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs ..."
	g++ -I./DefineTypesD -c run_mgs.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs_kernelsT.o run_mgs.o gqd_qd_util.o mgs_host.o \
		$(QD_LIB)/libqd.a -o run_mgs_d \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs_dd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=dd -I./DefineTypesDD \
	     -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c mgs_kernelsT.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -Wall -c mgs_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -I./DefineTypesDD -o mgs_host.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs ..."
	g++ -c run_mgs.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -I./DefineTypesDD
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs_kernelsT.o run_mgs.o gqd_qd_util.o mgs_host.o \
		$(QD_LIB)/libqd.a -I./DefineTypesDD -o run_mgs_dd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs_qd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=qd -I./DefineTypesQD \
	     -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c mgs_kernelsT_qd.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling code for execution on the host ..."
	g++ -I./DefineTypesQD \
            -Wall -c mgs_host.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o mgs_host.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs ..."
	g++ -I./DefineTypesQD -c run_mgs.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I./DefineTypesQD -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs_kernelsT_qd.o run_mgs.o gqd_qd_util.o mgs_host.o \
		$(QD_LIB)/libqd.a -o run_mgs_qd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

clean_norm:
	test -z "run_norm_d" || rm -f run_norm_d
	test -z "run_norm_dd" || rm -f run_norm_dd
	test -z "run_norm_qd" || rm -f run_norm_qd
	rm -f *.o

clean_gram:
	test -z "run_gram_d" || rm -f run_gram_d
	test -z "run_gram_dd" || rm -f run_gram_dd
	test -z "run_gram_qd" || rm -f run_gram_qd
	rm -f *.o

clean_mgs3:
	test -z "run_mgs3_d" || rm -f run_mgs3_d
	test -z "run_mgs3_dd" || rm -f run_mgs3_dd
	test -z "run_mgs3_qd" || rm -f run_mgs3_qd
	rm -f *.o

clean_mgs2:
	test -z "run_mgs2_d" || rm -f run_mgs2_d
	test -z "run_mgs2_dd" || rm -f run_mgs2_dd
	test -z "run_mgs2_qd" || rm -f run_mgs2_qd
	rm -f *.o

clean_mgs:
	test -z "run_mgs_d" || rm -f run_mgs_d
	test -z "run_mgs_dd" || rm -f run_mgs_dd
	test -z "run_mgs_qd" || rm -f run_mgs_qd
	rm -f *.o

clean_chan:
	test -z "chandra_d" || rm -f chandra_d
	test -z "chandra_dd" || rm -f chandra_dd
	test -z "chandra_qd" || rm -f chandra_qd
	rm -f *.o

clean: clean_norm clean_gram clean_mgs3 clean_mgs2 clean_mgs clean_chan
