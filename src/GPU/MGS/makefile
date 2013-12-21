QD_ROOT=/usr/local/qd-2.3.9
GQD_HOME=/usr/local/gqd_1_2
SDK_HOME=/usr/local/cuda/sdk

all: run_mgs_d run_mgs_dd run_mgs_qd

run_mgs_d:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=d -I./DefineTypesD \
             -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c mgs_kernelsT.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling a separate library ..."
	g++ -I./DefineTypesD \
            -Wall -c separate.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o separate.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs ..."
	g++ -I./DefineTypesD -c run_mgs.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs_kernelsT.o run_mgs.o gqd_qd_util.o separate.o \
		$(QD_ROOT)/src/libqd.a -o run_mgs_d \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs_dd1:
	@-echo ">>> compiling kernels ..."
	nvcc -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c mgs_kernelsT_dd.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling a separate library ..."
	g++ -Wall -c separate.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o separate.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs ..."
	g++ -c run_mgs.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs_kernelsT_dd.o run_mgs.o gqd_qd_util.o separate.o \
		$(QD_ROOT)/src/libqd.a -o run_mgs_dd1 \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs_dd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=dd -I./DefineTypesDD \
	     -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c mgs_kernelsT.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling a separate library ..."
	g++ -Wall -c separate.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -I./DefineTypesDD -o separate.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs ..."
	g++ -c run_mgs.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -I./DefineTypesDD
	@-echo ">>> linking ..."
	g++ -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs_kernelsT.o run_mgs.o gqd_qd_util.o separate.o \
		$(QD_ROOT)/src/libqd.a -I./DefineTypesDD -o run_mgs_dd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

run_mgs_qd:
	@-echo ">>> compiling kernels ..."
	nvcc --define-macro precision=qd -I./DefineTypesQD \
	     -I$(GQD_HOME)/inc -I$(SDK_HOME)/C/common/inc \
             -I/usr/local/cuda/include \
             -c mgs_kernelsT_qd.cu -arch=sm_13 --ptxas-options=-v
	@-echo "compiling a separate library ..."
	g++ -I./DefineTypesQD \
            -Wall -c separate.cpp -I/usr/local/gqd_1_2/inc \
            -I/usr/local/cuda/include -I/usr/local/include \
            -o separate.o
	@-echo ">>> compiling utilities ..."
	g++ -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include -c gqd_qd_util.cpp
	@-echo ">>> compiling run_mgs ..."
	g++ -I./DefineTypesQD -c run_mgs.cpp \
            -I/usr/local/cuda/include -I$(GQD_HOME)/inc \
            -I$(QD_ROOT)/include
	@-echo ">>> linking ..."
	g++ -I./DefineTypesQD -I$(GQD_HOME)/inc -I$(QD_ROOT)/include \
		mgs_kernelsT_qd.o run_mgs.o gqd_qd_util.o separate.o \
		$(QD_ROOT)/src/libqd.a -o run_mgs_qd \
            -lcuda -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

clean:
	test -z "run_mgs_d" || rm -f run_mgs_d
	test -z "run_mgs_dd" || rm -f run_mgs_dd
	test -z "run_mgs_qd" || rm -f run_mgs_qd
	rm -f *.o
