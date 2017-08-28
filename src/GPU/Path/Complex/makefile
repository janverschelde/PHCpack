# Location of the QD library and its GPU version:

QD_ROOT=/usr/local/qd-2.3.17
QD_LIB=/usr/local/lib
GQD_HOME=/usr/local/gqd_1_2
SDK_HOME=/usr/local/cuda/sdk

QDFLAG=-I$(QD_ROOT)/include
GQDFLAG=-I$(QD_ROOT)/include -I/usr/local/cuda/include -I$(GQD_HOME)/inc
SDKFLAG=-I$(SDK_HOME)/C/common/inc

gpp=/usr/bin/g++

test_complexH:
	$(gpp) -o /tmp/test_complexH -I$(QD_ROOT)/include test_complexH.cpp \
                  /usr/local/lib/libqd.a

test_complexsqrt_d:
	nvcc -o /tmp/test_complexsqrt_d $(GQDFLAG) test_complexsqrt_d.cu

test_complexdsqrt:
	nvcc -ccbin=$(gpp) $(GQDFLAG) $(SDKFLAG) -c test_complexdsqrt.cu
	$(gpp) test_complexdsqrt.o $(QD_LIB)/libqd.a -o /tmp/test_complexdsqrt \
            -lcutil_x86_64 -lcudart \
            -L/usr/local/cuda/lib64 -L$(SDK_HOME)/C/lib

test_complexdsqrt_dd:
	nvcc $(GQDFLAG) $(SDKFLAG) -c test_complexdsqrt_dd.cu

gqd_qd_util.o:
	$(gpp) -O2 $(GQDFLAG) -c gqd_qd_util.cpp 

gqd2qd_convertors.o:
	$(gpp) $(GQDFLAG) -c gqd2qd_convertors.cpp

clean:
	/bin/rm -f -r *.o
