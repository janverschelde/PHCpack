# Location of the QD library and its GPU version:

QD_ROOT=/usr/local/qd-2.3.17
QD_LIB=/usr/local/lib
GQD_HOME=/usr/local/gqd_1_2
SDK_HOME=/usr/local/cuda/sdk

QDFLAG=-I$(QD_ROOT)/include
GQDFLAG=-I$(QD_ROOT)/include -I/usr/local/cuda/include -I$(GQD_HOME)/inc

gpp=/usr/bin/g++

test_complexH:
	$(gpp) -o /tmp/test_complexH -I$(QD_ROOT)/include test_complexH.cpp \
                  /usr/local/lib/libqd.a

test_complexsqrt_d:
	nvcc -o /tmp/test_complexsqrt_d -I$(GQDFLAG) test_complexsqrt_d.cu
