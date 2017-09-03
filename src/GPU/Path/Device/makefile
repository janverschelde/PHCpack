# Location of the QD library and its GPU version:

QD_ROOT=/usr/local/qd-2.3.17
QD_LIB=/usr/local/lib
GQD_HOME=/usr/local/gqd_1_2
SDK_HOME=/usr/local/cuda/sdk

QDFLAG=-I$(QD_ROOT)/include
GQDFLAG=-I$(QD_ROOT)/include -I/usr/local/cuda/include -I$(GQD_HOME)/inc
SDKFLAG=-I$(SDK_HOME)/C/common/inc

gpp=/usr/bin/g++

mgs_large_norm.o:
	nvcc $(GQDFLAG) $(SDKFLAG) -c mgs_large_norm.cu

mgs_large_row_reduce.o:
	nvcc $(GQDFLAG) $(SDKFLAG) -c mgs_large_row_reduce.cu

mgs_large_backsub.o:
	nvcc $(GQDFLAG) $(SDKFLAG) -c mgs_large_backsub.cu

mgs_large_row.o:
	nvcc $(GQDFLAG) $(SDKFLAG) -c mgs_large_row.cu

mgs_large_block_reduce.o:
	nvcc $(GQDFLAG) $(SDKFLAG) -c mgs_large_block_reduce.cu

mgs_large_block.o:
	nvcc $(GQDFLAG) $(SDKFLAG) -c mgs_large_block.cu

mgs.o:
	nvcc $(GQDFLAG) $(SDKFLAG) -c mgs.cu

clean:
	/bin/rm -f -r *.o
