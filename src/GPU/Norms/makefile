# tested with nvcc version 5.5.0 and gcc 4.4.7 on C2050
# Note that the -arch=sm_13 is needed on the C2050.
# tested with nvcc version 10.1 and gcc 4.8.5 on Volta

gpp=/usr/bin/g++
gcc=/usr/bin/gcc

dbl_norm_kernels.o:
	@-echo ">>> compiling vector norm kernels for doubles ..."
	# nvcc -ccbin=$(gcc) -arch=sm_13 -c dbl_norm_kernels.cu
	nvcc -ccbin=$(gcc) -c dbl_norm_kernels.cu

dbl_norm_host.o:
	@-echo ">>> compile code for execution on the host ..."
	$(gpp) -O3 -c dbl_norm_host.cpp
	
run_norm_d.o:
	@-echo ">>> compiling the main program ..."
	$(gpp) -I/usr/local/cuda/include -O3 -c run_norm_d.cpp

run_norm_d: dbl_norm_kernels.o dbl_norm_host.o run_norm_d.o
	@-echo ">>> linking ..."
	$(gpp) -o run_norm_d \
               dbl_norm_kernels.o dbl_norm_host.o run_norm_d.o \
            -lcuda -lcudart -L/usr/local/cuda/lib64 

clean:
	/bin/rm -f -r run_norm_d
	/bin/rm -f -r dbl_norm_kernels.o dbl_norm_host.o run_norm_d.o
