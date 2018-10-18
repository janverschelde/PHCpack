QD_ROOT=/usr/local/qd-2.3.17
QD_LIB=/usr/local/lib

gpp=/usr/bin/g++

CMPLX = ../Complex
POLYS = ../Poly3

test_mgs:
	$(gpp) -c -I$(CMPLX) test_mgs.cpp -o test_mgs.o
	$(gpp) test_mgs.o $(QD_LIB)/libqd.a -o /tmp/test_mgs

parameter.o:
	$(gpp) -c parameter.cpp -o parameter.o

test_parameter: parameter.o
	$(gpp) -c -I$(CMPLX) test_parameter.cpp -o test_parameter.o
	$(gpp) test_parameter.o parameter.o $(QD_LIB)/libqd.a \
               -o /tmp/test_parameter

test_pathdata:
	$(gpp) -c -I$(CMPLX) -I$(POLYS) test_pathdata.cpp -o test_pathdata.o
	$(gpp) test_pathdata.o $(QD_LIB)/libqd.a -o /tmp/test_pathdata

test_workspace:
	$(gpp) -c -I$(CMPLX) test_workspace.cpp -o test_workspace.o
	$(gpp) test_workspace.o $(QD_LIB)/libqd.a -o /tmp/test_workspace

test_eval:
	$(gpp) -c -I$(CMPLX) -I$(POLYS) test_eval.cpp -o test_eval.o
	$(gpp) test_eval.o $(QD_LIB)/libqd.a -o /tmp/test_eval

test_newton: parameter.o
	$(gpp) -c -I$(CMPLX) -I$(POLYS) test_newton.cpp -o test_newton.o
	$(gpp) parameter.o test_newton.o $(QD_LIB)/libqd.a -o /tmp/test_newton

test_jobs:
	$(gpp) -pthread -o /tmp/test_jobs test_jobs.cpp -lpthread

test_worker:
	$(gpp) -pthread -o /tmp/test_worker test_worker.cpp -lpthread

test_start:
	$(gpp) -c -I$(CMPLX) -I$(POLYS) test_start.cpp -o test_start.o
	$(gpp) test_start.o $(QD_LIB)/libqd.a -o /tmp/test_start

test_path: parameter.o
	$(gpp) -pthread -c -I$(CMPLX) -I$(POLYS) test_path.cpp -o test_path.o
	$(gpp) test_path.o parameter.o $(QD_LIB)/libqd.a -lpthread \
               -o /tmp/test_path

clean:
	/bin/rm -f -r *.o
