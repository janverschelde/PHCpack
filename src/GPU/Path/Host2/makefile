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

clean:
	/bin/rm -f -r *.o
