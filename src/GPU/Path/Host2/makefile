QD_ROOT=/usr/local/qd-2.3.17
QD_LIB=/usr/local/lib

gpp=/usr/bin/g++

CMPLX = ../Complex
POLYS = ../Poly2

test_mgs:
	$(gpp) -c -I$(CMPLX) test_mgs.cpp -o test_mgs.o
	$(gpp) test_mgs.o $(QD_LIB)/libqd.a -o /tmp/test_mgs

clean:
	/bin/rm -f -r *.o
