QD_ROOT=/usr/local/qd-2.3.17
QD_LIB=/usr/local/lib

gpp=/usr/bin/g++

CFLAG   = ../Complex

test_polymon: 
	$(gpp) -c -I$(CFLAG) test_polymon.cpp -o test_polymon.o
	$(gpp) test_polymon.o $(QD_LIB)/libqd.a -o /tmp/test_polymon

test_polyeq:
	$(gpp) -c -I$(CFLAG) test_polyeq.cpp -o test_polyeq.o
	$(gpp) test_polyeq.o $(QD_LIB)/libqd.a -o /tmp/test_polyeq

test_polysys:
	$(gpp) -c -I$(CFLAG) test_polysys.cpp -o test_polysys.o
	$(gpp) test_polysys.o $(QD_LIB)/libqd.a -o /tmp/test_polysys

test_poly:
	$(gpp) -c -I$(CFLAG) test_poly.cpp -o test_poly.o
	$(gpp) test_poly.o $(QD_LIB)/libqd.a -o /tmp/test_poly

test_polysol:
	$(gpp) -c -I$(CFLAG) test_polysol.cpp -o test_polysol.o
	$(gpp) test_polysol.o $(QD_LIB)/libqd.a -o /tmp/test_polysol

clean:
	/bin/rm -f -r *.o
