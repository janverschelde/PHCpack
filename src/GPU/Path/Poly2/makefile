QD_ROOT=/usr/local/qd-2.3.17
QD_LIB=/usr/local/lib

gpp=/usr/bin/g++

CFLAG   = ../Complex

test_poly:
	$(gpp) -c -I$(CFLAG) test_poly.cpp -o test_poly.o
	$(gpp) test_poly.o $(QD_LIB)/libqd.a -o /tmp/test_poly

clean:
	/bin/rm -f -r *.o
