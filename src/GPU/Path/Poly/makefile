gpp=/usr/bin/g++

CFLAG   = ../Complex
CDFLAG  = ../DefineTypesD

test_poly:
	$(gpp) -c -I$(CDFLAG) utilities.cpp -o utilities_d.o
	$(gpp) -c -I$(CDFLAG) poly.cpp -o poly_d.o
	$(gpp) -c -I$(CDFLAG) test_poly.cpp -o test_poly.o
	$(gpp) utilities_d.o poly_d.o test_poly.o -o /tmp/test_poly

clean:
	/bin/rm -f -r *.o
