BIN=../bin

# first part of the library: the basics

ts_dcmplx: dcmplx.c ts_dcmplx.c
	gcc -o $(BIN)/ts_dcmplx dcmplx.c ts_dcmplx.c -lm

ts_roots: dcmplx.c dc_roots.c ts_roots.c
	gcc -o $(BIN)/ts_roots dcmplx.c dc_roots.c ts_roots.c -lm

ts_dc_determinant: dcmplx.c dc_matrix.c ts_dc_determinant.c
	gcc -o $(BIN)/ts_dc_determinant dcmplx.c dc_matrix.c \
	       ts_dc_determinant.c -lm

ts_timer: timer.c
	gcc -o $(BIN)/ts_timer ts_timer.c timer.c

lib.a: 
	gcc -c dcmplx.c  dc_matrix.c  dc_inverse.c 
	ar ruv lib.a dcmplx.o  dc_matrix.o  dc_inverse.o

ts_dc_inverse: ts_dc_inverse.c lib.a
	gcc -o $(BIN)/ts_dc_inverse ts_dc_inverse.c lib.a -lm

# second part: numerical GCD and Smith form

ts_gcd: dcmplx.c dc_roots.c dc_interpolation.c poly_gcd.c \
        poly_dcmplx.c ts_gcd.c
	gcc -o $(BIN)/ts_gcd dcmplx.c dc_roots.c dc_interpolation.c \
               poly_gcd.c poly_dcmplx.c ts_gcd.c -lm

ts_interpolation: dcmplx.c dc_roots.c dc_interpolation.c \
                  poly_dcmplx.c poly_gcd.c ts_interpolation.c
	 gcc -o $(BIN)/ts_interpolation dcmplx.c dc_roots.c \
                poly_dcmplx.c poly_gcd.c \
	        dc_interpolation.c ts_interpolation.c -lm 

ts_hermite: dcmplx.c dc_roots.c dc_interpolation.c \
            poly_dcmplx.c poly_gcd.c poly_hermite.c poly_matrix.c \
            poly_smith.c ts_hermite.c dc_matrix.c
	gcc -o $(BIN)/ts_hermite dcmplx.c dc_roots.c dc_interpolation.c \
               poly_dcmplx.c poly_gcd.c poly_hermite.c poly_matrix.c \
               poly_smith.c ts_hermite.c dc_matrix.c -lm

ts_smith: dcmplx.c dc_roots.c dc_interpolation.c poly_dcmplx.c \
          poly_gcd.c poly_smith.c poly_matrix.c dc_matrix.c ts_smith.c
	gcc -o $(BIN)/ts_smith dcmplx.c dc_roots.c dc_interpolation.c \
               poly_dcmplx.c poly_gcd.c \
               poly_smith.c poly_matrix.c dc_matrix.c ts_smith.c -lm

ts_multiply: dcmplx.c dc_roots.c dc_interpolation.c \
             poly_gcd.c poly_dcmplx.c poly_smith.c poly_matrix.c \
             ts_multiply.c
	gcc -o $(BIN)/ts_multiply dcmplx.c dc_roots.c \
               dc_interpolation.c poly_gcd.c poly_dcmplx.c \
               poly_matrix.c poly_smith.c ts_multiply.c -lm

ts_divide: dcmplx.c dc_roots.c dc_interpolation.c poly_gcd.c poly_dcmplx.c \
           poly_matrix.c poly_smith.c ts_divide.c
	gcc -o $(BIN)/ts_divide dcmplx.c dc_roots.c dc_interpolation.c \
               poly_gcd.c poly_dcmplx.c poly_matrix.c poly_smith.c \
               ts_divide.c -lm

ts_poly_inverse: dcmplx.c poly_dcmplx.c poly_matrix.c poly_smith.c \
                 dc_roots.c dc_interpolation.c poly_gcd.c ts_poly_inverse.c
	gcc -o $(BIN)/ts_poly_inverse dcmplx.c poly_dcmplx.c \
               poly_matrix.c poly_smith.c dc_roots.c dc_interpolation.c \
               poly_gcd.c ts_poly_inverse.c -lm

# third part: realization algorithms

ts_c2ada: dcmplx.c dc_roots.c dc_interpolation.c poly_gcd.c poly_dcmplx.c \
          poly_matrix.c poly_smith.c c2ada_poly_matrix.c ts_c2ada.c 
	gcc -o $(BIN)/ts_c2ada dcmplx.c dc_roots.c dc_interpolation.c \
               poly_dcmplx.c poly_gcd.c poly_matrix.c poly_smith.c \
               c2ada_poly_matrix.c ts_c2ada.c -lm

read_input: dcmplx.c dc_matrix.c dc_inverse.c c2ada_dc_matrix.c read_input.c
	gcc -o $(BIN)/read_input dcmplx.c dc_matrix.c dc_inverse.c \
               c2ada_dc_matrix.c read_input.c -lm

ts_realization: dcmplx.c dc_matrix.c dc_roots.c dc_interpolation.c \
                dc_inverse.c poly_dcmplx.c poly_gcd.c poly_matrix.c \
                poly_hermite.c poly_smith.c append_polymatrix.c \
                realization.c ts_realization.c
	gcc -o $(BIN)/ts_realization dcmplx.c dc_matrix.c dc_roots.c \
               dc_interpolation.c dc_inverse.c poly_dcmplx.c poly_gcd.c \
               poly_matrix.c poly_hermite.c poly_smith.c \
	       append_polymatrix.c realization.c ts_realization.c -lm

clean:
	/bin/rm -f *.o lib.a *.exe *~
