with text_io;                          use text_io;
with Standard_Complex_Polynomials;     use Standard_Complex_Polynomials;

procedure Black_Box_Factorization
            ( infilename : in string; file : in file_type; p : in Poly );

-- DESCRIPTION :
--   Black-box factorization of one polynomial in several variables,
--   invoked as phc -b on an input which contains p as system.
