with text_io;                          use text_io;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;

package Black_Box_Factorization is

-- DESCRIPTION :
--   This package offers the factorization routines in blackbox mode.

  procedure Standard_Black_Box_Factorization
              ( infilename : in string; file : in file_type;
                p : in Standard_Complex_Polynomials.Poly );
  procedure DoblDobl_Black_Box_Factorization
              ( infilename : in string; file : in file_type;
                p : in DoblDobl_Complex_Polynomials.Poly );
  procedure QuadDobl_Black_Box_Factorization
              ( infilename : in string; file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Black-box factorization of one polynomial in several variables,
  --   in standard double, double double, or quad double precision,
  --   invoked as phc -b on an input which contains p as system.

end Black_Box_Factorization;
