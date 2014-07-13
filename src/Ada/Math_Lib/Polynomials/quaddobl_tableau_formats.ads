with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;     use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;     use QuadDobl_Complex_Laur_Systems;

package QuadDobl_Tableau_Formats is

-- DESCRIPTION :
--   A tableau format of a polynomial system is a simple way to denote
--   multivariate polynomials and systems of equations, written in their
--   fully expanded form with quad double complex coefficients.
--
-- TABLEAU FORMAT DEFINITION :
--   Every term occupies one line, first the real and then the imaginary
--   part of the coefficient, followed by the exponents of the variables.
--   The first item in the tableau form of a polynomial is the number of
--   terms, then followed by the terms, each term on a separate line.
--   followed by as many lines as that number.
--   The first item in the tableau form of a polynomial system is the 
--   number of polynomials and variables, followed by the polynomials
--   in tableau format.
--
-- EXAMPLE :
--   the system in symbolic format represented as
--
--   2
--    x**2 + 4*y**2 - 4;
--           2*y**2 - x;
--
--   is in tableau format represented as
--
--   2
--   3
--    1.00000000000000E+00 0.00000000000000E+00 2 0
--    4.00000000000000E+00 0.00000000000000E+00 0 2
--   -4.00000000000000E+00 0.00000000000000E+00 0 0
--   2
--    2.00000000000000E+00 0.00000000000000E+00 0 2
--   -1.00000000000000E+00 0.00000000000000E+00 1 0

  procedure Write_Tableau 
              ( file : in file_type;
                t : in QuadDobl_Complex_Polynomials.Term );
  procedure Write_Tableau 
              ( file : in file_type;
                t : in QuadDobl_Complex_Laurentials.Term );

  -- DESCRIPTION :
  --   Writes the term t in tableau format to file.

  procedure Write_Tableau 
              ( file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly );
  procedure Write_Tableau 
              ( file : in file_type;
                p : in QuadDobl_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Writes the polynomial p in tableau format to file.

  procedure Write_Tableau ( file : in file_type; p : in Poly_Sys );
  procedure Write_Tableau ( file : in file_type; p : in Laur_Sys );

  -- DESCRIPTION :
  --   Writes the polynomial system p in tableau format to file.

  procedure Convert_Polynomial_into_Tableau_Format;
  procedure Convert_Laurent_into_Tableau_Format;

  -- DESCRIPTION :
  --   Prompts the user for an input system and the name of the output file
  --   to write the input system in a tableau format.

  procedure Read_Tableau
              ( file : in file_type; n : in natural32;
                t : in out QuadDobl_Complex_Polynomials.Term );
  procedure Read_Tableau
              ( file : in file_type; n : in natural32;
                t : in out QuadDobl_Complex_Laurentials.Term );

  -- DESCRIPTION :
  --   Reads one line from file, for a term in n variables.

  -- REQUIRED :
  --   The term t has allocated space for n exponents in t.dg.

  procedure Read_Tableau
              ( file : in file_type; n : in natural32;
                p : out QuadDobl_Complex_Polynomials.Poly );
  procedure Read_Tableau
              ( file : in file_type; n : in natural32;
                p : out QuadDobl_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Reads a polynomial in n variables in tableau format from file.

  procedure Read_Tableau
              ( file : in file_type; n : in natural32; p : out Poly_Sys );
  procedure Read_Tableau
              ( file : in file_type; n : in natural32; p : out Laur_Sys );

  -- DESCRIPTION :
  --   Reads an n-dimensional polynomial system p in tableau format from file.

  procedure Read_Tableau
              ( file : in file_type; p : out Link_to_Poly_Sys );
  procedure Read_Tableau
              ( file : in file_type; p : out Link_to_Laur_Sys );

  -- DESCRIPTION :
  --   Reads a polynomial system p in tableau format from file.

  procedure Extract_Coefficients_and_Exponents
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                c : out QuadDobl_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec );
  procedure Extract_Coefficients_and_Exponents_Copies
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                c : out QuadDobl_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec );
  procedure Extract_Coefficients_and_Exponents
              ( p : in QuadDobl_Complex_Laurentials.Poly;
                c : out QuadDobl_Complex_Vectors.Vector;
                e : out Standard_Integer_VecVecs.VecVec );
  procedure Extract_Coefficients_and_Exponents_Copies
              ( p : in QuadDobl_Complex_Laurentials.Poly;
                c : out QuadDobl_Complex_Vectors.Vector;
                e : out Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Extracts the coefficients and exponent vectors of the polynomial p.
  --   Without the suffix "_Copies", the exponent vectors in e on return
  --   are shared pointers to the exponents in p: no memory allocation.
  --   With the suffix "_Copies", deep copies of the exponents are made.

  -- REQUIRED : c'range = e'range = 1..Number_of_Terms(p).

  -- ON ENTRY :
  --   p        a polynomial with complex coefficients.

  -- ON RETURN :
  --   c        the coefficients of the terms in p;
  --   e        corresponding coefficients of the terms in p.

  procedure Extract_Coefficients_and_Exponents
              ( p : in Poly_Sys;
                c : out QuadDobl_Complex_VecVecs.VecVec;
                e : out Standard_Natural_VecVecs.Array_of_VecVecs );
  procedure Extract_Coefficients_and_Exponents_Copies
              ( p : in Poly_Sys;
                c : out QuadDobl_Complex_VecVecs.VecVec;
                e : out Standard_Natural_VecVecs.Array_of_VecVecs );
  procedure Extract_Coefficients_and_Exponents
              ( p : in Laur_Sys;
                c : out QuadDobl_Complex_VecVecs.VecVec;
                e : out Standard_Integer_VecVecs.Array_of_VecVecs );
  procedure Extract_Coefficients_and_Exponents_Copies
              ( p : in Laur_Sys;
                c : out QuadDobl_Complex_VecVecs.VecVec;
                e : out Standard_Integer_VecVecs.Array_of_VecVecs );

  -- DESCRIPTION :
  --   Extracts the coefficients and exponents of the system p.
  --   Without the suffix "_Copies", the exponent vectors in e on return
  --   are shared pointers to the exponents in p: no memory allocation.
  --   With the suffix "_Copies", deep copies of the exponents are made.

  -- REQUIRED : c'range = e'range = p'range.

  -- ON ENTRY :
  --   p        a (Laurent) polynomial system with complex coefficients.

  -- ON RETURN :
  --   c        c(i) has the coefficients of the terms in p(i);
  --   e        corresponding coefficients of the terms in p.

  procedure Convert_Tableau_into_Symbolic_Format;

  -- DESCRIPTION :
  --   Prompts the user for the names of input and output files.
  --   The input file should contain a system in tableau format,
  --   which will be written in symbolic form to the output file.

  procedure Main_Interactive_Driver;

  -- DESCRIPTION :
  --   This procedure implements an interactive driver
  --   to the format conversions.

end QuadDobl_Tableau_Formats;
