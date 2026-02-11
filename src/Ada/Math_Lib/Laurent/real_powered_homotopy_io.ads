with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_Laurentials;

package Real_Powered_Homotopy_IO is

-- DESCRIPTION :
--   A real powered homotopy is a system of Laurent polynomials where
--   for every coefficient there is a corresponding real powered series.
--   In the output procedures below, the coefficients of the Laurent
--   polynomial are ignored, as they are replaced by the series.

-- WRITE OUTPUT :

  function to_string ( q : Standard_Complex_Laurentials.Poly;
                       c : Standard_Complex_VecVecs.VecVec;
                       p : Standard_Floating_VecVecs.VecVec;
                       t : character := 't'; vrblvl : integer32 := 0 )
                     return string;

  -- DESCRIPTION :
  --   Returns the string representation of a Laurent homotopy polynomial,
  --   defined by q and for the k-th coefficient in q there is the
  --   real powered series with coefficients in c(k) and powers in p(k),
  --   using the symbol t to write the series.
  --   The verbose level is given in vrblvl.

  -- ON ENTRY :
  --   q         a Laurent polynomial with complex coefficients;
  --   c         for the k-th coefficient in q, c(k) contains
  --             the coefficients of a real powered series;
  --   p         p(k) contains the powers of the series coefficient
  --             of the k-th term in q;
  --   t         symbol used to represent the series parameter;
  --   vrblvl    is the verbose level.

  procedure put ( q : in Standard_Complex_Laurentials.Poly;
                  c : in Standard_Complex_VecVecs.VecVec;
                  p : in Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 );
  procedure put ( file : in file_type;
                  q : in Standard_Complex_Laurentials.Poly;
                  c : in Standard_Complex_VecVecs.VecVec;
                  p : in Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Writes the Laurent homotopy polynomial with real powered series
  --   to standard output or to file, using the symbol t.

  procedure put_line ( q : in Standard_Complex_Laurentials.Poly;
                       c : in Standard_Complex_VecVecs.VecVec;
                       p : in Standard_Floating_VecVecs.VecVec;
                       t : in character := 't'; vrblvl : in integer32 := 0 );
  procedure put_line ( file : in file_type;
                       q : in Standard_Complex_Laurentials.Poly;
                       c : in Standard_Complex_VecVecs.VecVec;
                       p : in Standard_Floating_VecVecs.VecVec;
                       t : in character := 't'; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Writes the Laurent homotopy polynomial with real powered series
  --   to standard output or to file, using the symbol t,
  --   in new line format, with a new line for every term.

-- PARSE INPUT :

  function number_of_terms
             ( s : string; vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Assuming all coefficient series are enclosed within round brackets,
  --   that is: inside (), returns the number of terms in the string
  --   representation of a real powered Laurent homotopy polynomial.
  --   Returns -1 if the formatting appears to be wrong.
  --   The verbose level is given by the value of vrblvl.

  procedure parse_series
              ( s : in string;
                c : out Standard_Complex_VecVecs.VecVec;
                p : out Standard_Floating_VecVecs.VecVec;
                t : in character := 't'; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Assuming all coefficient series are enclosed within round brackets,
  --   that is: inside (), parses the string for the coefficients and
  --   the powers of the series, using t as the symbol.

  function extract_polynomial_string
             ( s : in string; m : in integer32;
               vrblvl : in integer32 := 0 ) return string;

  -- DESCRIPTION :
  --   Given in s is a string representation of a Laurent homotopy
  --   polynomial with m terms, returns the string representation
  --   of a Laurent polynomial with all series coefficients removed.

  procedure parse_string
              ( s : in string; n : in integer32;
                q : out Standard_Complex_Laurentials.Poly;
                c : out Standard_Complex_VecVecs.VecVec;
                p : out Standard_Floating_VecVecs.VecVec;
                t : in character := 't'; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Parses a string for a real powered series Laurent polynomial.

  -- ON ENTRY :
  --   s        string representation of a Laurent polynomial
  --            with real powered series as coefficients;
  --   n        number of variables in the string representations;
  --   t        symbol used in the series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   q        a Laurent polynomial;
  --   c        coefficients of the power series;
  --   p        powers of the series.

  procedure get ( n,size : in integer32;
                  q : out Standard_Complex_Laurentials.Poly;
                  c : out Standard_Complex_VecVecs.VecVec;
                  p : out Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 );
  procedure get ( file : in file_type; n,size : in integer32;
                  q : out Standard_Complex_Laurentials.Poly;
                  c : out Standard_Complex_VecVecs.VecVec;
                  p : out Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Reads a real powered Laurent polynomial from standard input
  --   or from file.

  -- ON ENTRY :
  --   file     must be opened for input;
  --   n        number of variables in the string representations;
  --   size     is (an upper bound on) the size of series coefficients;
  --   t        symbol used in the series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   q        a Laurent polynomial;
  --   c        coefficients of the power series;
  --   p        powers of the series.

end Real_Powered_Homotopy_IO;
