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

end Real_Powered_Homotopy_IO;
