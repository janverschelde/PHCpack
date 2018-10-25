with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;

package Jacobian_Rabinowitsch_Trick is

-- DESSCRIPTION :
--   Formulating a polynomial system along the lines of the trick of 
--   Rabinowitsch allows to move singular solutions to infinity.

  function Identity_Matrix
              ( n : integer32 ) return Standard_Complex_Matrices.Matrix;
  function Identity_Matrix
              ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix;
  function Identity_Matrix
              ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the n-by-n identity matrix.

  procedure Add_Last_Multiplier
              ( p : in out Standard_Complex_Polynomials.Poly;
                d : in integer32 );
  procedure Add_Last_Multiplier
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                d : in integer32 );
  procedure Add_Last_Multiplier
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                d : in integer32 );

  -- DESCRIPTION :
  --   Updates the last equation as p*y - 1,
  --   where y is the variable at index d.

  function Jacobian_Rabinowitsch
              ( p : Standard_Complex_Poly_Systems.Poly_Sys )
              return Standard_Complex_Poly_Systems.Poly_Sys;
  function Jacobian_Rabinowitsch
              ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
              return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Jacobian_Rabinowitsch
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
              return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the system augmented with the Jacobian matrix which
  --   puts singular solutions at infinity.

end Jacobian_Rabinowitsch_Trick;
