with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Matrices;
with DoblDobl_CSeries_Poly_Systems;

package Test_DoblDobl_CSeries_Systems is

-- DESCRIPTION :
--   Tests systems of series polynomials in double double precision.

  procedure Read_Series_Vector
              ( v : out DoblDobl_Complex_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 );

  -- DESCRIPTION :
  --   Prompts for a vector of series.  The dim is the number of variables
  --   in the system where the series will be evaluated.
  --   The idx is the index of the variable used as series variable.

  procedure Test_Evaluation
              ( p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 );

  -- DESCRIPTION :
  --   Prompts for a vector of series and evaluates the system p.
  --   The idx is the index of the variable used as series variable.

  procedure Write ( A : DoblDobl_Complex_Series_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the series in the matrix A to screen.

  procedure Test_Newton_Step
              ( p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 );

  -- DESCRIPTION :
  --   Prompts for a vector of series and applies one Newton step on p.
  --   The idx is the index of the variable used as series variable.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a system of series polynomials,
  --   displays the menu of test operations, asks for a choice,
  --   and then launches the corresponding test.

end Test_DoblDobl_CSeries_Systems;
