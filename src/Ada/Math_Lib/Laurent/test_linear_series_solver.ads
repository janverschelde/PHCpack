with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_Real_Power_Series;
with Double_rpSeries_Vectors;
with Double_rpSeries_Matrices;

package Test_Linear_Series_Solver is

-- DESCRIPTION :
--   Tests solving a linear system of real power series.

  function Random_Series
             ( size : integer32 )
             return Double_Real_Power_Series.Link_to_Series;

  -- DESCRIPTION :
  --   Generates a random series with truncation index equal to size.

  function Random_Series_Vector
             ( dim,size : integer32 )
             return Double_rpSeries_Vectors.Vector;

  -- DESCRIPTION :
  --   Generates a vector of dimension dim of random series
  --   with truncation index equal to size.

  function Random_Series_Matrix
             ( dim,size : integer32 )
             return Double_rpSeries_Matrices.Matrix;

  -- DESCRIPTION :
  --   Generates a matrix of dimension dim of random series
  --   with truncation index equal to size.

  procedure Write ( x : in Double_rpSeries_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes all components of x.

  procedure Write ( A : in Double_rpSeries_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes all elements of A.

  function Right_Hand_Side
             ( A : Double_rpSeries_Matrices.Matrix;
               x : in Double_rpSeries_Vectors.Vector ) 
             return Double_rpSeries_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the product of A with x.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimensions, generates a random problem,
  --   and then tests the solver.

end Test_Linear_Series_Solver;
