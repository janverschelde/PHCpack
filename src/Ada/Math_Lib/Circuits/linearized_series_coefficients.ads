with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;

package Linearized_Series_Coefficients is

-- DESCRIPTION :
--   The procedures in this package convert linearized representations
--   into vectors of coefficient vectors of the power series.

  procedure Delinearize ( vy,yv : in Standard_Complex_VecVecs.VecVec );
  procedure Delinearize ( vy,yv : in DoblDobl_Complex_VecVecs.VecVec );
  procedure Delinearize ( vy,yv : in QuadDobl_Complex_VecVecs.VecVec );

  --  DESCRIPTION :
  --    Converts the linearized representation in vy into the vector yv
  --    of coefficient vectors of the series.
  --    This conversion is convenient for the difference computation.

  -- REQUIRED :
  --   if vy'range = 0..degree and vy(k)'range = 1..dimension,
  --   then yv'range = 1..dimension and yv(k)'range = 0..degree.

end Linearized_Series_Coefficients;
