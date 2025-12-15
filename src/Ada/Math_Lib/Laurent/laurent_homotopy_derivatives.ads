with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;

package Laurent_Homotopy_Derivatives is

-- DESCRIPTION :
--   Provides functions to compute the derivatives of a Laurent homotopy
--   evaluated at some value for the continuation parameter t.

  function Eval ( dim : integer32;
                  cff : Standard_Complex_Vectors.Vector;
                  ctp : Standard_Floating_Vectors.Vector;
                  deg : Standard_Integer_VecVecs.VecVec; t : double_float )
                return Standard_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Evaluates one polynomial in the Laurent homotopy at a value 
  --   for the continuation parameter t.

  -- ON ENTRY :
  --   dim        number of variables in the polynomials;
  --   cff        coefficients of a polynomial in a Laurent homotopy;
  --   ctp        powers of t coefficients;
  --   deg        supports of one polynomial in the homotopy;
  --   t          value for the continuation parameter.

  function Diff ( p : Standard_Complex_Laurentials.Poly;
                  idx : Standard_Integer_Vectors.Vector;
                  vrblvl : integer32 := 0 )
                return Standard_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns the derivative of p according to the indices in idx.
  --   The verbose level is defined by vrblvl.

  function Eval ( cff : Standard_Complex_VecVecs.VecVec;
                  ctp : Standard_Floating_VecVecs.VecVec;
                  deg : Standard_Integer_VecVecs.Array_of_VecVecs;
                  t : double_float )
                return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Evaluates the Laurent homotopy at a value 
  --   for the continuation parameter t.

  -- ON ENTRY :
  --   cff        coefficients of the Laurent homotopy;
  --   ctp        powers of t coefficients;
  --   deg        supports of the polynomials in the homotopy;
  --   t          value for the continuation parameter.

  function Diff ( cff : Standard_Complex_VecVecs.VecVec;
                  ctp : Standard_Floating_VecVecs.VecVec;
                  deg : Standard_Integer_VecVecs.Array_of_VecVecs;
                  idx : Standard_Integer_Vectors.Vector;
                  z : Standard_Complex_Vectors.Vector; t : double_float;
                  vrblvl : integer32 := 0 )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the homotopy at t and computes the value of the derivative
  --   according to the indices at idx, evaluated at z.

  -- ON ENTRY :
  --   cff        coefficients of the Laurent homotopy;
  --   ctp        powers of t coefficients;
  --   deg        supports of the polynomials in the homotopy;
  --   idx        differential indices;
  --   z          values for the variables in the homotopy;
  --   t          value for the continuation parameter;
  --   vrblvl     is the verbose level.

end Laurent_Homotopy_Derivatives;
