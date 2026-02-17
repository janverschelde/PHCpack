with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;

package Real_Powered_Homotopy is

-- DESCRIPTION :
--   A real powered homotopy is represented by a 3-tuple
--   1) a system of Laurent polynomials,
--   2) for every monomial we have a vector of coefficients, and
--   3) for every monomial we have a vector of real powers.

  function Is_Linear ( p : Standard_Complex_Laurentials.Poly )
                     return boolean;
  function Is_Linear ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                     return boolean;

  -- DESCRIPTION :
  --   Returns true if the monomials in p define a linear system.
  --   Returns false otherwise.

  procedure Get_Constant_Coefficients
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                c : in Standard_Complex_VecVecs.Array_of_VecVecs;
                A : out Standard_Complex_Matrices.Matrix;
                b : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given the Laurent monomials in q and the corresponding
  --   coefficients in c, returns the linear system (A, b) defined
  --   by the constant coefficients of the real powered homotopy.
  --   The verbose level is given in vrblvl.

end Real_Powered_Homotopy;
