with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;

package Affine_Transformations is

-- DESCRIPTION :
--   Given a polynomial system in 1-homogeneous coordinates,
--   replaces in every polynomial the homogeneous variable by one
--   and removes the last linear equation in the system.

  function Make_Affine
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Make_Affine
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Make_Affine
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the affine version of the homogeneous system p,
  --   in double, double double, and quad double precision.

end Affine_Transformations;
