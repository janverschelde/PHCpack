with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;

package Homotopy_Membership_Target is

-- DESCRIPTION :
--   This package provides operations to make the target system
--   in a homotopy membership test.

  function Adjusted_Slices 
             ( sli : Standard_Complex_VecVecs.VecVec;
               sol : Standard_Complex_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec;
  function Adjusted_Slices 
             ( sli : DoblDobl_Complex_VecVecs.VecVec;
               sol : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_VecVecs.VecVec;
  function Adjusted_Slices 
             ( sli : QuadDobl_Complex_VecVecs.VecVec;
               sol : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Adjusts the constant coefficient of each slice such that
  --   the solution vector satisfies the equations,
  --   for double, double double, and quad double precision.

  function Adjusted_Target
             ( ep : Standard_Complex_Poly_Systems.Poly_Sys;
               dim : natural32;
               pnt : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Adjusted_Target
             ( ep : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               dim : natural32;
               pnt : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Adjusted_Target
             ( ep : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               dim : natural32;
               pnt : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Given an embedded polynomial system for a witness set
  --   of dimension dim and a point pnt,
  --   returns a copy of ep, where the last dim linear equations
  --   are adjusted so the point pnt satisfies those equations,
  --   in double, double double, and quad double precision.

  function Adjusted_Target
             ( ep : Standard_Complex_Laur_Systems.Laur_Sys;
               dim : natural32;
               pnt : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Laur_Systems.Laur_Sys;
  function Adjusted_Target
             ( ep : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               dim : natural32;
               pnt : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Adjusted_Target
             ( ep : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               dim : natural32;
               pnt : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Given an embedded Laurent Laurnomial system for a witness set
  --   of dimension dim and a point pnt,
  --   returns a copy of ep, where the last dim linear equations
  --   are adjusted so the point pnt satisfies those equations,
  --   in double, double double, and quad double precision.

end Homotopy_Membership_Target;
