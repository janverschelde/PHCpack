with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;

package Multiplicity_Sampling_Machine is

-- DESCRIPTION :
--   This package implements a sampling machine for positive dimensional
--   solutions components with multiplicity higher than one.

  procedure Initialize ( sp : in Standard_Complex_Poly_Systems.Poly_Sys;
                         mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                         dim,ind : in natural );

  -- DESCRIPTION :
  --   Initializes the sampling machine with the sliced and embedded
  --   system with coefficients in standard and multi-precision.
  --   The format of the embedding must allow two types of moving
  --   equations: the usual slices to get the generic points, and an
  --   equation of type zz = z(eps) for the tube around the component.

  -- ON ENTRY :
  --   sp       sliced and embedded system with standard coefficients;
  --   mp       sliced and embedded system with multi-precision coefficients;
  --   dim      dimension of the component;
  --   ind      index of the moving equation that regulates the distance
  --            to the component, typically this is zz = z(eps).

  procedure Clear;

  -- DESCRIPTION :
  --   Destroys the internal data structures.

end Multiplicity_Sampling_Machine;
