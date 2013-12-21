with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_VecVecs;
with Multprec_Complex_VecVecs;
with Standard_Complex_Solutions;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;

package Multprec_Breakup_Components is

-- DESCRIPTION :
--   This package offers three strategies to implement the break up of
--   equidimensional components into irreducibles:
--    1) Massive interpolate constructs a full grid and returns for each
--       component as many interpolating polynomials as its degree.
--       While this is expensive, it is interesting to check the accuracy.
--    2) Incremental interpolate only returns one interpolating polynomial
--       for each component.
--    3) Dynamic interpolate also works incrementally but constructs
--       the linear subspaces spanned by the components and uses central
--       projections to random hyperplanes.
--   This third strategy is the most efficient one.

  function Massive_Interpolate
                ( file : file_type;
                  embsys : Standard_Complex_Poly_Systems.Poly_Sys;
                  orgsys : Multprec_Complex_Poly_Systems.Poly_Sys;
                  sols : Standard_Complex_Solutions.Solution_List;
                  hyp : Multprec_Complex_VecVecs.VecVec;
                  level,size : natural32 )
                return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the interpolating polynomials whose zeros determine the
  --   hypersurface where the generic points are taken from.
  --   This massive version starts with doing the sampling all at once
  --   and does a huge oversampling in this way, but is interesting.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results on;
  --   embsys     embedded system;
  --   sols       generic points are solutions with last component(s) = 0;
  --   hyp        random hyperplanes to project the solutions on;
  --   level      number of random hyperplanes added to the original system;
  --   size       size of the coefficients in the interpolator.

  function Incremental_Interpolate
                ( file : file_type;
                  embsys : Standard_Complex_Poly_Systems.Poly_Sys;
                  orgsys : Multprec_Complex_Poly_Systems.Poly_Sys;
                  sols : Standard_Complex_Solutions.Solution_List;
                  hyp : Multprec_Complex_VecVecs.VecVec;
                  level,size : natural32 )
                return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the interpolating polynomials whose zeros determine the
  --   hypersurface where the generic points are taken from.
  --   The algorithm proceeds incrementally in sampling the points.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results on;
  --   embsys     embedded system;
  --   sols       generic points are solutions with last component(s) = 0;
  --   hyp        random hyperplanes to project the solutions on;
  --   level      number of random hyperplanes added to the original system;
  --   size       size of the coefficients in the interpolator.

  procedure Dynamic_Interpolate
                ( file : in file_type;
                  embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                  orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_List;
                  level,size : in natural32;
                  hyp : in Multprec_Complex_VecVecs.VecVec;
                  centproj : in boolean;
                  subspaces
                    : out Multprec_Complex_Poly_Systems.Array_of_Poly_Sys;
                  pivots : out Standard_Integer_VecVecs.VecVec;
                  basepts : out Multprec_Complex_VecVecs.Array_of_VecVecs;
                  basecard : out Standard_Natural_Vectors.Vector;
                  filters
                    : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Performs projections from base points before the interpolation.
  --   This reduces the degrees of the filtering polynomials.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results on;
  --   embsys     embedded polynomial system;
  --   orgsys     original system with multi-precision coefficients;
  --   sols       generic points as solutions of embsys;
  --   level      number of random hyperplanes added to original system;
  --   size       size of the multi-precision numbers;
  --   hyp        general hyperplanes used in the projections;
  --   centproj   if true, then central projections will be used,
  --              otherwise only subspace restrictions are done.

  -- ON RETURN :
  --   subspaces  polynomials describing linear subspaces that contain
  --              the components;
  --   pivots     remaining variables after restriction to the subspaces;
  --   basepts    base points used in skew line projections;
  --   basecard   number of base points used for each component.
  --   filters    interpolating polynomials.

end Multprec_Breakup_Components;
