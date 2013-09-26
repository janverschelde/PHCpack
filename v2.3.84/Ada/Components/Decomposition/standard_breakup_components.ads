with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Standard_Breakup_Components is

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
                ( file : file_type; embsys : Poly_Sys; sols : Solution_List;
                  hyp : Standard_Complex_VecVecs.VecVec; level : natural32 )
                return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the interpolating polynomials whose zeros determine the
  --   hypersurface where the generic points are taken from.
  --   All sampling is done at once, which leads to massive oversampling,
  --   but with interesting comparisons.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results on;
  --   embsys     embedded polynomial system;
  --   sols       generic points are solutions with last component(s) = 0;
  --   hyp        random hyperplanes to project the solutions on;
  --   level      number of random hyperplanes added to the original system.

  function Incremental_Interpolate
                ( file : file_type; embsys : Poly_Sys; sols : Solution_List;
                  hyp : Standard_Complex_VecVecs.VecVec; level : natural32 )
                return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the interpolating polynomials whose zeros determine the
  --   hypersurface where the generic points are taken from.
  --   The construction of the interpolating polynomials is done
  --   incrementally, with a more economic sampling.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results on;
  --   embsys     embedded polynomial system;
  --   sols       generic points are solutions with last component(s) = 0;
  --   hyp        random hyperplanes to project the solutions on;
  --   level      number of random hyperplanes added to the original system.

  procedure Dynamic_Interpolate
                ( file : in file_type; embsys : in Poly_Sys;
                  level : in natural32; sols : in Solution_List;
                  hyp : in Standard_Complex_VecVecs.VecVec;
                  centproj : in boolean; subspaces : out Array_of_Poly_Sys;
                  pivots : out Standard_Integer_VecVecs.VecVec;
                  basepts : out Standard_Complex_VecVecs.Array_of_VecVecs;
                  basecard : out Standard_Natural_Vectors.Vector;
                  filters : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Performs projections from base points before the interpolation.
  --   This reduces the degrees of the filtering polynomials.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results on;
  --   embsys     embedded polynomial system;
  --   level      number of random hyperplanes added to original system;
  --   sols       generic points as solutions of embsys;
  --   hyp        random hyperplanes used in projections;
  --   centproj   true if central projections are allowed to be used.

  -- ON RETURN :
  --   subspaces  polynomials describing linear subspaces that contain
  --              the components;
  --   pivots     remaining variables after restriction to subspaces;
  --   basepts    based points used in skew line projections;
  --   basecard   number of base points used for each component;
  --   filters    interpolating polynomials.

end Standard_Breakup_Components;
