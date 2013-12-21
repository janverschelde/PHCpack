with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers; 
with Standard_Natural_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;
with Multprec_Complex_Poly_Systems;

package Multprec_Membership_Tests is

-- DESCRIPTION :
--   This package offers implementations of the membership test to decide
--   whether a point belongs to the component or not.

  function In_Subspace
                ( file : in file_type;
                  p : Multprec_Complex_Poly_Systems.Poly_Sys;
                  ind : integer32; x : Standard_Complex_Vectors.Vector;
                  tol : double_float ) return boolean;
  function In_Subspace
                ( file : in file_type;
                  p : Multprec_Complex_Poly_Systems.Poly_Sys;
                  ind : integer32; x : Multprec_Complex_Vectors.Vector;
                  tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Given the polynomial representations of the linear container
  --   subspaces, the function "In_Subspace" returns true if the point x
  --   evaluates within the prescribed tolerance to lie inside the space.
  --   Otherwise false is returned.  Writes one diagnostic line on file.

  -- ON ENTRY :
  --   file       to write intermediate diagnostics on;
  --   p          linear equations for the span of the component;
  --   ind        index to the current component;
  --   x          point to evaluate in the equations;
  --   tol        tolerance on residual to decide true or false.

  function On_Component
                ( file : file_type;
                  p : Multprec_Complex_Poly_Systems.Poly_Sys; ind : integer32;
                  sols : Standard_Complex_Solutions.Solution_List;
                  hyp : Multprec_Complex_VecVecs.VecVec;
                  level : integer32; size : natural32; tol : double_float )
                return boolean;
  function On_Component
                ( file : file_type;
                  p : Multprec_Complex_Poly_Systems.Poly_Sys; ind : integer32;
                  sols : Multprec_Complex_Solutions.Solution_List;
                  hyp : Multprec_Complex_VecVecs.VecVec;
                  level : integer32; size : natural32; tol : double_float )
                return boolean;

  -- DESCRIPTION :
  --   Returns true if the solution in sols at the indicated position ind
  --   satisfies one of the polynomial in p(1..ind-1).

  -- ON ENTRY :
  --   file       to write diagnostics;
  --   p          interpolating polynomials;
  --   ind        current index of the polynomials and the solution;
  --   sols       checks only the solution at position given by ind;
  --   hyp        general hyperplanes used in the projection;
  --   level      number of added slices;
  --   tol        tolerance on the residual.

  function On_Component
                ( file : file_type;
                  p : Multprec_Complex_Poly_Systems.Poly_Sys; ind : integer32;
                  subspaces : Multprec_Complex_Poly_Systems.Array_of_Poly_Sys;
                  pivots : Standard_Integer_VecVecs.VecVec;
                  basepts : Multprec_Complex_VecVecs.Array_of_VecVecs;
                  basecard : Standard_Natural_Vectors.Vector;
                  sols : Standard_Complex_Solutions.Solution_List;
                  hyp : Multprec_Complex_VecVecs.VecVec;
                  level : integer32; size : natural32; tol : double_float )
                return boolean;
  function On_Component
                ( file : file_type;
                  p : Multprec_Complex_Poly_Systems.Poly_Sys; ind : integer32;
                  subspaces : Multprec_Complex_Poly_Systems.Array_of_Poly_Sys;
                  pivots : Standard_Integer_VecVecs.VecVec;
                  basepts : Multprec_Complex_VecVecs.Array_of_VecVecs;
                  basecard : Standard_Natural_Vectors.Vector;
                  sols : Multprec_Complex_Solutions.Solution_List;
                  hyp : Multprec_Complex_VecVecs.VecVec;
                  level : integer32; size : natural32; tol : double_float )
                return boolean;

  -- DESCRIPTION :
  --   Returns true if the solution in sols at the indicated position ind
  --   satisfies one of the polynomials in p(1..ind-1). 

  -- ON ENTRY :
  --   file       to write diagnostics on;
  --   p          filtering polynomial, one for each component;
  --   ind        current index of the generic point to be tested;
  --   subspaces  polynomials for the linear container space;
  --   pivots     indicate remaining variables after subspace restriction;
  --   basepts    base points used in the skew line projections;
  --   basecard   cardinalities of the set of base points.

end Multprec_Membership_Tests;
