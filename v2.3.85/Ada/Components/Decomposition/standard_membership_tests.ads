with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers; 
with Standard_Natural_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Standard_Membership_Tests is

-- DESCRIPTION :
--   This package offers implementations of the membership test to decide
--   whether a point belongs to the component or not.

  function In_Subspace
                ( file : in file_type;
                  p : Poly_Sys; ind : integer32; 
                  x : Standard_Complex_Vectors.Vector;
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
                ( file : file_type; intpols : Poly_Sys; ind : integer32;
                  sols : Solution_List; hyp : Standard_Complex_VecVecs.VecVec;
                  level : integer32; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the solution in sols at the indicated position ind
  --   satisfies one of the polynomials in intpols(1..ind-1).

  -- ON ENTRY :
  --   file       to write diagnostics;
  --   intpols    interpolating polynomials;
  --   ind        current index of the polynomials and the solution;
  --   sols       checks only the solution at position given by ind;
  --   hyp        general hyperplanes used in the projection;
  --   level      number of added slices;
  --   tol        tolerance on the residual.

  function On_Component
                ( file : file_type; intpols : Poly_Sys; ind : integer32;
                  pivots : Standard_Integer_VecVecs.VecVec;
                  sols : Solution_List; hyp : Standard_Complex_VecVecs.VecVec;
                  level : integer32; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the solution in sols at the indicated position ind
  --   satisfies one of the polynomials in intpols(1..ind-1). 
  --   This is an extended version with extra argument pivots: these
  --   pivots indicate which variables remain after restriction to a
  --   linear subspace that contains the component.

  function On_Component
                ( file : file_type;
                  p : Poly_Sys; ind : integer32;
                  subspaces : Array_of_Poly_Sys;
                  pivots : Standard_Integer_VecVecs.VecVec;
                  basepts : Standard_Complex_VecVecs.Array_of_VecVecs;
                  basecard : Standard_Natural_Vectors.Vector;
                  sols : Solution_List;
                  hyp : Standard_Complex_VecVecs.VecVec;
                  level : integer32; tol : double_float ) return boolean;

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

end Standard_Membership_Tests;
