with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;          use Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;         use Multprec_Complex_Matrices;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;     use Multprec_Complex_Poly_Systems;
with Multprec_Complex_Jaco_Matrices;    use Multprec_Complex_Jaco_Matrices;
with Multprec_Complex_Solutions;        use Multprec_Complex_Solutions;

package Multprec_Deflate_Singularities is

-- DESCRIPTION :
--   This package implements deflation techniques for singularities.

  function Deflate ( f : Poly_Sys; m : natural32;
                     a : Matrix; c : Vector ) return Poly_Sys;

  -- DESCRIPTION :
  --   The deflation of a polynomial system with m multipliers consists
  --   of the original polynomial system augmented with m random 
  --   combinations of the columns of the jacobian matrix.

  -- REQUIRED :
  --   a'range(1) = 1..Number_of_Variables(f(f'first)),
  --   a'range(2) = 1..m, c'range = 1..m.

  -- ON ENTRY :
  --   f         system of n polynomial equations in N variables;
  --   m         number of multipliers to add;
  --   a         matrix to multiply Jacobian matrix with;
  --   c         coefficients of the last linear equation.

  -- ON RETURN :
  --   The system on return has N+m variables and 2*n+1 equations,
  --   where the last equation is linear in the m multipliers.

  function Deflate ( f : Poly_Sys; jm : Jaco_Mat; m : natural32;
                     a : Matrix; c : Vector ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the deflation of f, with in jm the Jacobian matrix of f.
  --   Same parameters as the above Deflate.

  function Deflate ( f : Poly_Sys; m,size : natural32 ) return Poly_Sys;
  function Deflate ( f : Poly_Sys; jm : Jaco_Mat;
                     m,size : natural32 ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the deflation of the system f with m multipliers,
  --   generating random constants of the given size.

  function Deflate_Corank_One
             ( f : Poly_Sys; size : natural32 ) return Poly_Sys;
  function Deflate_Corank_One
             ( f : Poly_Sys; jm : Jaco_Mat; size : natural32 ) return Poly_Sys;
  function Deflate_Corank_One ( f : Poly_Sys; c : Vector ) return Poly_Sys;
  function Deflate_Corank_One ( f : Poly_Sys; jm : Jaco_Mat; c : Vector )
                              return Poly_Sys;

  -- DESCRIPTION :
  --   If the corank of the Jacobian matrix jm of f is already one,
  --   then there is no need to right-multiply the Jacobian matrix jm
  --   with a random matrix, and the number of multipliers is n-1,
  --   where n equals the number of variables.

  -- REQUIRED : c'range = 1..n, n = Number_of_Variables(f(f'first)).

  procedure Multipliers ( f : in Poly_Sys; z : in Vector; m : in natural32;
                          lambda : out Vector; resid : out Floating_Number );

  -- DESCRIPTION :
  --   Computes the values for the multipliers using the last m equations
  --   in the system f and the current approximation of the root in z.

  -- ON ENTRY :
  --   f          polynomial system, deflated with m multipliers;
  --   z          current approximation of the root;
  --   m          number of multipliers added in deflation.

  -- ON RETURN :
  --   lambda     vector of range 1..m with values of multipliers;
  --   resid      residual of overconstrained system solved.

  function Strip_Multipliers ( t : Term; nv : natural32 ) return Term;
  function Strip_Multipliers ( p : Poly; nv : natural32 ) return Poly;

  -- DESCRIPTION :
  --   Removes all multipliers from the term t or polynomial p,
  --   leaving only the original first nv variables.

  -- REQUIRED : Number_of_Unknowns(p) >= nv.

  function Strip_Multipliers
              ( f : Poly_Sys; nq,nv : natural32 ) return Poly_Sys;

  -- DESCRIPTION :
  --   Removes all multipliers from the system f, leaving only the original
  --   first nq equations in the original first nv variables.

  -- REQUIRED : f'last >= nq and Number_of_Unknowns(f(i)) >= nv.

  function Strip_Multipliers ( s : Solution; nv : natural32 ) return Solution;
  function Strip_Multipliers
             ( s : Solution_List; nv : natural32 ) return Solution_List;

  -- DESCRIPTION :
  --   Returns the first original nv components of the solutions.

end Multprec_Deflate_Singularities;
