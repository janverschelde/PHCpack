with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors;
with Generic_Matrices;
with Generic_Polynomials;
with Generic_Polynomial_Functions;
with Generic_Polynomial_Systems;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Matrices is new Generic_Matrices(Ring,Vectors);
  with package Polynomials is new Generic_Polynomials(Ring);
  with package Polynomial_Functions is
         new Generic_Polynomial_Functions(Ring,Vectors,Polynomials);
  with package Polynomial_Systems is
         new Generic_Polynomial_Systems(Ring,Polynomials);

package Generic_Hessian_Matrices is 

-- DESCRIPTION :
--   This package provides data structures and evaluation functions
--   for Hessian matrices of polynomials.

  use Ring,Vectors,Matrices;
  use Polynomials,Polynomial_Functions,Polynomial_Systems;

-- DATA STRUCTURES :

  type Hessian is array ( integer32 range <>, integer32 range <> ) of Poly;
  type Link_to_Hessian is access Hessian;

  type Array_of_Hessians is array ( integer32 range <> ) of Link_to_Hessian;
  type Link_to_Array_of_Hessians is access Array_of_Hessians;

-- CREATORS :

  function Create ( p : Poly ) return Hessian;
  function Create ( p : Poly ) return Link_to_Hessian;

  -- DESCRIPTION :
  --   Returns the Hessian matrix of p,
  --   which has as many rows and columns as the number of variables in p.

  function Create ( p : Poly_Sys ) return Array_of_Hessians;

  -- DESCRIPTION :
  --   The array on return has p'range 
  --   with in its i-th entry the Hessian for p(i).

  function Create ( p : Poly; k : integer32 ) return Hessian;
  function Create ( p : Poly; k : integer32 ) return Link_to_Hessian;

  -- DESCRIPTION :
  --   Returns the Hessian of p, skipping the k-th variable,
  --   which is typically the homotopy continuation parameter.
  --   The polynomials in the Hessian on return will have as many
  --   variables as p, including the k-th variable,
  --   but the number of rows and columns in the Hessian
  --   is one less than the number of variables of p.

  function Create ( p : Poly_Sys; k : integer32 ) return Array_of_Hessians;

  -- DESCRIPTION :
  --   The array on return has p'range 
  --   with in its i-th entry the Hessian for p(i),
  --   with k as its skipped variable.

-- EVALUATORS :

  function Eval ( h : Hessian; x : Vector ) return Matrix;
  function Eval ( h : Link_to_Hessian; x : Vector ) return Matrix;

  -- DESCRIPTION :
  --   Returns the evaluated Hessian matrix at x.
  --   If the Hessian is created with the index k, then the vector x 
  --   should have a value for the k-th parameter at index k.

  -- REQUIRED : x'range = 1..Number_of_Unknowns(h(i,j)) for all i and j.
  
-- DESTRUCTORS :

  procedure Clear ( h : in out Hessian );
  procedure Clear ( h : in out Link_to_Hessian );
  procedure Clear ( h : in out Array_of_Hessians );
  procedure Clear ( h : in out Link_to_Array_of_Hessians );

  -- DESCRIPTION :
  --   Deallocation of the occupied memory.

end Generic_Hessian_Matrices;
