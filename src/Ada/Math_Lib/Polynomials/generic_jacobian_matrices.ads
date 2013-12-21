with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors,Generic_VecVecs;
with Generic_Matrices;
with Generic_Polynomials;
with Generic_Polynomial_Functions;
with Generic_Polynomial_Systems;
with Generic_Poly_System_Functions;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);
  with package Matrices is new Generic_Matrices(Ring,Vectors);
  with package Polynomials is new Generic_Polynomials(Ring);
  with package Poly_Functions is
         new Generic_Polynomial_Functions(Ring,Vectors,Polynomials);
  with package Poly_Systems is
         new Generic_Polynomial_Systems(Ring,Polynomials);
  with package Poly_SysFun is
         new Generic_Poly_System_Functions(Ring,Vectors,VecVecs,Polynomials,
                                           Poly_Functions,Poly_Systems);

package Generic_Jacobian_Matrices is 

-- DESCRIPTION :
--   This package provides data structures and evaluation functions for
--   Jacobian matrices of systems of polynomials.

  use Ring,Vectors,VecVecs,Matrices;
  use Polynomials,Poly_Functions,Poly_Systems,Poly_SysFun;

-- FUNCTION TYPE :

  type Evaluator is access function ( x : Vector ) return Matrix;

-- DATA STRUCTURES :

  type Jaco_Mat is array ( integer32 range <>, integer32 range <> ) of Poly;
  type Link_to_Jaco_Mat is access Jaco_Mat;
  type Eval_Jaco_Mat is 
    array ( integer32 range <>, integer32 range <> ) of Eval_Poly;
  type Link_to_Eval_Jaco_Mat is access Eval_Jaco_Mat;

  type Array_of_Jaco_Mat is array ( integer32 range <> ) of Link_to_Jaco_Mat;
  type Array_of_Eval_Jaco_Mat is
    array ( integer32 range <> ) of Link_to_Eval_Jaco_Mat;

  type Eval_Coeff_Jaco_Mat is
    array ( integer32 range <>, integer32 range <> ) of Eval_Coeff_Poly;
  type Link_to_Eval_Coeff_Jaco_Mat is access Eval_Coeff_Jaco_Mat;
  type Mult_Factors is
    array ( integer32 range <>, integer32 range <> ) of Link_to_Vector;
  type Link_to_Mult_Factors is access Mult_Factors;

  -- USAGE :
  --   p : Poly_Sys(1..n);
  --   j : Jaco_Mat(1..n,1..n) := Create(p);
  --  =>  j(i,j) = Diff(p(i),j)

-- CREATORS :

  function Create ( p : Poly_Sys ) return Jaco_Mat;

  -- REQUIRED :
  --   The number of the unknowns of each polynomial must be the same

  function Create ( j : Jaco_Mat ) return Eval_Jaco_Mat;

  procedure Create ( p : Poly_Sys; 
                     j : out Eval_Coeff_Jaco_Mat; m : out Mult_Factors );

-- EVALUATORS :

  function Eval ( j : Jaco_Mat;      x : Vector ) return Matrix; -- return j(x);
  function Eval ( j : Eval_Jaco_Mat; x : Vector ) return Matrix; -- return j(x);

  function Eval ( j : Eval_Coeff_Jaco_Mat; m : Mult_Factors; 
                  c : VecVec; x : Vector ) return Matrix;

    -- returns j(c,x) with c the coefficients of the original polynomials
  
-- DESTRUCTORS :

  procedure Clear ( j : in out Jaco_Mat );
  procedure Clear ( j : in out Link_to_Jaco_Mat );
  procedure Shallow_Clear ( j : in out Link_to_Jaco_Mat );
  procedure Clear ( j : in out Eval_Jaco_Mat );
  procedure Clear ( j : in out Link_to_Eval_Jaco_Mat );
  procedure Shallow_Clear ( j : in out Link_to_Eval_Jaco_Mat );
  procedure Clear ( j : in out Eval_Coeff_Jaco_Mat );
  procedure Clear ( j : in out Link_to_Eval_Coeff_Jaco_Mat );
  procedure Shallow_Clear ( j : in out Link_to_Eval_Coeff_Jaco_Mat );

  procedure Clear ( m : in out Mult_Factors );
  procedure Clear ( m : in out Link_to_Mult_Factors );
  procedure Shallow_Clear ( m : in out Link_to_Mult_Factors );

  -- DESCRIPTION :
  --   Deallocation of the occupied memory.  A shallow clear releases only
  --   the pointer and not the underlying content.

end Generic_Jacobian_Matrices;
