with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring,Abstract_Ring.Field;
with Generic_Vectors,Generic_VecVecs;
with Generic_Matrices;
with Generic_Laurent_Polynomials;
with Generic_Laur_Poly_Functions;
with Generic_Laur_Poly_Systems;
with Generic_Laur_System_Functions;

generic

  with package Ring is new Abstract_Ring(<>);
  with package FField is new Ring.Field(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);
  with package Matrices is new Generic_Matrices(Ring,Vectors);
  with package Polynomials is new Generic_Laurent_Polynomials(Ring);
  with package Poly_Functions is
         new Generic_Laur_Poly_Functions(Ring,FField,Vectors,Polynomials);
  with package Poly_Systems is
         new Generic_Laur_Poly_Systems(Ring,Polynomials);
  with package Poly_SysFun is
         new Generic_Laur_System_Functions
               (Ring,FField,Vectors,VecVecs,
                Polynomials,Poly_Functions,Poly_Systems);

package Generic_Laur_Jaco_Matrices is 

-- DESCRIPTION :
--   This package provides data structures and evaluation functions for
--   Jacobian matrices of systems of Laurent polynomials.

  use Ring,FField,Vectors,VecVecs,Matrices;
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
  --   p : Laur_Sys;
  --   j : Jaco_Mat := Create(p);
  --  =>  j(i,j) = Diff(p(i),j)

-- CREATORS :

  function Create ( p : Laur_Sys ) return Jaco_Mat;

  -- REQUIRED :
  --   The number of the unknowns of each polynomial must be the same

  function Create ( j : Jaco_Mat ) return Eval_Jaco_Mat;

  procedure Create ( p : Laur_Sys; 
                     j : out Eval_Coeff_Jaco_Mat; m : out Mult_Factors );

-- EVALUATORS :

  function Eval ( j : Jaco_Mat;      x : Vector ) return Matrix; -- return j(x);
  function Eval ( j : Eval_Jaco_Mat; x : Vector ) return Matrix; -- return j(x);

  function Eval ( j : Eval_Coeff_Poly; m,c,x : Vector ) return number;

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
  --   Deallocation of the occupied memory.

end Generic_Laur_Jaco_Matrices;
