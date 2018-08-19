with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Standard_Natural_Vectors;
with Generic_Polynomials;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Polynomials is new Generic_Polynomials(Ring);

package Generic_Polynomial_Systems is 

-- DESCRIPTION :
--   This package provides polynomial systems as an array of polynomials
--   with some arithmetic operations.

  use Ring,Polynomials;

-- DATA STRUCTURES :

  type Poly_Sys is array ( integer32 range <> ) of Poly;
  type Link_to_Poly_Sys is access Poly_Sys;
  type Array_of_Poly_Sys is array ( integer32 range <> ) of Link_to_Poly_Sys;
  type Link_to_Array_of_Poly_Sys is access Array_of_Poly_Sys;

-- COPYING :

  procedure Copy ( p : in Poly_Sys; q : in out Poly_Sys );

-- SELECTORS :

  function Variables_in_Support
             ( p : Poly_Sys ) return Standard_Natural_Vectors.Vector;
  -- returns 0/1 vector with 1 at place of a variable with nonzero exponent

  function Size_of_Support ( p : Poly_Sys ) return natural32;
  -- returns the number of variables with nonzero exponent

-- ARITHMETIC OPERATIONS :

  function "+" ( p,q : Poly_Sys ) return Poly_Sys;            -- return p+q
  function "-" ( p,q : Poly_Sys ) return Poly_Sys;            -- return p-q
  function "-" ( p : Poly_Sys ) return Poly_Sys;              -- return -p
  function "*" ( a : number; p : Poly_Sys ) return Poly_Sys;  -- return a*p
  function "*" ( p : Poly_Sys; a : number) return Poly_Sys;   -- return p*a

  procedure Add ( p : in out Poly_Sys; q : in Poly_Sys );     -- p := p+q
  procedure Sub ( p : in out Poly_Sys; q : in Poly_Sys );     -- p := p-q
  procedure Min ( p : in out Poly_Sys );                      -- p := -p
  procedure Mul ( p : in out Poly_Sys; a : in number );       -- p := a*p

-- DIFFERENTIATORS :

  function  Diff ( p : Poly_Sys; i : integer32 ) return Poly_Sys;
  procedure Diff ( p : in out Poly_Sys; i : in integer32 );

-- DESTRUCTORS :

  procedure Clear ( p : in out Poly_Sys );
  procedure Clear ( p : in out Link_to_Poly_Sys );
  procedure Shallow_Clear ( p : in out Link_to_Poly_Sys );
  procedure Clear ( p : in out Array_of_Poly_Sys );
  procedure Clear ( p : in out Link_to_Array_of_Poly_Sys );

  -- DESCRIPTION :
  --   Deallocation of memory, a shallow clear only releases the pointer
  --   to the system; the system is not destroyed.

end Generic_Polynomial_Systems;
