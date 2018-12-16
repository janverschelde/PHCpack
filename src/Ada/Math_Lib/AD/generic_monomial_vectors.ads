with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors;
with Generic_Matrices;
with Generic_Monomials;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Matrices is new Generic_Matrices(Ring,Vectors);
  with package Monomials is new Generic_Monomials(Ring,Vectors,Matrices);

package Generic_Monomial_Vectors is

-- DESCRIPTION :
--   A monomial vector is a vector of monomials in several variables.

  type Monomial_Vector is 
    array ( integer32 range <> ) of Monomials.Link_to_Monomial;

  type Link_to_Monomial_Vector is access Monomial_Vector;

  procedure Clear ( v : in out Monomial_Vector );
  procedure Clear ( v : in out Link_to_Monomial_Vector );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the monomial vector m.

end Generic_Monomial_Vectors;
