with Generic_Laurent_Polynomials;
with Standard_Complex_Ring;

package Standard_Complex_Laurentials is 
  new Generic_Laurent_Polynomials(Standard_Complex_Ring);

-- DESCRIPTION :
--   Defines the Laurent polynomials over the ring of standard complex 
--   numbers.  The "Laurential" is a contraction of "Laurent polynomial".
