with Standard_Complex_Ring;
with Standard_Complex_Polynomials;
with Generic_Lists_of_Terms;

package Standard_Complex_Term_Lists is
  new Generic_Lists_of_Terms(Standard_Complex_Ring,
                             Standard_Complex_Polynomials);

-- DESCRIPTION :
--   Defines lists of terms over the standard complex numbers.
