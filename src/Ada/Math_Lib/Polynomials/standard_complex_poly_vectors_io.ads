with Standard_Complex_Poly_Ring_io;
with Standard_Complex_Poly_Vectors;
with Generic_Vectors_io;

package Standard_Complex_Poly_Vectors_io is 
  new Generic_Vectors_io(Standard_Complex_Poly_Ring_io,
                         Standard_Complex_Poly_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of standard complex polynomials.
