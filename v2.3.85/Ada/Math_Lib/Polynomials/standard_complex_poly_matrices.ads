with Standard_Complex_Poly_Ring;         use Standard_Complex_Poly_Ring;
with Standard_Complex_Poly_Vectors;
with Generic_Matrices;

package Standard_Complex_Poly_Matrices is
  new Generic_Matrices(Standard_Complex_Poly_Ring,
                       Standard_Complex_Poly_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of standard complex polynomials.
