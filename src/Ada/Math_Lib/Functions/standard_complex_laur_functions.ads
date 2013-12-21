with Standard_Complex_Ring;
with Standard_Complex_Ring.FField;
with Standard_Complex_Vectors;
with Standard_Complex_Laurentials;
with Generic_Laur_Poly_Functions;

package Standard_Complex_Laur_Functions is
  new Generic_Laur_Poly_Functions(Standard_Complex_Ring,
                                  Standard_Complex_Ring.FField,
                                  Standard_Complex_Vectors,
                                  Standard_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial functions for standard complex numbers.
