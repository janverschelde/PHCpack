with Standard_Complex_Ring;
with Standard_Complex_Laurentials;
with Generic_Laur_Poly_Systems;

package Standard_Complex_Laur_Systems is
  new Generic_Laur_Poly_Systems(Standard_Complex_Ring,
                                Standard_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial systems over the standard complex numbers.
