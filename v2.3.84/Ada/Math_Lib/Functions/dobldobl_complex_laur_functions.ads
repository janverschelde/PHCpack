with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Ring.FField;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Laurentials;
with Generic_Laur_Poly_Functions;

package DoblDobl_Complex_Laur_Functions is
  new Generic_Laur_Poly_Functions(DoblDobl_Complex_Ring,
                                  DoblDobl_Complex_Ring.FField,
                                  DoblDobl_Complex_Vectors,
                                  DoblDobl_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial functions for complex double doubles.
