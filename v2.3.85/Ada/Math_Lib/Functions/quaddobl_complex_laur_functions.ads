with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Ring.FField;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Laurentials;
with Generic_Laur_Poly_Functions;

package QuadDobl_Complex_Laur_Functions is
  new Generic_Laur_Poly_Functions(QuadDobl_Complex_Ring,
                                  QuadDobl_Complex_Ring.FField,
                                  QuadDobl_Complex_Vectors,
                                  QuadDobl_Complex_Laurentials);

-- DESCRIPTION :
--   Defines Laurent polynomial functions for complex quad doubles.
