with Standard_Complex_Ring;              use Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Generic_Matrices;

package Standard_Complex_Matrices is
  new Generic_Matrices(Standard_Complex_Ring,
                       Standard_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of standard complex numbers.
