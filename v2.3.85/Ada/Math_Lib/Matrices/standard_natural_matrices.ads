with Standard_Natural_Ring;              use Standard_Natural_Ring;
with Standard_Natural_Vectors;
with Generic_Matrices;

package Standard_Natural_Matrices is
  new Generic_Matrices(Standard_Natural_Ring,
                       Standard_Natural_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of standard natural numbers.
