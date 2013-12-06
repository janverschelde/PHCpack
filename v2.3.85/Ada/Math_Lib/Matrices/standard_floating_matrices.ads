with Standard_Floating_Ring;              use Standard_Floating_Ring;
with Standard_Floating_Vectors;
with Generic_Matrices;

package Standard_Floating_Matrices is
  new Generic_Matrices(Standard_Floating_Ring,
                       Standard_Floating_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of standard floating-point numbers.
