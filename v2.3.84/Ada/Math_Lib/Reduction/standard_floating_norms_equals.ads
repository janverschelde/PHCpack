with Standard_Floating_Ring;
with Standard_Floating_Ring.FField;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Generic_Norms_Equals;

package Standard_Floating_Norms_Equals is
  new Generic_Norms_Equals(Standard_Floating_Ring,
                           Standard_Floating_Ring.FField,
                           Standard_Floating_Vectors,
                           Standard_Floating_Matrices);

-- DESCRIPTION :
--   Defines norms and equalities for standard floating-point numbers.
