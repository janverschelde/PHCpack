with Standard_Floating_Ring;
with Standard_Floating_Ring.FField;
with Generic_Complex_Numbers;

package Standard_Complex_Numbers is 
  new Generic_Complex_Numbers(Standard_Floating_Ring,
                              Standard_Floating_Ring.FField);

-- DESCRIPTION :
--   Defines the complex numbers based on standard floating-point numbers.
