with Standard_Floating_Numbers_io;
with Standard_Floating_Ring;
with Abstract_Ring_io;

package Standard_Floating_Ring_io is
  new Abstract_Ring_io(Standard_Floating_Ring,
                       Standard_Floating_Numbers_io.get,
                       Standard_Floating_Numbers_io.put,
                       Standard_Floating_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for standard floating-point numbers.
