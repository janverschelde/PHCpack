with Standard_Dense_Series2_io;
with Standard_Dense_Series2_Ring;
with Abstract_Ring_io;

package Standard_Dense_Series2_Ring_io is
  new Abstract_Ring_io(Standard_Dense_Series2_Ring,
                       Standard_Dense_Series2_io.get,
                       Standard_Dense_Series2_io.put,
                       Standard_Dense_Series2_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in standard double precision.
