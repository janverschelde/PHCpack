with Standard_Dense_Series_io;
with Standard_Dense_Series_Ring;
with Abstract_Ring_io;

package Standard_Dense_Series_Ring_io is
  new Abstract_Ring_io(Standard_Dense_Series_Ring,
                       Standard_Dense_Series_io.get,
                       Standard_Dense_Series_io.put,
                       Standard_Dense_Series_io.put);

-- DESCRIPTION :
--   Defines the input/output for elements of the ring of truncated power
--   series with complex coefficients in standard double precision.
