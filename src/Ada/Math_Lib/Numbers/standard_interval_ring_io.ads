with Standard_Interval_Numbers_io;
with Standard_Interval_Ring;
with Abstract_Ring_io;

package Standard_Interval_Ring_io is
  new Abstract_Ring_io(Standard_Interval_Ring,
                       Standard_Interval_Numbers_io.get,
                       Standard_Interval_Numbers_io.put,
                       Standard_Interval_Numbers_io.put);

-- DESCRIPTION :
--   Defines the input/output for standard interval numbers.
