with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package Partitions_of_Sets_Strings is

-- DESCRIPTION :
--   Provides operations to write partitions of sets of unknowns
--   to strings and to parse strings into partitions.

  function to_String ( z : Partition ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the partition.

  function Number_of_Sets ( s : string ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of opening and closing braces in the string.

  function Parse ( s : string; n : natural32 ) return Partition;

  -- DESCRIPTION :
  --   Parses the string into a partition of a set of n unknowns.
  --   The sets are separated by curly braces { } and
  --   the symbols of the variables in a set are separated by spaces.

end Partitions_of_Sets_Strings;
