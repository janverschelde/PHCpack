with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package Partitions_of_Sets_of_Unknowns_io is

-- DESCRIPTION :
--   This package provides i/o operations for partitions of
--   set of unknowns.

  procedure get ( p : in out Partition );
  procedure get ( file : in file_type; p : in out Partition );

  -- DESCRIPTION :
  --   A partition is read from standard input or from file.

  procedure put ( p : in Partition );
  procedure put ( file : in file_type; p : in Partition );
 
  -- DESCRIPTION :
  --   Writes a partition on standard output or on file.

-- INTERACTIVE INPUT :

  function iget ( m : natural32 ) return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a partition of the symbols into m sets,
  --   prompting the user for each symbol in the table for its set.
  --   On return is a vector of range 1..n, where n is the number of symbols.

  function Make_Partition
             ( n,m : natural32; p : Standard_Natural_Vectors.Vector )
             return Partition;

  -- DESCRIPTION :
  --   Returns the m-partition p, represented as a partition of
  --   a set of unknowns, a subset of the total set of n variables.
  --   The input vector p has range 1..n and p(i) defines the set to which
  --   the i-th variable belongs in the returned partition.

  function iget ( m : natural32 ) return Partition;

  -- DESCRIPTION :
  --   Returns a partition of the symbols into m sets,
  --   prompting the user for each symbol in the table for its set.

end Partitions_of_Sets_of_Unknowns_io;
