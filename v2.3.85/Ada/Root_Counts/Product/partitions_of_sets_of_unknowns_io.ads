with text_io;                            use text_io;
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

end Partitions_of_Sets_of_Unknowns_io;
