with text_io;                            use text_io;
with Degree_Sets_Tables;                 use Degree_Sets_Tables;

package Degree_Sets_Tables_io is

-- DESCRIPTION :
--   This package provides output for the date structure to compute
--   a generalized permanent based on a set structure.

  procedure put ( ase : in Array_of_Sets );
  procedure put ( file : in file_type; ase : in Array_of_Sets );

  -- DESCRIPTION :
  --   Writes the array of sets ase to standard output or to file.

  procedure put ( dst : in Degree_Sets_Table );
  procedure put ( file : in file_type; dst : in Degree_Sets_Table );

  -- DESCRIPTION :
  --   Writes the degree sets table to standard output or to file.

end Degree_Sets_Tables_io;
