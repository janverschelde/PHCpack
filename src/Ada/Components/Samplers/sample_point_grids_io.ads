with text_io;                             use text_io;
with Sample_Point_Grids;                  use Sample_Point_Grids;

package Sample_Point_Grids_io is

-- DESCRIPTION :
--   This package provides input-output routines for grids
--   of sampled points.

  procedure get ( grid,grid_last : in out Standard_Sample_Grid );
  procedure get ( file : in file_type;
                  grid,grid_last : in out Standard_Sample_Grid );
  procedure get ( grid,grid_last : in out Multprec_Sample_Grid );
  procedure get ( file : in file_type;
                  grid,grid_last : in out Multprec_Sample_Grid );

  -- DESCRIPTION :
  --   Reads a grid from standard input or from file.  The parameter grid
  --   will point to the first element in the list, while grid_last will
  --   be a pointer to the last element in the list.
  --   The input format must be the same as the output format, see below.

  procedure put ( grid : in Standard_Sample_Grid );
  procedure put ( file : in file_type; grid : in Standard_Sample_Grid );
  procedure put ( grid : in Multprec_Sample_Grid );
  procedure put ( file : in file_type; grid : in Multprec_Sample_Grid );

  -- DESCRIPTION :
  --   Writes a grid on standard output or on file.
  --   The output format consists of the number of lists, followed
  --   by the lists, each in their standard format, first the three
  --   numbers with their dimensions (length,#variables,#slices),
  --   followed by the actual samples.

end Sample_Point_Grids_io;
