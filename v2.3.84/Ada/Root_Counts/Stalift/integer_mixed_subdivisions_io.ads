with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;  

package Integer_Mixed_Subdivisions_io is

-- DESCRIPTION :
--   This package provides some routines for i/o of mixed subdivisions.

  procedure get ( n,m : in natural32; mic : out Mixed_Cell );
  procedure get ( file : in file_type;
                  n,m : in natural32; mic : out Mixed_Cell );

  -- DESCRIPTION :
  --   Reads the normal and for each list of points, the length
  --   of the list and the list itself from standard input or from file.

  procedure get ( n,m : out natural32; mixed_type : out Link_to_Vector;
                  mixsub : out Mixed_Subdivision );
  procedure get ( file : in file_type; n,m : out natural32;
                  mixed_type : out Link_to_Vector;
                  mixsub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Reads the dimension, the number of different supports,
  --   the type of mixture, the number of mixed cells and the mixed
  --   subdivision from standard input or from file.

  procedure put ( n : in natural32; mix : in Vector; mic : in Mixed_Cell );
  procedure put ( n : in natural32; mix : in Vector; mic : in out Mixed_Cell;
                  mv : out natural32 );
  procedure put ( file : in file_type;
                  n : in natural32; mix : in Vector; mic : in Mixed_Cell );
  procedure put ( file : in file_type;
                  n : in natural32; mix : in Vector; mic : in out Mixed_Cell;
                  mv : out natural32 );

  -- DESCRIPTION :
  --   Puts the normal, the length of each point list and the
  --   points belonging to the cell on standard output or on file.
  --   More text banners are provided when the parameter `mv' is supplied.
  --   This `mv' contains the mixed volume of the cell on return.
  --   When the mixed volume is computed, eventually the cell is refined.

  procedure put ( n : in natural32; mix : in Vector;
                  mixsub : in Mixed_Subdivision );
  procedure put ( n : in natural32; mix : in Vector;
                  mixsub : in out Mixed_Subdivision; mv : out natural32 );
  procedure put ( file : in file_type; n : in natural32; mix : in Vector;
                  mixsub : in Mixed_Subdivision );
  procedure put ( file : in file_type; n : in natural32; mix : in Vector;
                  mixsub : in out Mixed_Subdivision; mv : out natural32 );

  -- DESCRIPTION :
  --   Puts the dimension, the type of mixture, the number of mixed
  --   cells and the mixed subdivision on file or on standard output.
  --   More text banners are provided when the parameter `mv' is supplied.
  --   This `mv' contains the mixed volume of the cell on return.
  --   When the mixed volume is computed, eventually the cells are refined.

end Integer_Mixed_Subdivisions_io;
