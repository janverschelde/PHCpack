with text_io;                            use text_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_Maple_Solutions_io is

-- DESCRIPTION :
--   This package converts a list of solutions into a list format,
--   accessible to Maple. 

  -- FORMAT : We use a list of equations:
  --
  --       [[x = 1, y = 0],          or     [[x = 1 + 0*I,
  --        [x = 1, y = 0],                   y = 1 + 0*I ],
  --        [x = -1, y = 0],                 [x = -1 + 0*I,
  --        [x = -1, y = 0]];                 y = -1 + 0*I ]];
  --
  --   where each solution is put on a separate line for the space
  --   floating-point numbers usually require, the semi-colon concludes 
  --   an input given to Maple.  I = sqrt(-1).

  procedure put ( s : in Solution );
  procedure put ( file : in file_type; s : in Solution );
  procedure put ( sols : in Solution_List );
  procedure put ( file : in file_type; sols : in Solution_List );

  -- DESCRIPTION :
  --   Write the Maple format (see above) of the solution (list)
  --   to standard output or to file.

 -- function Write ( s : Solution ) return string;
 -- function Write ( sols : Solution_List ) return string;

  -- DESCRIPTION :
  --   Returns a string writing the solution(s) in the Maple format.

  procedure put_pairs ( sols : in Solution_List );
  procedure put_pairs ( file : in file_type; sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes the solution list to standard output or to file,
  --   in another format accessible to Maple.

  procedure get ( sols : out Solution_List );
  procedure get ( file : in file_type; sols : out Solution_List );

  -- DESCRIPTION :
  --   Reads solutions in Maple format from standard input or from file,
  --   in the format produced by put.

end Standard_Maple_Solutions_io;
