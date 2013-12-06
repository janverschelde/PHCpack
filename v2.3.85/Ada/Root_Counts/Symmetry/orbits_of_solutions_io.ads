with text_io;                            use text_io;
with Orbits_of_Solutions;                use Orbits_of_Solutions;

package Orbits_of_Solutions_io is

-- DESCRIPTION :
--   This package contains output routines for the orbit structure of the
--   solutions set of a symmetric polynomial system.

  procedure put ( lorb : in List_of_Orbits );

  -- DESCRIPTION :
  --   Produces the output of the orbit structure of a solution list,
  --   on standard output.

  procedure put ( file : in file_type; lorb : in List_of_Orbits );

  -- DESCRIPTION :
  --   Produces the output of the orbit structure of a solution list,
  --   on a file.

end Orbits_of_Solutions_io;
