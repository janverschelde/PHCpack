with text_io;                         use text_io;
with Standard_Bracket_Systems;        use Standard_Bracket_Systems;

package Standard_Bracket_Systems_io is

-- DESCRIPTION :
--   This package provides output routines for systems of bracket polynomials.

  procedure put ( s : in Bracket_System );
  procedure put ( file : in file_type; s : in Bracket_System );

  -- DESCRIPTION :
  --   Writes the bracket system to standard output or to file.

  procedure put_line ( s : in Bracket_System );
  procedure put_line ( file : in file_type; s : in Bracket_System );

  -- DESCRIPTION :
  --   In writing the bracket system to standard output or to file,
  --   a new line is started for every monomial.

end Standard_Bracket_Systems_io;
