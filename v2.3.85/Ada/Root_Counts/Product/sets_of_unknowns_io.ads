with text_io;                            use text_io;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;

package Sets_of_Unknowns_io is

-- DESCRIPTION :
--   This package contains routines for input and output of set of unknows.

  procedure get ( s : in out Set );
  procedure get ( file : in file_type; s : in out Set );

  -- REQUIRED :
  --   The symbol table must be initialized!
  --   A set must begin with '{' and end with '}.

  procedure put ( s : in Set );
  procedure put ( file : in file_type; s : in Set );

end Sets_of_Unknowns_io;
