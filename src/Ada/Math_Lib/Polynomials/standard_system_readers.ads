with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;

package Standard_System_Readers is

-- DESCRIPTION :
--   Wrappers around routines to read a system from file,
--   as needed in various drivers.

  procedure Read_System
              ( file : in out file_type; filename : in string;
                p : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Read_System
              ( file : in out file_type; filename : in string;
                p : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys );

  -- DESCRIPTION :
  --   Tries to open the file with given name in the string filename
  --   with intent to read a polynomial system.  If this fails, then
  --   a message is written to the user and p will be null on return.
  --   Otherwise, the p on return will contain a polynomial system
  --   and the file on return is open for input.

end Standard_System_Readers;
