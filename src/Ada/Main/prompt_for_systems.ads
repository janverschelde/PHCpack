with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;

package Prompt_for_Systems is

-- DESCRIPTION :
--   By default, the input system is on the file, with the given file name,
--   but otherwise the user will be prompted to provide a file name,
--   or given the possibility to type in the polynomials in the system.

  procedure Scan_System
              ( file : in out file_type; name : in string;
                lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                onfile : out boolean );
  procedure Scan_System
              ( file : in out file_type; name : in string;
                lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                onfile : out boolean );
  procedure Scan_System
              ( file : in out file_type; name : in string;
                lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                onfile : out boolean );
  procedure Scan_System
              ( file : in out file_type; name : in string;
                lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                onfile : out boolean );
  procedure Scan_System
              ( file : in out file_type; name : in string;
                lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                onfile : out boolean );
  procedure Scan_System
              ( file : in out file_type; name : in string;
                lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                onfile : out boolean );

  -- DESCRIPTION :
  --   Checks whether the given file name corresponds to a file with
  --   a polynomial system in a correct format.
  --   If this is the case, then onfile is true on return and lp
  --   contains the system.
  --   The scan_system procedures are called in the read_system procedures.

  procedure Read_System
              ( file : in out file_type; name : in string;
                lp : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                onfile : out boolean );
  procedure Read_System
              ( file : in out file_type; name : in string;
                lp : in out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                onfile : out boolean );
  procedure Read_System
              ( file : in out file_type; name : in string;
                lp : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                onfile : out boolean );
  procedure Read_System
              ( file : in out file_type; name : in string;
                lp : in out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                onfile : out boolean );
  procedure Read_System
              ( file : in out file_type; name : in string;
                lp : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                onfile : out boolean );
  procedure Read_System
              ( file : in out file_type; name : in string;
                lp : in out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                onfile : out boolean );

  -- DESCRIPTION :
  --   Reads a system from file or from standard input.
  --   Coefficients are parsed in standard double, double double,
  --   or quad double precision.

  -- ON ENTRY :
  --   name     name of the file where a system may be.

  -- ON RETURN :
  --   file     file opened for input, if onfile,
  --            this file may contain the solutions of the system;
  --   lp       points to a system read from file or type in;
  --   onfile   is true if the system lp was read from file.

end Prompt_for_Systems;
