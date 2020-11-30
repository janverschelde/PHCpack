with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;

package Main_Scaling is

-- DESCRIPTION :
--   Defines the main procedures to scale polynomial systems.

  procedure Standard_Read_System 
              ( file : in out file_type; filename : in string;
                dim : out integer32;
                lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Attempts to open the file defined by the name in filename
  --   and to read a polynomial system in standard double precision.

  procedure DoblDobl_Read_System 
              ( file : in out file_type; filename : in string;
                dim : out integer32;
                lp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Attempts to open the file defined by the name in filename
  --   and to read a polynomial system in double double precision.

  procedure QuadDobl_Read_System 
              ( file : in out file_type; filename : in string;
                dim : out integer32;
                lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Attempts to open the file defined by the name in filename
  --   and to read a polynomial system in quad double precision.

  procedure Standard_Separate_File 
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scalvec : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Prompts the user if the system, scaled with standard double
  --   precision arithmetic, needs to be written to a separate file.

  procedure DoblDobl_Separate_File 
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scalvec : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Prompts the user if the system, scaled with double double
  --   precision arithmetic, needs to be written to a separate file.

  procedure QuadDobl_Separate_File 
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scalvec : in QuadDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Prompts the user if the system, scaled with quad double
  --   precision arithmetic, needs to be written to a separate file.

  procedure Standard_Rescale
                ( dim : in integer32;
                  infile,outfile : in out file_type;
                  sysonfile : in boolean ); 

  -- DESCRIPTION :
  --   Rescales the solutions with respect to given scaling coefficients.
  --   The user is either prompted for solutions and scaling coefficients
  --   or these are read from the input file infile.
  --   The rescaling of the solutions happens in standard double precision.

  -- ON ENTRY :
  --   dim        number of polynomial equations in the system;
  --   infile     if sysonfile, then infile is the input file;
  --   outfile    for writing extra output information
  --   sysonfile  is true if the input system and scaling coefficients
  --              are available on the infile.

  -- ON RETURN :
  --   infile     closed for input;
  --   outfile    closed for outpute.
  
  procedure DoblDobl_Rescale
                ( dim : in integer32;
                  infile,outfile : in out file_type;
                  sysonfile : in boolean );

  -- DESCRIPTION :
  --   Rescales the solutions with respect to given scaling coefficients.
  --   The user is either prompted for solutions and scaling coefficients
  --   or these are read from the input file infile.
  --   The rescaling of the solutions happens in double double precision.

  -- ON ENTRY :
  --   dim        number of polynomial equations in the system;
  --   infile     if sysonfile, then infile is the input file;
  --   outfile    for writing extra output information
  --   sysonfile  is true if the input system and scaling coefficients
  --              are available on the infile.

  -- ON RETURN :
  --   infile     closed for input;
  --   outfile    closed for outpute.
  
  procedure QuadDobl_Rescale
                ( dim : in integer32;
                  infile,outfile : in out file_type;
                  sysonfile : in boolean );

  -- DESCRIPTION :
  --   Rescales the solutions with respect to given scaling coefficients.
  --   The user is either prompted for solutions and scaling coefficients
  --   or these are read from the input file infile.
  --   The rescaling of the solutions happens in quad double precision.

  -- ON ENTRY :
  --   dim        number of polynomial equations in the system;
  --   infile     if sysonfile, then infile is the input file;
  --   outfile    for writing extra output information
  --   sysonfile  is true if the input system and scaling coefficients
  --              are available on the infile.

  -- ON RETURN :
  --   infile     closed for input;
  --   outfile    closed for outpute.
  
  procedure Standard_Display_and_Dispatch_Menu
               ( infile,outfile : in out file_type; dim : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 sysonfile : in boolean );

  -- DESCRIPTION :
  --   Displays the menu and returns a choice, corresponding to one of the
  --   three available scaling procedures in standard double precision.

  procedure DoblDobl_Display_and_Dispatch_Menu
               ( infile,outfile : in out file_type; dim : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sysonfile : in boolean );

  -- DESCRIPTION :
  --   Displays the menu and returns a choice, corresponding to one of the
  --   three available scaling procedures in double double precision.

  procedure QuadDobl_Display_and_Dispatch_Menu
               ( infile,outfile : in out file_type; dim : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sysonfile : in boolean );

  -- DESCRIPTION :
  --   Displays the menu and returns a choice, corresponding to one of the
  --   three available scaling procedures in quad double precision.

  procedure Standard_Main ( infilename,outfilename : in string );

  -- DESCRIPTION :
  --   Parses the given system into standard double precision
  --   and then calls the procedure to display the scaling menu.

  procedure DoblDobl_Main ( infilename,outfilename : in string );

  -- DESCRIPTION :
  --   Parses the given system into double double precision
  --   and then calls the procedure to display the scaling menu.

  procedure QuadDobl_Main ( infilename,outfilename : in string );

  -- DESCRIPTION :
  --   Parses the given system into quad double precision
  --   and then calls the procedure to display the scaling menu.

  procedure Main ( infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines what phc -s does.
  --   The first two arguments are the names of the input and output files.
  --   The last argument is the verbose level.

end Main_Scaling;
