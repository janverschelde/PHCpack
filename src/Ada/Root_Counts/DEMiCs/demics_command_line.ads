with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Poly_Systems;

package DEMiCs_Command_Line is

-- DESCRIPTION :
--   This package offers functions and procedures to define a basic
--   file based command line interface to DEMiCs,
--   to compute all mixed cells by dynamic enumeration.
--   DEMiCs was developed by Tomohiko Mizutani, Akiko Takeda,
--   and Masakazu Kojima and is licensed under GNU GPL Version 2.

-- IMPORTANT NOTE :
--   For this interface to work, the source of DEMiCs 0.951 must
--   be modified to write the lifting values with precision 16.
--   In simplex.cpp of the folder SRC of the DEMiCs code, the statement
--   cout << fixed << setprecision(16);
--   is added before the writing of the lifting values.

  function Random_Name ( prefix : string ) return string;

  -- DESCRIPTION :
  --   Returns the prefix string (which may be empty)
  --   followed by a random 8-digit integer.

  procedure Prepare_Input
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                iptname : in string;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Prepares the input file for the polynomial system p,
  --   extracting the supports.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   iptname  name of the input file;
  --   verbose  flag to indicate if extra output is wanted.

  -- ON RETURN :
  --   supports are the supports of p, assuming fully mixed.

  procedure Call_DEMiCs
              ( infile,outfile : in string;
                execfile : in string := "/tmp/demics";
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Calls DEMiCs on the input file with name in infilename.
  --   The name of the output file is returned in rndname.

  -- ON ENTRY :
  --   infile   name of the input file with the input for DEMiCs;
  --   outfile  name of the output file to write the output of DEMiCs;
  --   execfile is the name of the absolute path to the executable;
  --   verbose  is a flag to write the command to screen.

  function Extract_Lifting_Values
             ( vals : string ) return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Given in vals a string with floats, separated by one space each,
  --   extracts the lifting values.

  procedure Parse_Lifting
              ( file : in file_type; dim : in integer32;
                lifting : out Standard_Floating_VecVecs.VecVec;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Extracts the lifting values from file.

  -- NOTE : the precision of the printing of the lifting values
  --   is increased to 16, changed from the original version. 

  -- REQUIRED : lifting'range = 1..dim.

  -- ON ENTRY :
  --   file     output file of demics, opened for input;
  --   dim      number of supports;
  --   verbose  flag to indicate that more output is wanted.

  -- ON RETURN :
  --   lifting  lifting values for each support.

  function Extract_Cell_Indices
              ( dim : integer32; vals : string;
                verbose : boolean := true )
              return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Given in vals are the indices to the points that span
  --   the mixed cells, extracts those indices and returns
  --   the vector of index vectors, for each component.
  --   The vector on return has range 1..2*dim,
  --   assuming the supports are fully mixed.

  procedure Parse_Cells
              ( file : in file_type; dim : in integer32;
                cells : out Lists_of_Integer_Vectors.List;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Extracts the mixed cells from file.

  -- ON ENTRY :
  --   file     output file of demics, opened for output;
  --   dim      dimension of the problem;
  --   verbose  flag to indicate for extra diagnostic output.

  -- ON RETURN :
  --   cells    indices to the points that span the mixed cells.

  procedure Process_Output
              ( dim : in integer32; filename : in string;
                mv : out natural32;
                lif : out Standard_Floating_VecVecs.VecVec;
                cells : out Lists_of_Integer_Vectors.List;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Writes to screen the mixed volume on the output file
  --   in the string filename.

  -- ON ENTRY :
  --   dim      dimension of the problem;
  --   filename is the name of the output file;
  --   verbose  flag to indicate is more output is wanted.

  -- ON RETURN :  
  --   mv       the mixed volume computed by DEMiCs;
  --   lif      random lifting generated for the supports;
  --   cells    indices to the points that span the mixed cells.

end DEMiCs_Command_Line;
