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
                mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
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
  --   mix      type of mixture of the supports;
  --   perm     permutation of the original position of supports;
  --   supports are the supports of p, permuted so the same supports
  --            appear jointly one after the other.

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
              ( file : in file_type;
                lifting : out Standard_Floating_VecVecs.VecVec;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Extracts the lifting values from file.

  -- NOTE : the precision of the printing of the lifting values by demics
  --   is increased to 16, changed from the original version. 

  -- REQUIRED : lifting'range corresponds to the type of mixture.

  -- ON ENTRY :
  --   file     output file of demics, opened for input;
  --   verbose  flag to indicate that more output is wanted.

  -- ON RETURN :
  --   lifting  lifting values for each support.

  function Offset_for_Index
              ( mix : Standard_Integer_Vectors.Vector;
                idx : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the offset for the index of the component idx,
  --   relative to the type of mixture in mix.

  function Extract_Cell_Indices
              ( nbrpts : integer32;
                mix : Standard_Integer_Vectors.Vector;
                vals : string; verbose : boolean := true )
              return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Given in vals are the indices to the points that span
  --   a fine mixed cell, extracts those indices and returns
  --   the vector of indices, parsed as integer numbers.
  --   The vector on return has range 1..nbrpts,
  --   where nbrpts is the number of points in each fine mixed cell.
  --   The type of mixture is given by mix.

  procedure Line2Cell_Indices
              ( line : in string; nbrpts : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                indices : out Standard_Integer_Vectors.Vector;
                verbose : boolean := true );

  -- DESCRIPTION :
  --   Given a line in the output of DEMiCs, converts the characters
  --   into a vector with integer indices.

  -- ON ENTRY :
  --   line     the output line starts with the cell number;
  --   nbrpts   the number of points in a fine mixed cell;
  --   mix      type of mixture;
  --   verbose  flag to indicate if extra output is wanted.

  -- ON RETURN :
  --   indices  indices to the points which span a fine mixed cell,
  --            it range must be 1..nbrpts.

  function Number_of_Points_in_Cell
              ( mix : Standard_Integer_Vectors.Vector ) 
              return integer32;

  -- DESCRIPTION :
  --   Returns the number of points in a fine mixed cell,
  --   according to the type of mixture defined by mix.

  procedure Parse_Cells
              ( file : in file_type; dim : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                cells : out Lists_of_Integer_Vectors.List;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Extracts the mixed cells from file.

  -- ON ENTRY :
  --   file     output file of demics, opened for output;
  --   dim      dimension of the problem;
  --   mix      type of mixture;
  --   verbose  flag to indicate for extra diagnostic output.

  -- ON RETURN :
  --   cells    indices to the points that span the mixed cells.

  procedure Process_Output
              ( dim : in integer32; filename : in string;
                mix : in Standard_Integer_Vectors.Vector;
                mv : out natural32;
                lif : out Standard_Floating_VecVecs.VecVec;
                cells : out Lists_of_Integer_Vectors.List;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Writes to screen the mixed volume on the output file
  --   in the string filename.

  -- ON ENTRY :
  --   dim      ambient dimension of the problem;
  --   filename is the name of the output file;
  --   mix      type of mixture;
  --   verbose  flag to indicate is more output is wanted.

  -- ON RETURN :  
  --   mv       the mixed volume computed by DEMiCs;
  --   lif      random lifting generated for the supports;
  --   cells    indices to the points that span the mixed cells.

end DEMiCs_Command_Line;
