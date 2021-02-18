with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Main_Solution_Filters is

-- DESCRIPTION :
--   This package provides some utilities and the main procedures
--   to filter lists of solutions, subject to user given criteria.
--   The drivers support double, double double, and quad double precision
--   and work on solution lists that are entirely in main memory
--   and on solution lists that are directly read from file.

  procedure Write_Symbols;

  -- DESCRIPTION :
  --   Writes the current content of the symbol table to screen.

  function Prompt_Symbol return natural32;

  -- DESCRIPTION :
  --   Prompts the user for a symbol and returns 
  --   the position of the symbol in the symbol table.

  procedure Read_Double_Double ( f : out double_double );
  procedure Read_Quad_Double ( f : out quad_double );

  -- DESCRIPTION :
  --   Reads a double double or a quad double.
  --   In case of an exception, the user is invited to try again.

  procedure Read_Double_Complex
              ( c : out Standard_Complex_Numbers.Complex_Number );
  procedure Read_DoblDobl_Complex
              ( c : out DoblDobl_Complex_Numbers.Complex_Number );
  procedure Read_QuadDobl_Complex
              ( c : out QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Prompts the user to give the real and imaginary part
  --   of a complex number, which will be returned in c.

  function Show_Menu_and_Prompt_Answer return character;

  -- DESCRIPTION :
  --   Shows the user the menu with available options and
  --   returns the selected answer.

  procedure Main_Solution_Filter
               ( file : in file_type;
                 sols : in Standard_Complex_Solutions.Solution_List );
  procedure Main_Solution_Filter
               ( file : in file_type;
                 sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Main_Solution_Filter
               ( file : in file_type;
                 sols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Creates new solution lists from filtering a given list subject
  --   to the following exclusive criteria :
  --     1) vanishing or spurious;
  --     2) zero or free k-th component;
  --     3) regular or singular;
  --     4) reached a target value.
  --   The criteria are set up with respect to a given tolerance.

  procedure Standard_Solution_Filter
               ( infile,outfile : in file_type; len,dim : in natural32 );
  procedure DoblDobl_Solution_Filter
               ( infile,outfile : in file_type; len,dim : in natural32 );
  procedure QuadDobl_Solution_Filter
               ( infile,outfile : in file_type; len,dim : in natural32 );

  -- DESCRIPTION :
  --   Reads solutions from infile and writes a selection of those
  --   solutions to outfile, subject to user given criteria.

  -- ON ENTRY :
  --   infile    input file, must be properly positioned right after
  --             reading the dimension and length of the solution list;
  --   outfile   output file to write the selection to;
  --   len       length of the solution list;
  --   dim       length of the vectors in the solution list.

  procedure Read_Dimensions
              ( infile : out file_type; len,dim : out natural32 );

  -- DESCRIPTION :
  --   Prompts the user for the name of a file for the solutions,
  --   eventually preceeded by a polynomial system.  Only the dimensions 
  --   of the solution list are read, which is fine for huge lists.

  procedure Main ( infilename,outfilename : in string );

  -- DESCRIPTION :
  --   Filters solution lists on file, subject to criteria selected
  --   by the user.
  --   The arguments are the names of the input and output files.

end Main_Solution_Filters;
