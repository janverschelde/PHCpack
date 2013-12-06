with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Drivers_for_Solution_Filters is

-- DESCRIPTION :
--   This package provides some utilities and the main driver to
--   filter lists of solutions, subject to user given criteria.

  procedure Write_Symbols;

  -- DESCRIPTION :
  --   Writes the current content of the symbol table to screen.

  function Prompt_Symbol return natural32;

  -- DESCRIPTION :
  --   Prompts the user for a symbol and returns 
  --   the position of the symbol in the symbol table.

  procedure Read_Double_Complex ( c : out Complex_Number );

  -- DESCRIPTION :
  --   Prompts the user to give the real and imaginary part
  --   of a complex number, which will be returned in c.

  function Show_Menu_and_Prompt_Answer return character;

  -- DESCRIPTION :
  --   Shows the user the menu with available options and
  --   returns the selected answer.

  procedure Driver_for_Solution_Filters
               ( file : in file_type; sols : in Solution_List );

  -- DESCRIPTION :
  --   Creates new solution lists from filtering a given list subject
  --   to the following exclusive criteria :
  --     1) vanishing or spurious;
  --     2) zero or free k-th component;
  --     3) regular or singular;
  --     4) reached a target value.
  --   The criteria are set up with respect to a given tolerance.

  procedure Driver_for_Solution_Filters
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

end Drivers_for_Solution_Filters;
