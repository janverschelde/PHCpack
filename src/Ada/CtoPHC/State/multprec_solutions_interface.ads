with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Multprec_Solutions_Interface is

-- DESCRIPTION :
--   The functions below define the interface to the container of
--   solutions in multiprecision.  The integer returned by all
--   functions should be zero if the job was successful,
--   otherwise the nonzero job code is returned.

  function Multprec_Solutions_Read 
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRPTION :
  --   Prompts the user for solutions in multiprecision
  --   and initializes the container with the user input.
  --   The verbose level is given in vrblvl.

  function Multprec_System_Solutions_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION 
  --   Reads a system with solutions in quad double precision from file.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the file name,
  --           in a[1] is the number of decimal places;
  --   b       the file name to read the system and solutions;
  --   vrblvl  is the verbose level.

  function Multprec_Solutions_Write
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the solutions in multiprecision
  --   to the defined output file or to screen.

  function Multprec_Solutions_Write_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the solutions of the start system in multiprecision
  --   to the defined output file or to screen.

  function Multprec_Solutions_Write_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the solutions of the target system in multiprecision
  --   to the defined output file or to screen.

  function Multprec_Solutions_Size 
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of solutions in multiprecision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the number of solutions in multiprecision.

  function Multprec_Solutions_Dimension
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the solutions in multiprecision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the dimension of the solutions in multiprecision.

  function Multprec_Solutions_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the string representation of
  --   a solution stored in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of characters in the string representation
  --           of the solution with the given index.

  function Multprec_Solutions_Get_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation of
  --   a solution stored in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution,
  --           in a[1] is the number of characters in the string
  --           representation of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the string representation of the solution 
  --           with the given index.

  function Multprec_Solutions_Add_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Appends a solution given in its string representation
  --   to the list of solutions stored in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of characters in the string b;
  --   b       the string representation of a solution;
  --   vrblvl  is the verbose level.

  function Multprec_Solutions_Set_Pointer
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the pointer to the start of the list
  --   of solutions stored in multiprecision.

  function Multprec_Solutions_Move_Pointer
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Moves the pointer in the list to the next solution
  --   stored in multiprecision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the index of the current solution,
  --           after the moving of the pointer.

  function Multprec_Solutions_Current_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the string representation
  --   of the current solution stored in multiprecision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the index of the current solution;
  --   b       in b[0] is the size of the string representation.

  function Multprec_Solutions_Current_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation
  --   of the current solution stored in multiprecision.

  -- ON ENTRY :
  --   a       in a[0] is the length of the string representation,
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the string representation of the current solution.

  function Multprec_Solutions_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates all solutions in multiprecision.

end Multprec_Solutions_Interface;
