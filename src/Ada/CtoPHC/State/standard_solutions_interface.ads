with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Standard_Solutions_Interface is

-- DESCRIPTION :
--   The functions below define the interface to the container of
--   solutions in double precision.  The integer returned by all
--   functions should be zero if the job was successful,
--   otherwise the nonzero job code is returned.

  function Standard_Solutions_Read 
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRPTION :
  --   Prompts the user for solutions in double precision
  --   and initializes the container with the user input.
  --   The verbose level is given in vrblvl.

  function Standard_Solutions_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION 
  --   Reads solutions in double precision from file.

  -- ON ENTRY :
  --   a       the number of characters in the file name;
  --   b       the file name to read the solutions from;
  --   vrblvl  is the verbose level.

  function Standard_System_Solutions_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION 
  --   Reads a system with solutions in double precision from file.

  -- ON ENTRY :
  --   a       the number of characters in the file name;
  --   b       the file name to read the system and solutions from;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Write
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the solutions in double precision
  --   to the defined output file or to screen.

  function Standard_Solutions_Size 
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of solutions in double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the number of solutions in double precision.

  function Standard_Solutions_Dimension
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the solutions in double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the dimension of the solutions in double precision.

  function Standard_Solutions_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
           
  -- DESCRIPTION :
  --   Retrieves a solution in double precision.

  -- ON ENTRY :
  --   a       index of the solution to be retrieved;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       integer attributes of the solution;
  --   c       double attributes of the solutions.

  function Standard_Solutions_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
           
  -- DESCRIPTION :
  --   Updates a solution in double precision.

  -- ON ENTRY :
  --   a       index of the solution to be retrieved;
  --   b       integer attributes of the solution;
  --   c       double attributes of the solutions.
  --   vrblvl  the verbose level.

  function Standard_Solutions_Add
             ( b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
           
  -- DESCRIPTION :
  --   Appends a solution in double precision.

  -- ON ENTRY :
  --   b       integer attributes of the solution;
  --   c       double attributes of the solutions.
  --   vrblvl  the verbose level.

  function Standard_Solutions_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the string representation of
  --   a solution stored in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of characters in the string representation
  --           of the solution with the given index.

  function Standard_Solutions_Get_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation of
  --   a solution stored in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution,
  --           in a[1] is the number of characters in the string
  --           representation of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the string representation of the solution 
  --           with the given index.

  function Standard_Solutions_Drop_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Drops from the solutions stored in double precision
  --   one of the coordinates.
  --   The coordinate to be dropped is given by its index.

  -- ON ENTRY :
  --   a       in a[0] is the index of the coordinated to be dropped;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Drop_by_Name
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Drops from the solutions stored in double precision
  --   one of the coordinates.
  --   The coordinate to be dropped is given by its name.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       is the name of the coordinate to be dropped;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Make_Homogeneous
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the solutions store in double precision 1-homogeous,
  --   adding one extra coordinate equal to one to every solution.

  function Standard_Solutions_Multi_Homogeneous
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the solutions store in double precision m-homogeous,
  --   adding m extra coordinates equal to one to every solution.

  -- ON ENTRY :
  --   a       in a[0] is the value for m;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_1Hom2Affine
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Turns the solutions stored in double precision
  --   from 1-homogeneous to affine coordinates.
  --   The verbose level is given in vrblvl.

  function Standard_Solutions_mHom2Affine
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Turns the solutions stored in double precision
  --   from m-homogeneous to affine coordinates.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of sets in the partition;
  --   b       contains the index representation of the partition;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates all solutions in double precision.

end Standard_Solutions_Interface;
