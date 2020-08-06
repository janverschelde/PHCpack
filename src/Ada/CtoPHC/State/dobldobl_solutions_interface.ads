with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package DoblDobl_Solutions_Interface is

-- DESCRIPTION :
--   The functions below define the interface to the container of
--   solutions in double double precision.  The integer returned by all
--   functions should be zero if the job was successful,
--   otherwise the nonzero job code is returned.

  function DoblDobl_Solutions_Read 
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRPTION :
  --   Prompts the user for solutions in double double precision
  --   and initializes the container with the user input.
  --   The verbose level is given in vrblvl.

  function DoblDobl_Solutions_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION 
  --   Reads solutions in double double precision from file.

  -- ON ENTRY :
  --   a       the number of characters in the file name;
  --   b       the file name to read the solutions from;
  --   vrblvl  is the verbose level.

  function DoblDobl_System_Solutions_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION 
  --   Reads a system with solutions in double double precision from file.

  -- ON ENTRY :
  --   a       the number of characters in the file name;
  --   b       the file name to read the system and solutions from;
  --   vrblvl  is the verbose level.

  function DoblDobl_Solutions_Write
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the solutions in double double precision
  --   to the defined output file or to screen.

  function DoblDobl_Solutions_Size 
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of solutions in double double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the number of solutions in double double precision.

  function DoblDobl_Solutions_Dimension
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the solutions in double double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the dimension of the solutions in double double precision.

  function DoblDobl_Solutions_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
           
  -- DESCRIPTION :
  --   Retrieves a solution in double double precision.

  -- ON ENTRY :
  --   a       index of the solution to be retrieved;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       integer attributes of the solution;
  --   c       double attributes of the solutions.

  function DoblDobl_Solutions_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
           
  -- DESCRIPTION :
  --   Updates a solution in double double precision.

  -- ON ENTRY :
  --   a       index of the solution to be retrieved;
  --   b       integer attributes of the solution;
  --   c       double attributes of the solutions.
  --   vrblvl  the verbose level.

  function DoblDobl_Solutions_Add
             ( b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;
           
  -- DESCRIPTION :
  --   Appends a solution in double double precision.

  -- ON ENTRY :
  --   b       integer attributes of the solution;
  --   c       double attributes of the solutions.
  --   vrblvl  the verbose level.

  function DoblDobl_Solutions_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the string representation of
  --   a solution stored in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of characters in the string representation
  --           of the solution with the given index.

  function DoblDobl_Solutions_Get_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation of
  --   a solution stored in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution,
  --           in a[1] is the number of characters in the string
  --           representation of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the string representation of the solution 
  --           with the given index.

  function DoblDobl_Solutions_Add_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Appends a solution given in its string representation
  --   to the list of solutions stored in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of characters in the string b;
  --   b       the string representation of a solution;
  --   vrblvl  is the verbose level.

  function DoblDobl_Solutions_Retrieve_Next
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the pointer to the start of the list of solutions in
  --   double double precision or returns the next solution.

  -- ON ENTRY :
  --   a       if a[0] = 0, then the pointer with be initialized,
  --           otherwise the next solution is returned;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       if a[0] /= 0, the a[0] contains the index of the
  --           returned solution;
  --   b       the multiplicity of the next solution;
  --   c       4*n+10 doubles, first for the two double doubles for t,
  --           then 4*n doubles for the coefficients of the solution vector,
  --           followed by one double for forward error,
  --           one double for the estimate of the inverse condition number,
  --           and one double for the backward error.

  function DoblDobl_Solutions_Move_Pointer
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Moves the pointer in the list to the next solution
  --   stored in double double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the index of the current solution,
  --           after the moving of the pointer.

  function DoblDobl_Solutions_Current_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the string representation
  --   of the current solution stored in double double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the index of the current solution;
  --   b       in b[0] is the size of the string representation.

  function DoblDobl_Solutions_Current_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation
  --   of the current solution stored in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the length of the string representation,
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the string representation of the current solution.

  function DoblDobl_Solutions_Drop_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Drops from the solutions stored in double double precision
  --   one of the coordinates.
  --   The coordinate to be dropped is given by its index.

  -- ON ENTRY :
  --   a       in a[0] is the index of the coordinated to be dropped;
  --   vrblvl  is the verbose level.

  function DoblDobl_Solutions_Drop_by_Name
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Drops from the solutions stored in double double precision
  --   one of the coordinates.
  --   The coordinate to be dropped is given by its name.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       is the name of the coordinate to be dropped;
  --   vrblvl  is the verbose level.

  function DoblDobl_Solutions_Make_Homogeneous
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the solutions store in double double precision 1-homogeous,
  --   adding one extra coordinate equal to one to every solution.

  function DoblDobl_Solutions_Multi_Homogeneous
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the solutions store in double double precision m-homogeous,
  --   adding m extra coordinates equal to one to every solution.

  -- ON ENTRY :
  --   a       in a[0] is the value for m;
  --   vrblvl  is the verbose level.

  function DoblDobl_Solutions_1Hom2Affine
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Turns the solutions stored in double double precision
  --   from 1-homogeneous to affine coordinates.
  --   The verbose level is given in vrblvl.

  function DoblDobl_Solutions_mHom2Affine
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Turns the solutions stored in double double precision
  --   from m-homogeneous to affine coordinates.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of sets in the partition;
  --   b       contains the index representation of the partition;
  --   vrblvl  is the verbose level.

  function DoblDobl_Solutions_Tzero
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the continuation parameter t of all solutions in double
  --   double precision to zero.  The verbose level is in vrblvl.

  function DoblDobl_Solutions_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates all solutions in double double precision.

end DoblDobl_Solutions_Interface;
