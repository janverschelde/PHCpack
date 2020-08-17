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

  function Standard_Solutions_Add_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Appends a solution given in its string representation
  --   to the list of solutions stored in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of characters in the string b;
  --   b       the string representation of a solution;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Replace_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Replaces a solution stored in double precision by
  --   the solution given in its string representation.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution to be replace,
  --           in a[1] is the number of variables,
  --           in a[2] is the number of characters in the string b;
  --   b       the string representation of a solution;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Intro_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of characters of the string representation
  --           of the intro to the solution in double precision.

  function Standard_Solutions_Intro_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation of the intro
  --   to a solution in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution,
  --           in a[1] is the number of characters in the string
  --           representation of the diagnostics of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the string representation of the intro
  --           to the solution in double precision.

  function Standard_Solutions_Vector_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the string representation of all coordinates 
  --   of a solution in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of characters of the string representation
  --           of all coordinates of the solution in double precision.

  function Standard_Solutions_Vector_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation of all coordinates
  --   of a solution in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution,
  --           in a[1] is the number of characters in the string
  --           representation of the diagnostics of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the string representation of all coordinates
  --           of the solution in double precision.

  function Standard_Solutions_Diagnostics_String_Size
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the diagnostics string of a solution
  --   stored in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the number of characters in the diagnostics string.

  function Standard_Solutions_Diagnostics_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the diagnostics string of a solution in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution,
  --           in a[1] is the number of characters in the string
  --           representation of the diagnostics of the solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the string representation of the diagnostics
  --           of a solution in double precision.

  function Standard_Solutions_Retrieve_Next
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the pointer to the start of the list of solutions in
  --   double precision or returns the next solution.

  -- ON ENTRY :
  --   a       if a[0] = 0, then the pointer with be initialized,
  --           otherwise the next solution is returned;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       if a[0] /= 0, the a[0] contains the index of the
  --           returned solution;
  --   b       the multiplicity of the next solution;
  --   c       2*n+5 doubles, first for the two double for t,
  --           then 2*n doubles for the coefficients of the solution vector,
  --           followed by one double for forward error,
  --           one double for the estimate of the inverse condition number,
  --           and one double for the backward error.

  function Standard_Solutions_Move_Pointer
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Moves the pointer in the list to the next solution
  --   stored in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the index of the current solution,
  --           after the moving of the pointer.

  function Standard_Solutions_Current_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size of the string representation
  --   of the current solution stored in double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the index of the current solution;
  --   b       in b[0] is the size of the string representation.

  function Standard_Solutions_Current_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation
  --   of the current solution stored in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the length of the string representation,
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       the string representation of the current solution.

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

  function Standard_Solutions_Tzero
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the continuation parameter t of all solutions in double
  --   precision to zero.  The verbose level is in vrblvl.

  function Standard_Solutions_Read_Next
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads the next solution in double precision from file.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the multiplicity of the solution;
  --   c       2*n + 5 doubles, starting with two doubles for t,
  --           2*n doubles for the solution vectors,
  --           one double for the forward error,
  --           one double for the estimate for the inverse condition number,
  --           and one double for the backward error.

  function Standard_Solutions_Write_Next
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the next solution in double precision to file.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables;
  --   b       the multiplicity of the solution;
  --   c       2*n + 5 doubles, starting with two doubles for t,
  --           2*n doubles for the solution vectors,
  --           one double for the forward error,
  --           one double for the estimate for the inverse condition number,
  --           and one double for the backward error;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Next_to_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the next solution in double precision
  --   to the defined output file.

  -- ON ENTRY :
  --   a       in a[0] is the index of the solution;
  --   b       the multiplicity of the solution;
  --   c       2*n + 5 doubles, starting with two doubles for t,
  --           2*n doubles for the solution vectors,
  --           one double for the forward error,
  --           one double for the estimate for the inverse condition number,
  --           and one double for the backward error;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the count of the solutions written to file.

  function Standard_Solutions_Total_Degree
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the next total degree start solution in double precision.
  --   This start solution is the solution of a total degree start system.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the index of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the multiplicity of the solution;
  --   c       2*n + 5 doubles, starting with two doubles for t,
  --           2*n doubles for the solution vectors,
  --           one double for the forward error,
  --           one double for the estimate for the inverse condition number,
  --           and one double for the backward error.

  function Standard_Solutions_Next_Product
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the next linear-product start solution in double precision.
  --   This start solution is the solution of a linear-product start system.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the index of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the multiplicity of the solution;
  --   c       2*n + 5 doubles, starting with two doubles for t,
  --           2*n doubles for the solution vectors,
  --           one double for the forward error,
  --           one double for the estimate for the inverse condition number,
  --           and one double for the backward error.

  function Standard_Solutions_Lex_Product
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the linear-product start solution in double precision,
  --   indexed by the lexicographic order defined by the given index.
  --   This start solution is the solution of a linear-product start system.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the index of the start solution;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the multiplicity of the solution;
  --   c       2*n + 5 doubles, starting with two doubles for t,
  --           2*n doubles for the solution vectors,
  --           one double for the forward error,
  --           one double for the estimate for the inverse condition number,
  --           and one double for the backward error.

  function Standard_Solutions_Next_Witness
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads the next witness point solution from file in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the witness set,
  --           in a[1] is the number of variables.

  -- ON RETURN :
  --   b       the multiplicity of the solution;
  --   c       2*n + 5 doubles, starting with two doubles for t,
  --           2*n doubles for the solution vectors,
  --           one double for the forward error,
  --           one double for the estimate for the inverse condition number,
  --           and one double for the backward error;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Scan_Banner
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Scans the input file for the SOLUTIONS banner and returns 0
  --   if the banner is found, otherwise the job code is returned.
  --   The verbose level is given by the value of vrblvl.

  function Standard_Solutions_Read_Dimensions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads the length of the solution list and the dimension
  --   of the solution vectors from the input file.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the length of the solution list;
  --   b       in b[0] is the dimension of each solution vector.

  function Standard_Solutions_Write_Dimensions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the length of the solution list and the dimension
  --   of the solution vectors to the output file.

  -- ON ENTRY :
  --   a       in a[0] is the length of the solution list;
  --   b       in b[0] is the dimension of each solution vector;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Close_Input_File
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Closes the input file.

  -- ON ENTRY :
  --   a       in a[0] is the index of the file,
  --           if 0, then there is only one input file,
  --           otherwise, there may be many input files;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Close_Output_File
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Closes the output file.
  --   The verbose level is given by the value of vrblvl.

  function Standard_Solutions_Banner_to_Output
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the SOLUTIONS banner to the defined output file,
  --   or if the output is undefined, to standard output.

  function Standard_Solutions_Dimensions_to_Output
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the length of the solution list and the dimension
  --   of the solution vectors to the defined output file,
  --   or if undefined, to the standard output.

  -- ON ENTRY :
  --   a       in a[0] is the length of the solution list;
  --   b       in b[0] is the dimension of each solution vector;
  --   vrblvl  is the verbose level.

  function Standard_Solutions_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Deallocates all solutions in double precision.

end Standard_Solutions_Interface;
