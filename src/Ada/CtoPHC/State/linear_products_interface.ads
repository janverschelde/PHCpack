with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Linear_Products_Interface is

-- DESCRIPTION :
--   The functions below interface to the linear-product start systems,
--   based on set structures or on partitions of the set of variables.

  function Linear_Products_Structure_Make
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a supporting set structure for the polynomial system
  --   stored in double precision.
  --   The value of the verbose level is in vrblvl.
  
  function Linear_Products_Structure_Write
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes to screen the set structure made with the function
  --   Linear_Products_Structure_Make.

  function Linear_Products_Structure_Bound
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the set structure bound.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the bound based on the set structure.

  function Linear_Products_System_Make
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a random linear-product system based
  --   on the defined set structure.

  function Linear_Products_System_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a linear-product system from file, for given file name.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string b;
  --   b       contains the name of the file;
  --   vrblvl  is the verbose level.

  function Linear_Products_System_Solve
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Solves the define linear-product system and stores the
  --   solutions in double precision.

  function Linear_Products_Structure_String_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation of the structure
  --   to define a linear-product system.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] are the number of characters in the string;
  --   b       contains the string representation of the structure.

  function Linear_Products_Structure_String_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets the structure to define a linear-product system
  --   via the given string representation.

  -- ON ENTRY :
  --   a       in a[0] are the number of characters in the string;
  --   b       the string representation of the structure;
  --   vrblvl  is the verbose level.

  function Linear_Products_Structure_Check
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Checks if the defined structure supports the polynomial system
  --   stored in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is 1 for true, 0 for false.

  function Linear_Products_Partition_Make
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a partition for the set of variables of the polynomial
  --   system stored in double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the number of characters in the string b;
  --   b       the string representation for the partition.

  function Linear_Products_Partition_Bound
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the multi-homogeneous Bezout bound for the polynomial
  --   system stored in double precision and for the given partition.

  -- ON ENTRY :
  --   a       the number of characters in the string b;
  --   b       the string representation of a partition;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the multi-homogeneous Bezout bound for the polynomial
  --           system in double precision, for the given partition.

  function Linear_Products_Partition_System
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a linear-product system for the polynomial system stored
  --   in double precision and for the given partition.
  --   The system stored in double precision is replaced by
  --   the linear-product system.
  --   To solve the linear-product system, apply the function
  --   Linear_Products_System_Solve.

  -- ON ENTRY :
  --   a       the number of characters in the string b;
  --   b       the string representation of a partition;
  --   vrblvl  the verbose level.

  function Linear_Products_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the set structure and also all data
  --   used to define a linear product system.

end Linear_Products_Interface;
