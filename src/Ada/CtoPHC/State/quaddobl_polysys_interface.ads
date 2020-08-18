with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package QuadDobl_PolySys_Interface is

-- DESCRIPTION :
--   The functions below define the interface to the container of
--   polynomials in quad double precision.  The integer returned by all
--   functions should be zero if the job was successful,
--   otherwise the nonzero job code is returned.

  function QuadDobl_PolySys_Read
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system in quad double precision
  --   and initializes the container with the system provided by the user.
  --   The verbose level is given in vrblvl.

  function QuadDobl_PolySys_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Reads a system with quad double precision coefficients
  --   from file into the systems container.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the file name;
  --   b       the characters in the file name;
  --   vrblvl  the verbose level.

  function QuadDobl_PolySys_Write
             ( vrblvl : in integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the polynomial system stored in quad double precision
  --   to the defined output file or to screen.
  --   The verbose level is given in vrblvl.
 
  function QuadDobl_PolySys_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the dimension of the polynomial system stored
  --   in quad double precision.

  -- ON ENTRY :
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       the number of polynomials stored in quad double precision.

  function QuadDobl_PolySys_Set_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the number of polynomials in quad double precision
  --   and initializes the symbol table as well.

  -- ON ENTRY :
  --   a       the dimension of the polynomial system;
  --   vrblvl  the verbose level.

  function QuadDobl_PolySys_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of terms of a polynomial stored
  --   in quad double precision.

  -- ON ENTRY :
  --   a       a[1] is the index of a polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   a       a[0] is the number of terms of a polynomial 
  --           stored in quad double precision with index in a[1].

  function QuadDobl_PolySys_Degree
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the degree of a polynomial stored
  --   in double double precision.

  -- ON ENTRY :
  --   a       a[0] is the index of a polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       b[0] is the degree of a polynomial 
  --           stored in double double precision with index in a[0].

  function QuadDobl_PolySys_Get_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a term of a polynomial stored in quad double precision.

  -- ON ENTRY :
  --   a       a[1] contains the number of the polynomial,
  --           a[2] is the index to a term in that polynomial;
  --   vrblvl  the verbose level.

  -- ON RETURN :
  --   b       contains the exponent vector of the term;
  --   c       real and imaginary part of the coefficient.

  function QuadDobl_PolySys_Add_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Adds a term to a polynomial stored in quad double precision.

  -- ON ENTRY :
  --   a       a[1] contains the index of the polynomial;
  --   b       the exponents in the term;
  --   c       real and imaginary part of the coefficient;
  --   vrblvl  the verbose level.

  function QuadDobl_PolySys_String_Save
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Given the string representation of a polynomial,
  --   stores the polynomial in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the number of characters in the string,
  --           in a[1] is the number of variables in the polynomial,
  --           in a[2] is the index of the polynomial in the system;
  --   b       the string representation of the polynomial;
  --   vrblvl  is the verbose level.

  function QuadDobl_PolySys_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the size limit of the string representation of
  --   a polynomial in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the polynomial;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the size limit of the string representation of the
  --           polynomial in quad double precision with index a[0].

  function QuadDobl_PolySys_String_Load 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the string representation of a polynomial
  --   stored in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index of the polynomial;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the number of characters in b;
  --   b       the string representation of the polynomial.

  function QuadDobl_PolySys_Drop_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Drops a variable from the system stored in quad double precision.
  --   The variable to be dropped is given by its index.

  -- ON ENTRY :
  --   a      in a[0] is the index of the variable;
  --   vrblvl is the verbose level.

  function QuadDobl_PolySys_Drop_by_Name
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Drops a variable from the system 
  --   stored in double double precision.
  --   The variable to be dropped is given by its name.

  -- ON ENTRY :
  --   a      in a[0] is the number of characters in b;
  --   b      contains the name of the variable to be dropped;
  --   vrblvl is the verbose level.

  function QuadDobl_PolySys_Random_System
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes a random system with coefficients in quad double precision
  --   and stores the system.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the number of equations;
  --   b       in b[0] is the number of monomials per equation,
  --           in b[1] is the degree bound on the monomials,
  --           in b[2] is 0, 1, or 2, to set the type of coefficient:
  --           if 0, then complex coefficients on the unit circle,
  --           if 1, then all coefficients are equal to one,
  --           if 2, then real coefficients in [-1,+1].

  function QuadDobl_PolySys_Make_Homogeneous
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the system stored in quad double precision 1-homogeneous,
  --   eventually adding one random hyperplane.

  -- ON ENTRY :
  --   a       1 or 0, depending whether a random hyperplane
  --           needs to be added (if 1);
  --   vrblvl  is the verbose level.

  function QuadDobl_PolySys_Multi_Homogeneous
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Makes the system stored in quad double precision multihomogeneous,
  --   eventually adding one random hyperplane.

  -- ON ENTRY :
  --   a       in a[0] is the number of variables,
  --           in a[1] is the value of m, the size of the partition;
  --           if a[2] = 0, then random linear equations are added,
  --           if a[2] = 1, then euations of the for z0 - 1 = 0 are added;
  --   b       is the index representation of the partition;
  --   vrblvl  is the verbose level.

  function QuadDobl_PolySys_1Hom2Affine
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Turns the stored system in quad double precision
  --   from 1-homogeneous to affine coordinates.

  function QuadDobl_PolySys_mHom2Affine
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Turns the stored system in quad double precision
  --   from m-homogeneous to affine coordinates.

  -- ON ENTRY :
  --   a       in a[0] is the value for m;
  --   vrblvl  is the verbose level.

  function QuadDobl_PolySys_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears all polynomials stored in quad double precision.
  --   The verbose level is given in vrblvl.

  function QuadDobl_PolySys_Prompt_for_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Prompts the user for a target system in quad double precision.
  --   No solutions are read.  The verbose level is in vrblvl.

end QuadDobl_PolySys_Interface;
