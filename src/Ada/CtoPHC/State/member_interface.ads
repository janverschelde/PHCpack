with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Member_Interface is

-- DESCRIPTION :
--   The functions below decide whether a given point belongs to a
--   witness set, in double, double double, or quad double precision,
--   for polynomial or Laurent polynomial systems.
--   The test point is either given as a sequence of doubles,
--   or as a string representation of a solution in symbolic format,
--   with the symbols for the variables of the coordinates.

  function Member_Standard_Polynomial_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in double precision.
  --   The test point is given as a sequence of numbers.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --   b       in b[0] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in b[1] is the dimension of the witness set,
  --           in b[2] is the number of tasks, 0 for no tasking;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test,
  --           followed by 2*N doubles for the coordinates
  --           of the test point;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_DoblDobl_Polynomial_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in double double precision.
  --   The test point is given as a sequence of numbers.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --   b       in b[0] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in b[1] is the dimension of the witness set,
  --           in b[2] is the number of tasks, 0 for no tasking;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test,
  --           followed by 4*N doubles for the coordinates
  --           of the test point;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_QuadDobl_Polynomial_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in quad double precision.
  --   The test point is given as a sequence of numbers.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --   b       in b[0] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in b[1] is the dimension of the witness set,
  --           in b[2] is the number of tasks, 0 for no tasking;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test,
  --           followed by 8*N doubles for the coordinates
  --           of the test point;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_Standard_Laurent_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in double precision.
  --   The test point is given as sequence of numbers.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --   b       in b[0] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in b[1] is the dimension of the witness set,
  --           in b[2] is the number of tasks, 0 for no tasking;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test,
  --           followed by 2*N doubles for the coordinates
  --           of the test point;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_DoblDobl_Laurent_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in double double precision.
  --   The test point is given as a sequence of numbers.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --   b       in b[0] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in b[1] is the dimension of the witness set,
  --           in b[2] is the number of tasks, 0 for no tasking;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test,
  --           followed by 4*N doubles for the coordinates
  --           of the test point;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_QuadDobl_Laurent_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in quad double precision.
  --   The test point is given as a sequence of numbers.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --   b       in b[0] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in b[1] is the dimension of the witness set,
  --           in b[2] is the number of tasks, 0 for no tasking;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test,
  --           followed by 8*N doubles for the coordinates
  --           of the test point;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_Standard_Polynomial_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on polynomials in double precision.
  --   The test point is given as a string.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --           in a[1] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in a[2] is the dimension of the witness set,
  --           in a[3] is the number of characters in the string
  --           representation of the test point, as a solution;
  --           in a[4] is the number of tasks, 0 for no tasking;
  --   b       the characters, as many as the value of a[3],
  --           in the string representation of the test point;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_DoblDobl_Polynomial_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in double double precision.
  --   The test point is given as a string.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --           in a[1] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in a[2] is the dimension of the witness set,
  --           in a[3] is the number of characters in the string
  --           representation of the test point, as a solution;
  --           in a[4] is the number of tasks, 0 for no tasking;
  --   b       the characters, as many as the value of a[3],
  --           in the string representation of the test point;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_QuadDobl_Polynomial_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in quad double precision.
  --   The test point is given as a string.

  function Member_Standard_Laurent_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in double precision.
  --   The test point is given as a string.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --           in a[1] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in a[2] is the dimension of the witness set,
  --           in a[3] is the number of characters in the string
  --           representation of the test point, as a solution;
  --           in a[4] is the number of tasks, 0 for no tasking;
  --   b       the characters, as many as the value of a[3],
  --           in the string representation of the test point;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_DoblDobl_Laurent_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in double double precision.
  --   The test point is given as a string.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --           in a[1] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in a[2] is the dimension of the witness set,
  --           in a[3] is the number of characters in the string
  --           representation of the test point, as a solution;
  --           in a[4] is the number of tasks, 0 for no tasking;
  --   b       the characters, as many as the value of a[3],
  --           in the string representation of the test point;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

  function Member_QuadDobl_Laurent_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in quad double precision.
  --   The test point is given as a string.

  -- ON ENTRY :
  --   a       in a[0] is the verbose flag;
  --           in a[1] is the dimension N of the test point,
  --           equal to its number of complex coordinates,
  --           in a[2] is the dimension of the witness set,
  --           in a[3] is the number of characters in the string
  --           representation of the test point, as a solution;
  --           in a[4] is the number of tasks, 0 for no tasking;
  --   b       the characters, as many as the value of a[3],
  --           in the string representation of the test point;
  --   c       in c[0] is the tolerance on the residual for evaluation,
  --           in c[1] is the tolerance on the membership test;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the result of the evaluation test, 0 for fail, 1 for success;
  --   b       the result of the membership test, 0 for fail, 1 for success.

end Member_Interface;
