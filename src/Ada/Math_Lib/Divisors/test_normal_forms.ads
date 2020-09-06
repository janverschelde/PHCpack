with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Multprec_Integer_Matrices;

package Test_Normal_Forms is

-- DESCRIPTION :
--   Tests Hermite and Smith Normal Form of an integer matrix.

  procedure Test_Hermite_Normal_Form 
              ( mat : in Standard_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Test on the Hermite normal form of a matrix.

  procedure Interactive_Test_Hermite_Normal_Form ( n,m : in integer32 );

  -- DESCRIPTION :
  --   Prompts for an n-by-m matrix and tests the Hermite normal form.

  procedure Test_Standard32_Smith_Normal_Form
              ( mat : in Standard_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean );

  -- DESCRIPTION :
  --   Tests the Smith normal form of a matrix,
  --   using standard 32-bit integer arithmetic.
  --   If output is switched on, then all intermediate information is printed,
  --   otherwise this is only the case when a bug is detected;

  procedure Test_Standard64_Smith_Normal_Form
              ( mat : in Standard_Integer64_Matrices.Matrix;
                output : in boolean; bug : out boolean );

  -- DESCRIPTION :
  --   Tests the Smith normal form of a matrix,
  --   using standard 64-bit integer arithmetic.
  --   If output is switched on, then all intermediate information is printed,
  --   otherwise this is only the case when a bug is detected;

  procedure Fix_Empty_Entries
              ( mat : in out Multprec_Integer_Matrices.Matrix;
                output : in boolean );

  -- DESCRIPTION :
  --   Sets every empty entry of the matrix to zero.
  --   Also -0 is turned into 0.

  procedure Test_Multprec_Smith_Normal_Form
              ( mat : in Multprec_Integer_Matrices.Matrix;
                output : in boolean; bug : out boolean );

  -- DESCRIPTION :
  --   Tests the Smith normal form of a matrix,
  --   using multiprecision integer arithmetic.
  --   If output is switched on, then all intermediate information is printed,
  --   otherwise this is only the case when a bug is detected;

  procedure Interactive_Test_Standard_Smith_Normal_Form
              ( n,m : in integer32 );

  -- DESCRIPTION :
  --   Tests the Smith normal form on a given n-by-m matrix,
  --   using standard 32-bit integer arithmetic.

  procedure Interactive_Test_Multprec_Smith_Normal_Form
              ( n,m : in integer32 );

  -- DESCRIPTION :
  --   Tests the Smith normal form on a given n-by-m matrix,
  --   using multiprecision integer arithmetic.

  procedure Random_Test_Standard32_Smith_Normal_Form ( n,m : in integer32 );

  -- DESCRIPTION :
  --   Tests the Smith normal form on randomly generated n-by-m matrices,
  --   using standard 32-bit integer arithmetic.

  procedure Random_Test_Standard64_Smith_Normal_Form ( n,m : in integer32 );

  -- DESCRIPTION :
  --   Tests the Smith normal form on randomly generated n-by-m matrices,
  --   using standard 32-bit integer arithmetic.

  procedure Random_Test_Standard_Smith_Normal_Form ( n,m : in integer32 );

  -- DESCRIPTION :
  --   Tests the Smith normal form on randomly generated n-by-m matrices,
  --   using standard 32-bit or 64-bit integer arithmetic.

  procedure Random_Test_Multprec_Smith_Normal_Form ( n,m : in integer32 );

  -- DESCRIPTION :
  --   Tests the Smith normal form on randomly generated n-by-m matrices,
  --   using multiprecision integer arithmetic.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Normal_Forms;
