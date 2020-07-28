with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Pieri_Interface is

-- DESCRIPTION :
--   The functions below give access to Pieri homotopies.

  function Pieri_Write_Menu ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the menu of the available jobs in the interface.
  --   The verbose level is given in the value of vrblvl.

  function Pieri_Initialize_Dimensions
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the dimension with the values of m, p, and q.

  -- ON ENTRY :
  --   a       in a[0] is the value for m,
  --           in a[1] is the value for p,
  --           in a[2] is the value for q;
  --   vrblvl  is the verbose level.

  function Pieri_Initialize_Input_Planes
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the input planes for the Pieri homotopy.

  -- ON ENTRY :
  --   a       in a[0] is the value for m,
  --           in a[1] is the value for p,
  --           in a[2] is the value for q,
  --   b       the value m*p + q*(m+p);
  --   c       coefficients, real and imaginary parts,
  --           of the complex input planes, stored as row wise matrices;
  --   vrblvl  is the verbose level.

  function Pieri_Initialize_Interpolation_Points
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes the interpolation points for the Pieri homotopy.

  -- ON ENTRY :
  --   a       in a[0] is the value for m,
  --           in a[1] is the value for p,
  --           in a[2] is the value for q,
  --   b       the value m*p + q*(m+p);
  --   c       coefficients of the interpolation points,
  --           as real and imaginary parts of complex numbers;
  --   vrblvl  is the verbose level.

  function Pieri_Store_Start_Pivots
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Stores the pivots at the start.

  -- ON ENTRY :
  --   a       in a[0] is the value for p;
  --   b       the 2*p values of the top and bottom pivots;
  --   vrblvl  is the verbose level.

  function Pieri_Store_Target_Pivots 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Stores the pivots at the target.

  -- ON ENTRY :
  --   a       in a[0] is the value for p;
  --   b       the 2*p values of the top and bottom pivots;
  --   vrblvl  is the verbose level.

  function Pieri_Store_Start_Coefficients
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Stores the coefficients of a start solution.

  -- ON ENTRY :
  --   a       in a[0] is the degree of freedom;
  --   c       coefficients of a complex vectors in pairs of doubles,
  --           sequences of real and imaginary parts;
  --   vrblvl  is the verbose level.

  function Pieri_Get_Target_Solution
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a target solution.

  -- ON ENTRY :
  --   a       in a[0] is the degree of freedom;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coefficients of the target solution.

  function Pieri_Silent_Track
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one path without intermediate output.

  function Pieri_Report_Track
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks one path with intermediate output to screen.

  function Pieri_Silent_Verify
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Verifies the intersection conditions with no extra output.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       the sum of all residuals.
 
  function Pieri_Report_Verify
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Verifies the intersection conditions with extra output.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       the sum of all residuals.
 
  function Pieri_Root_Count
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Computes the Pieri root count.

  -- ON ENTRY :
  --   a       in a[0] is the value for m,
  --           in a[1] is the value for p,
  --           in a[2] is the value for q;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the root count.

  function Pieri_Localization_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTON :
  --   Return the location poset string.

  -- ON ENTRY :
  --   a       in a[0] is the value for m,
  --           in a[1] is the value for p,
  --           in a[2] is the value for q;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] are the characters in b;
  --   b       the characters in the string.

  function Pieri_Run_Homotopies
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the Pieri homotopies.

  -- ON ENTRY :
  --   a       in a[0] is the value for m,
  --           in a[1] is the value for p,
  --           in a[2] is the value for q;
  --   c       coefficients of the input planes;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       the Pieri root count.

  function Pieri_Real_Osculating_Planes
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Generates n real m-planes in d-space osculating a rational curve
  --   at n given points.

  -- ON ENTRY :
  --   a       in a[0] is the value for m,
  --           in a[1] is the value for p,
  --           in a[2] is the value for q;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       coefficients of the planes.

  function Pieri_Make_Target_System 
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Gets the dimensions and coefficients of the m-planes
  --   and puts the target system in the container.


  -- ON ENTRY :
  --   a       in a[0] is the value for m,
  --           in a[1] is the value for p,
  --           in a[2] is the value for q;
  --   c       coefficients of the planes.
  --   vrblvl  is the verbose level.

  function Pieri_Clear ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the Pieri homotopy.

end Pieri_Interface;
