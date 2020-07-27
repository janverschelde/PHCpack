with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Sweep_Interface is

-- DESCRIPTION :
--   The functions define the interface to run a sweep,
--   defined by a real-parameter or a complex-convex homotopy.

  function Sweep_Define_Parameters_Numerically
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Defines the parameters in the sweep numerically.

  -- ON ENTRY :
  --   a       in a[0] is the number of equations,
  --           in a[1] is the total number of variables,
  --           in a[2] is the number of parameters;
  --   b       in b is a sequence of a[2] integers in range 1..a[1],
  --           with the indices of the variables that are parameters;
  --   vrblvl  is the verbose level.

  function Sweep_Define_Parameters_Symbolically
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Defines the parameters in the sweep symbolically.

  -- ON ENTRY :
  --   a       in a[0] is the number of equations,
  --           in a[1] is the total number of variables,
  --           in a[2] is the number of parameters,
  --           in a[3] is the number of characters stored in b;
  --   b       the symbols that define the parameters are separated
  --           by one space;
  --   vrblvl  is the verbose level.

  function Sweep_Number_of_Equations
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of equations in the parameter homotopy.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of equations in the parameter homotopy.

  function Sweep_Number_of_Variables
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of variables in the parameter homotopy.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of variables in the parameter homotopy.

  function Sweep_Number_of_Parameters
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of parameters in the parameter homotopy.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of parameters in the parameter homotopy.

  function Sweep_Get_Parameters_Numerically
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the numerical definitions of the parameters.

  -- ON ENTRY :
  --   a       space allocated for the indices of the parameters;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the indices that define which variables are parameters.

  function Sweep_Get_Parameters_Symbolically
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the symbolic definitions of the parameters.

  -- ON ENTRY :
  --   b       space allocated for a string of symbols;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       the number of characters in b;
  --   b       symbols for the parameters each separated by one space.

  function Sweep_Parameters_Clear
             ( vrblvl : integer32 := 0 ) return integer32;
  
  -- DESCRIPTION :
  --   Clears the parameter definitions.

  function Sweep_Set_Parameter_Values
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets start or target parameter values,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively
  --           for double, double double, or quad double precision;
  --           in a[1] is 0 or 1, if values for start or target respectively;
  --   b       the number of coefficients to store the values of
  --           the parameters, as real and imaginary parts;
  --   c       are the values for the parameters;
  --   vrblvl  is the verbose level.

  function Sweep_Get_Parameter_Values
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Get start or target parameter values,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively
  --           for double, double double, or quad double precision;
  --           in a[1] is 0 or 1, if values for start or target respectively;
  --   c       space allocated for the parameter values;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       the values for the parameters;

  function Sweep_Complex_Convex_Parameter
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a complex sweep with a convex combination between parameters,
  --   in double, double double, or quad double precision.

  -- DESCRIPTION :
  --   Runs a real sweep on a natural-parameter homotopy,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2 for the precision level,
  --           respectively double, double double, or quad double,
  --           in a[1] is the parameter for random gamma,
  --           if a[1] = 2, then real and imaginary parts for gamma
  --           are expected to be given in c;
  --   c       in c[0] is the real part of the gamma constant,
  --           in c[1] is the imaginary part of the gamma;
  --   vrblvl  is the verbose level.

  function Sweep_Real_Natural_Parameter
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs a real sweep on a natural-parameter homotopy,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2 for the precision level,
  --           respectively double, double double, or quad double;
  --   vrblvl  is the verbose level.

end Sweep_Interface;
