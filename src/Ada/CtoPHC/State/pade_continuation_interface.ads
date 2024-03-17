with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

package Pade_Continuation_Interface is

-- DESCRIPTION :
--   Functions in this package define the interface to the new path
--   trackers with Pade approximants as predictors.

  function Pade_Continuation_Parameters_Set_Defaults
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets default values for the homotopy continuation parameters.
  --   The vrblvl is the verbose level.

  function Pade_Continuation_Parameters_Get_Value
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns in b or c one value of a homotopy continuation parameter.

  -- ON ENTRY :
  --   a       a natural number between 1 and 12;
  --   vrblvl  is the verbose level.

  -- RETURN :
  --   b       b[0] depends on the value of a[0]
  --           if a[0] = 2, then b[0] is the degree of Pade numerator,
  --           if a[0] = 3, then b[0] is the degree of Pade denominator,
  --           if a[0] = 12, then b[0] is the maximum #corrector steps,
  --           if a[0] = 13, then b[0] is the maximum #steps on a path.
  --   c       c[0] depends on the value of a[0]
  --           if a[0] = 1, then c[0] contains the real part of gamma,
  --           and c[1] contains the imaginary part of gamma;
  --           if a[0] = 4, then c[0] is the maximum step size,
  --           if a[0] = 5, then c[0] is the minimum step size,
  --           if a[0] = 6, then c[0] is the series step factor,
  --           if a[0] = 7, then c[0] is the pole radius step factor,
  --           if a[0] = 8, then c[0] is the curvature step factor,
  --           if a[0] = 9, then c[0] is the predictor residual tolerance,
  --           if a[0] = 10, then c[0] is the corrector residual tolerance,
  --           if a[0] = 11, then c[0] is the zero coefficient tolerance,

  function Pade_Continuation_Parameters_Set_Value
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets one value of the parameter.

  -- ON ENTRY :
  --   a       a natural number between 1 and 12;
  --   b       b[0] depends on the value of a[0]
  --           if a[0] = 2, then b[0] is the degree of Pade numerator,
  --           if a[0] = 3, then b[0] is the degree of Pade denominator,
  --           if a[0] = 12, then b[0] is the maximum #corrector steps,
  --           if a[0] = 13, then b[0] is the maximum #steps on a path;
  --   c       c[0] depends on the value of a[0]
  --           if a[0] = 1, then c[0] must contain the real part
  --           of gamma and c[1] the imaginary part of gamma,
  --           if a[0] = 4, then c[0] is the maximum step size,
  --           if a[0] = 5, then c[0] is the minimum step size,
  --           if a[0] = 6, then c[0] is the series step factor,
  --           if a[0] = 7, then c[0] is the pole radius step factor,
  --           if a[0] = 8, then c[0] is the curvature step factor,
  --           if a[0] = 9, then c[0] is the predictor residual tolerance,
  --           if a[0] = 10, then c[0] is the corrector residual tolerance,
  --           if a[0] = 11, then c[0] is the zero coefficient tolerance;
  --   vrblvl  is the verbose level.

  function Pade_Continuation_Parameters_Write
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Writes the homotopy continuation parameters to the defined
  --   output file.  The verbose level is in vrblvl.

  function Pade_Continuation_Parameters_Clear
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data to store the values of the parameters.
  --   The vrblvl is the verbose level.

  function Pade_Continuation_Parameters_Reset_Values
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Resets the values of the homotopy continuation parameters
  --   for the step-by-step path trackers, for the precision
  --   double if a[0] = 0, double double if a[0] = 1,
  --   and quad double if a[0] = 2.

  function Pade_Continuation_Track_Paths
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks path for a defined artificial-parameter homotopy,
  --   starting at the solutions in the container.

  -- ON ENTRY:
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --           in a[1] is the number of characters of the name of the
  --           output file, if a[1] is zero, then no output is written,
  --           otherwise, the characters of the output file name are in b,
  --           in a[2] is the value of the verbose flag;
  --           in a[3] is the value of the homogeneous coordinates flag,
  --           in a[4] is the flag to indicate whether the file is local;
  --   b       output file name if a[1] is not zero.

  function Pade_Continuation_Artificial_Homotopy
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes an artificial-parameter homotopy for 
  --   a step-by-step execution.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   b       in b[0] is the verbose flag,
  --           in b[1] is the flag for homogeneous coordinates
  --           (1 if on, 0 for affine);
  --   vrblvl  is the verbose level.

  function Pade_Continuation_Natural_Homotopy
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Initializes a natural-parameter homotopy for 
  --   a step-by-step execution.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   b       in b[0] is the verbose flag,
  --           in b[1] is the flag for homogeneous coordinates
  --           (1 if on, 0 for affine);
  --   vrblvl  is the verbose level.

  function Pade_Continuation_Set_Start_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Sets next start solution in the step-by-step tracker.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --           in a[1] is the index of the solution;
  --   b       in b[0] is the verbose flag;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the failure code of the retrieval.

  function Pade_Continuation_Next_Step
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Runs the next predictor-corrector step in the step-by-step tracker.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --           in a[1] is the index of the solution;
  --   b       in b[0] is the verbose flag;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the failure code of the retrieval.

  function Pade_Continuation_Set_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Places the current solution at the given index in the container.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --           in a[1] is the index of the solution;
  --   b       in b[0] is the verbose flag;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the failure code of the placement.

  function Pade_Continuation_Set_Predicted_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Places the predicted solution at the given index in the container.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --           in a[1] is the index of the solution;
  --   b       in b[0] is the verbose flag;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is the failure code of the placement.

  function Pade_Continuation_Pole_Radius
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Gets the current pole radius in the step-by-step tracker.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   vrblvl  is the verbose level.
  
  -- ON RETURN :
  --   c       in c[0] is the value of the pole radius.

  function Pade_Continuation_Closest_Pole
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Gets the closest pole in the step-by-step tracker.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the real part of the closest pole;
  --           in c[1] is the imaginary part of the closest pole.

  function Pade_Continuation_Get_Pole
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the real and imaginary parts of the pole with given index.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --           in a[1] is the index of the Pade approximant,
  --           in a[2] is the index of the pole;
  --   b       in b[0] is the verbose flag;
  --   vrblvl  is the verbose level.


  -- ON RETURN :
  --   c       in c[0] is the real part of the pole;
  --           in c[1] is the imaginary part of the pole.

  function Pade_Continuation_T_Value
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the current value for t.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the current value for t.

  function Pade_Continuation_Step_Size
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the current step size.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the current step size.

  function Pade_Continuation_Series_Step
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the current series step.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the value of the series step.

  function Pade_Continuation_Pole_Step
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the current pole step.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the value of the pole step.

  function Pade_Continuation_Estimated_Distance
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the estimated distance to the nearest solution.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the estimated distance.

  function Pade_Continuation_Hessian_Step
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the Hessian step.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the value of the Hessian step.

  function Pade_Continuation_Series_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a series coefficient for a given power.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --           in a[1] is the component of the series,
  --           in a[2] is the power for the coefficient;
  --   b       in b[0] is the verbose flag;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the real part of the coefficient,
  --           in c[1] is the imaginary part of the coefficient.

  function Pade_Continuation_Pade_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns a coefficient of a Pade approximant.

  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --           if a[1] = 0, then the denominator coefficient is returned,
  --           otherwise, on return is the nominator coefficient,
  --           in a[2] is the index of the approximant,
  --           in a[3] is the power for the coefficient;
  --   b       in b[0] is the verbose flag;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   c       in c[0] is the real part of the coefficient,
  --           in c[1] is the imaginary part of the coefficient.

  function Pade_Continuation_Clear_Data
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Clears the data used in the step-by-step trackers.
 
  -- ON ENTRY :
  --   a       in a[0] is 0, 1, or 2, respectively for double, double double,
  --           or quad double precision;
  --   vrblvl  is the verbose level.

end Pade_Continuation_Interface;
