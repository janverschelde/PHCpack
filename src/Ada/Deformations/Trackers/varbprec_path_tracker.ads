with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;

package Varbprec_Path_Tracker is

-- DESCRIPTION :
--   This package implements a path tracker with a next() method
--   in variable precision.

-- CONSTRUCTORS :

  procedure Init ( s : in Link_to_String; nv : in natural32 );

  -- DESRIPTION :
  --   Stores s as the current solution, leaves the homotopy as it is.
  --   This is useful for tracking the next path in the same homotopy.
  --   The number of variables in the solution must equal nv,
  --   as needed to initialize the standard path tracker.

  procedure Init ( p,q : in Link_to_Array_of_Strings;
                   fixed_gamma : in boolean );
  procedure Init ( p,q : in Link_to_Array_of_Strings;
                   fixed_gamma : in boolean; s : in Link_to_String );

  -- DESCRIPTION :
  --   Initializes the homotopy with target system p, start system q,
  --   and stores the start solution s (if given).
  --   If fixed_gamma is true, then the gamma will be a fixed default value,
  --   otherwise, a new random gamma will be generated. 
  --   The k constants in the homotopy is set to a good value as well.
  --   The condition of the continuation parameters is set to zero.

  procedure Init ( p,q : in Link_to_Array_of_Strings; k : in natural32;
                   gamma : in Standard_Complex_Numbers.Complex_Number );
  procedure Init ( p,q : in Link_to_Array_of_Strings; k : in natural32;
                   gamma : in Standard_Complex_Numbers.Complex_Number; 
                   s : in Link_to_String );

  -- DESCRIPTION :
  --   Initializes the homotopy with target system p, start system q,
  --   and stores the start solution s (if given).
  --   The given gamma and k will be used as constants in the homotopy.
  --   The condition of the continuation parameters is set to zero.

  procedure Init ( p,q : in Link_to_Array_of_Strings; k,cp : in natural32;
                   gamma : in Standard_Complex_Numbers.Complex_Number );
  procedure Init ( p,q : in Link_to_Array_of_Strings; k,cp : in natural32;
                   gamma : in Standard_Complex_Numbers.Complex_Number;
                   s : in Link_to_String );

  -- DESCRIPTION :
  --   Initializes the homotopy with target system p, start system q,
  --   and stores the start solution s (if given).
  --   The given gamma and k will be used as constants in the homotopy.
  --   The condition of the continuation parameters is set to cp.

  procedure Init ( h : in Link_to_Array_of_Strings; txk : in integer32 );
  procedure Init ( h : in Link_to_Array_of_Strings; txk : in integer32;
                   s : in Link_to_String );

  -- DESCRIPTION :
  --   Initializes the homotopy with the natural parameter homotopy h.
  --   The index txk points to x(k) = t, the variable that will be used
  --   as the continuation parameter.

  -- REQUIRED :
  --   For some n, h'range = 1..n and Number_of_Unknowns(h(i)) = n+1,
  --   for i in h'range.

-- SELECTORS :

  function get_current return Link_to_String;

  -- DESCRIPTION :
  --   Returns the current solution.

  function get_next ( want,maxprc,maxitr : natural32; output : boolean )
                    return Link_to_String;

  -- DESCRIPTION :
  --   Does one step of a predictor-corrector method.
  --   The default for the target of the continuation parameter t is one.

  -- ON ENTRY :
  --   want     wanted number of decimal places of accuracy;
  --   maxprc   maximum precision used in variable precision corrector;
  --   maxitr   maximum number of corrector steps.

  -- ON RETURN :
  --   string representation of a solution in PHCpack format,
  --   the triplet (err, rco, res) indicates the success of the correction.

  function get_next ( target_t : Standard_Complex_Numbers.Complex_Number;
                      want,maxprc,maxitr : natural32; output : boolean )
                    return Link_to_String;

  -- DESCRIPTION :
  --   Does one step of a predictor-corrector method.
  --   for the given target of the continuation parameter.

  -- ON ENTRY :
  --   target_t value of the homotopy continuation parameter
  --            at its destination, the get_next will not overshoot
  --            this target value for t; 
  --   want     wanted number of decimal places of accuracy;
  --   maxprc   maximum precision used in variable precision corrector;
  --   maxitr   maximum number of corrector steps.

  -- ON RETURN :
  --   string representation of a solution in PHCpack format,
  --   the triplet (err, rco, res) indicates the success of the correction.

-- DESCTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the homotopy data and resets all data.

end Varbprec_Path_Tracker;
