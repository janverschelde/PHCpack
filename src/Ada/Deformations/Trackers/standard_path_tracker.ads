with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Continuation_Data;         use Standard_Continuation_Data;

package Standard_Path_Tracker is

-- DESCRIPTION :
--   This package implements a path tracker with a next() method.

-- CONSTRUCTORS :

  procedure Init ( s : in Link_to_Solution );

  -- DESRIPTION :
  --   Stores s as the current solution, leaves the homotopy as it is.
  --   This is useful for tracking the next path in the same homotopy.

  procedure Init ( p,q : in Link_to_Poly_sys; fixed_gamma : in boolean );
  procedure Init ( p,q : in Link_to_Poly_sys; fixed_gamma : in boolean;
                   s : in Link_to_Solution );

  -- DESCRIPTION :
  --   Initializes the homotopy with target system p, start system q,
  --   and stores the start solution s (if given).
  --   If fixed_gamma is true, then the gamma will be a fixed default value,
  --   otherwise, a new random gamma will be generated. 
  --   The k constants in the homotopy is set to a good value as well.
  --   The condition of the continuation parameters is set to zero.

  procedure Init ( p,q : in Link_to_Poly_sys;
                   gamma : in Complex_Number; k : in natural32 );
  procedure Init ( p,q : in Link_to_Poly_sys; s : in Link_to_Solution;
                   gamma : in Complex_Number; k : in natural32 );

  -- DESCRIPTION :
  --   Initializes the homotopy with target system p, start system q,
  --   and stores the start solution s (if given).
  --   The given gamma and k will be used as constants in the homotopy.
  --   The condition of the continuation parameters is set to zero.

  procedure Init ( p,q : in Link_to_Poly_sys;
                   gamma : in Complex_Number; k,cp : in natural32 );
  procedure Init ( p,q : in Link_to_Poly_sys; s : in Link_to_Solution;
                   gamma : in Complex_Number; k,cp : in natural32 );

  -- DESCRIPTION :
  --   Initializes the homotopy with target system p, start system q,
  --   and stores the start solution s (if given).
  --   The given gamma and k will be used as constants in the homotopy.
  --   The condition of the continuation parameters is set to cp.

  procedure Init ( h : in Link_to_Poly_Sys; txk : in integer32 );
  procedure Init ( h : in Link_to_Poly_Sys; txk : in integer32;
                   s : in Link_to_Solution );

  -- DESCRIPTION :
  --   Initializes the homotopy with the natural parameter homotopy h.
  --   The index txk points to x(k) = t, the variable that will be used
  --   as the continuation parameter.

  -- REQUIRED :
  --   For some n, h'range = 1..n and Number_of_Unknowns(h(i)) = n+1,
  --   for i in h'range.

-- SELECTORS :

  function get_current return Link_to_Solution;
  function get_current return Solu_Info;

  -- DESCRIPTION :
  --   Returns the current solution.

  function get_next return Link_to_Solution;

  -- DESCRIPTION :
  --   Does one step of a predictor-corrector method.
  --   The default for the target of the continuation parameter t is one.

  function get_next ( target_t : Complex_Number ) return Link_to_Solution;

  -- DESCRIPTION :
  --   Does one step of a predictor-corrector method.
  --   for the given target of the continuation parameter.

-- DESCTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the homotopy data and resets all data.

end Standard_Path_Tracker;
