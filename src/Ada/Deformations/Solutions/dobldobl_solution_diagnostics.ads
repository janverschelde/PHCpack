with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Solutions;         use DoblDobl_Complex_Solutions;

package DoblDobl_Solution_Diagnostics is

-- DESCRIPTION :
--   This functions in this package provide diagnostics for solutions
--   of polynomial systems, using double double precision arithmetic.

  function Is_Real ( sol : Solution; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the imaginary part of every component of the 
  --   solution vector of sol is less than or equal to tol in absolute value.

  function Equal ( s1,s2 : Solution; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the difference in absolute value between every
  --   component of the solution vectors s1 and s2 is less than tol.

  function Is_Clustered ( sol : Solution; nb : natural32; sols : Solution_List;
                          tol : double_float ) return natural32;
  function Is_Clustered ( sol : Solution; nb : natural32; sols : Solution_Array;
                          tol : double_float ) return natural32;

  -- DESCRIPTION :
  --   Returns the index of the first other solution in sols,
  --   equal to the solution sol (up to the tolerance tol).
  --   If sol does not occur anywhere else than at position nb in sols,
  --   then that position nb is returned.

  -- ON ENTRY :
  --   sol      a solution occurring in the list sols at position nb;
  --   nb       position of the solution sols in the list sols;
  --   sols     list or arra of solutions;
  --   tol      tolerance to decide whether two solution vectors are equal.

  -- ON RETURN :
  --   Either the index nb if the solution sol is not clustered,
  --   or the index of the first other occurrence of the solution.

  function Multiplicity ( sol : Solution; sols : Solution_List; 
                          tol : double_float ) return natural32;
  function Multiplicity ( sol : Solution; sols : Solution_Array;
                          tol : double_float ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of occurrences of the solution sol
  --   in the list or array sols, using the tolerance tol.

  function At_Infinity ( sol : Solution; prj : boolean;
                         tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Decides whether the given solution lies at infinity.

  -- ON ENTRY :
  --   sol      solution in affine or projective coordinates;
  --   prj      true if sol is in projective coordinates,
  --            false if sol is in affine coordinates;
  --   tol      if prj then sol is at infinity if one of its
  --            coordinates has modulus larger than 1/tol,
  --            otherwise, sol is at infinity if one of its
  --            coordinates has modulas larger than tol.

  -- NOTICE :
  --   tol is typically a large number, e.g.: 1.0E+8,
  --   in contrast to the tolerances in other diagnostics.

end DoblDobl_Solution_Diagnostics;
