with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Solutions;         use Multprec_Complex_Solutions;

package Multprec_Solution_Diagnostics is

  function Is_Real ( sol : Solution; tol : Floating_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the imaginary part of every component of the 
  --   solution vector of sol is less than or equal to tol in absolute value.

  function Equal ( s1,s2 : Solution; tol : Floating_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the difference in absolute value between every
  --   component of the solution vectors s1 and s2 is less than tol.

  function Is_Clustered ( sol : Solution; nb : natural32;
                          sols : Solution_List; tol : Floating_Number )
                        return natural32;
  function Is_Clustered ( sol : Solution; nb : natural32;
                          sols : Solution_Array; tol : Floating_Number )
                        return natural32;

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
                          tol : Floating_Number ) return natural32;
  function Multiplicity ( sol : Solution; sols : Solution_Array;
                          tol : Floating_Number ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of occurrences of the solution sol
  --   in the list or array sols, using the tolerance tol.

end Multprec_Solution_Diagnostics;
