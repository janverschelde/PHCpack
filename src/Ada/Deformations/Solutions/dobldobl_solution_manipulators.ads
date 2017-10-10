with DoblDobl_Complex_Numbers;            use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Solutions;          use DoblDobl_Complex_Solutions;

package DoblDobl_Solution_Manipulators is

-- DESCRIPTION :
--   A solution manipulator operates on a solution and list of solutions.
--   The multitasked path trackers use the imaginary part of the target
--   for each solution to store the length of the path.
--   This imaginary part causes problems with the solution splitters.
--   The procedures in this package set the imaginary part of the target
--   to zero in each solution.

  procedure Remove_Imaginary_Part ( t : in out Complex_Number );

  -- DESCRIPTION :
  --   Removes the imaginary part of t, setting it to zero.

  procedure Remove_Imaginary_Target ( s : in out Link_to_Solution );

  -- DESCRIPTION :
  --   Removes the imaginary part of s.t, setting it to zero.

  procedure Remove_Imaginary_Target ( s : in out Solution_List );

  -- DESCRIPTION :
  --   Removes the imaginary part of every solution in s,
  --   setting it to zero.

end DoblDobl_Solution_Manipulators;
