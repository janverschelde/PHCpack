with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package BKK_Bound_Computations is

-- DESCRIPTION :
--   This package exports some routines for computing the BKK bound
--   and solving a random coefficient system by polyhedral continuation.
--   These function are black box routines: the user does not have to
--   worry about intermediate data structures.

  function BKK_by_Implicit_Lifting ( p : Poly_Sys ) return natural32;
  function BKK_by_Implicit_Lifting ( file : file_type; p : Poly_Sys )
                                   return natural32;

  function BKK_by_Static_Lifting ( p : Poly_Sys ) return natural32;
  function BKK_by_Static_Lifting ( file : file_type; p : Poly_Sys )
                                 return natural32;

  -- DESCRIPTION :
  --   If a file is specified, then the mixed subdivision will be written
  --   on that file.  Either implicit or random static lifting can be used.

  function Solve_by_Implicit_Lifting ( p : Poly_Sys ) return Solution_List;
  function Solve_by_Implicit_Lifting ( file : file_type; p : Poly_Sys ) 
                                     return Solution_List;

  function Solve_by_Static_Lifting ( p : Poly_Sys ) return Solution_List;
  function Solve_by_Static_Lifting ( file : file_type; p : Poly_Sys ) 
                                   return Solution_List;

  -- DESCRIPTION :
  --   If a file is specified, then intermediate results will be written
  --   on that file.  Either implicit or random static lifting can be used.

end BKK_Bound_Computations;
