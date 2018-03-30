with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Dobldobl_Complex_Poly_Systems;
with Dobldobl_Complex_Laur_Systems;
with Quaddobl_Complex_Laur_Systems;
with Quaddobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Multitasking_Root_Refiners is

-- DESCRIPTION :
--   Newton's method is executed with multitasking on a list of solutions
--   to verify the quality of the solution list.
--   The procedures come in different flavors:
--   (1) whether the systems is Laurent or simply polynomial;
--   (2) whether the refiner stays mute or whether the multitasking 
--       is silent or reporting its progress,
--   (3) whether the precision is double, double double, or quad double;
--   so therefore this package offers 2 x 3 x 3 = 18 procedures.

  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );

  -- DESCRIPTION :
  --   Applies the root refinement to the solutions in sols
  --   and the polynomial system in p with nt tasks.
  --   Mute root refiners do not write a summary table at the end.

  -- ON ENTRY :
  --   nt       number of tasks;
  --   p        a square polynomial system;
  --   sols     initial solutions to p;
  --   epsxa    maximum absolute error on the zero;
  --   epsfa    maximum absolute value for the residue;
  --   tolsing  tolerance on inverse condition number for singular solution;
  --   numit    the number of iterations, to be initialized with zero;
  --   max      maximum number of iterations per zero;
  --   deflate  if true, apply deflation to singular solutions.

  -- ON RETURN :
  --   sols     refined solutions to p.

  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );

  -- DESCRIPTION :
  --   Applies the root refinement to the solutions in sols
  --   and the polynomial system in p with nt tasks.
  --   No intermediate output is written to screen during the computations.

  -- ON ENTRY :
  --   file     output file for the results;
  --   nt       number of tasks;
  --   p        a square polynomial system;
  --   sols     initial solutions to p;
  --   epsxa    maximum absolute error on the zero;
  --   epsfa    maximum absolute value for the residue;
  --   tolsing  tolerance on inverse condition number for singular solution;
  --   numit    the number of iterations, to be initialized with zero;
  --   max      maximum number of iterations per zero;
  --   deflate  if true, apply deflation to singular solutions.

  -- ON RETURN :
  --   sols     refined solutions to p.

  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );

  -- DESCRIPTION :
  --   Applies the root refinement to the solutions in sols
  --   and the polynomial system in p with nt tasks.
  --   Intermediate output is written to screen during the computations.

  -- ON ENTRY :
  --   file     output file for the results;
  --   nt       number of tasks;
  --   p        a square polynomial system;
  --   sols     initial solutions to p;
  --   epsxa    maximum absolute error on the zero;
  --   epsfa    maximum absolute value for the residue;
  --   tolsing  tolerance on inverse condition number for singular solution;
  --   numit    the number of iterations, to be initialized with zero;
  --   max      maximum number of iterations per zero;
  --   deflate  if true, apply deflation to singular solutions.

  -- ON RETURN :
  --   sols     refined solutions to p.

  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );

  -- DESCRIPTION :
  --   Applies the root refinement to the solutions in sols
  --   and the Laurent polynomial system in p with nt tasks.
  --   Mute root refiners write no summary table at the end.

  -- ON ENTRY :
  --   p        a square Laurent polynomial system;
  --   sols     initial solutions to p;
  --   epsxa    maximum absolute error on the zero;
  --   epsfa    maximum absolute value for the residue;
  --   tolsing  tolerance on inverse condition number for singular solution;
  --   numit    the number of iterations, to be initialized with zero;
  --   max      maximum number of iterations per zero;
  --   deflate  if true, apply deflation to singular solutions.

  -- ON RETURN :
  --   sols     refined solutions to p;
  --   numit    the number of iterations;
  --   deflate  if set to false, then the system was too large to deflate.

  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );

  -- DESCRIPTION :
  --   Applies the root refinement to the solutions in sols
  --   and the Laurent polynomial system in p with nt tasks.
  --   No intermediate output is written to screen during the computations.

  -- ON ENTRY :
  --   file     output file for the results;
  --   nt       number of tasks;
  --   p        a square Laurent polynomial system;
  --   sols     initial solutions to p;
  --   epsxa    maximum absolute error on the zero;
  --   epsfa    maximum absolute value for the residue;
  --   tolsing  tolerance on inverse condition number for singular solution;
  --   numit    the number of iterations, to be initialized with zero;
  --   max      maximum number of iterations per zero;
  --   deflate  if true, apply deflation to singular solutions.

  -- ON RETURN :
  --   sols     refined solutions to p;
  --   numit    the number of iterations;
  --   deflate  if set to false, then the system was too large to deflate.

  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );
  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean );

  -- DESCRIPTION :
  --   Applies the root refinement to the solutions in sols
  --   and the Laurent polynomial system in p with nt tasks.
  --   Intermediate output is written to screen during the computations.

  -- ON ENTRY :
  --   file     output file for the results;
  --   nt       number of tasks;
  --   p        a square Laurent polynomial system;
  --   sols     initial solutions to p;
  --   epsxa    maximum absolute error on the zero;
  --   epsfa    maximum absolute value for the residue;
  --   tolsing  tolerance on inverse condition number for singular solution;
  --   numit    the number of iterations, to be initialized with zero;
  --   max      maximum number of iterations per zero;
  --   deflate  if true, apply deflation to singular solutions.

  -- ON RETURN :
  --   sols     refined solutions to p;
  --   numit    the number of iterations;
  --   deflate  if set to false, then the system was too large to deflate.

end Multitasking_Root_Refiners;
