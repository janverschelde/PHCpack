with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Main_Eqn_by_Eqn_Solvers is

-- DESCRIPTION :
--   This package provides access to the equation-by-equation solvers.
--   Since the order may have a dramatic impact on the efficiency,
--   the user has the opportunity to shuffle the equations in the system.

  function Shuffle ( file : file_type; p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Interactive Menu to decide which order to take.
  --   The permuted system is written to file and returned.

  procedure Shuffle_Polynomials_and_Solve
              ( file : in file_type; filename : in string;
                p : in Poly_Sys );

  -- DESCRIPTION :
  --   Presents the user with a menu to shuffle the polynomials
  --   and calls then the equation-by-equation solver.

  procedure Solver_using_Generics
              ( file : in file_type; filename : in string;
                p : in Poly_Sys );

  -- DESCRIPTION :
  --   Prepares the executable versions to instantiate the
  --   generic equation-by-equation solver.


  procedure Main ( infilename,outfilename : in string );

  -- DESCRIPTION :
  --   Defines phc -a, applies the equation-by-equation solver.

  -- ON ENTRY :
  --   infilename     file name for input for a polynomial system,
  --                  if empty, then the user will be prompted to supply
  --                  the name of an input file;
  --   outfilename    name of file for output, if empty, then the user will
  --                  be asked to supply the name of an output file.

end Main_Eqn_by_Eqn_Solvers;
