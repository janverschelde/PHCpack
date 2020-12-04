with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Symbol_Table;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;

package Main_Decomposition is

-- DESCRIPTION :
--   Exports the main procedures to compute a numerical irreducible
--   decomposition of the solution set of a polynomial system.

  procedure Read_Two_Witness_Sets
              ( lp1,lp2 : out Link_to_Poly_Sys;
                sols1,sols2 : out Solution_List;
                dim1,dim2 : out natural32; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for two files for two witness sets.

  procedure Read_Two_Witness_Sets
              ( lp1,lp2 : out Link_to_Poly_Sys;
                sols1,sols2 : out Solution_List;
                dim1,dim2 : out natural32;
                lsym1,lsym2 : out Symbol_Table.Link_to_Array_of_Symbols;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for two files for two witness sets,
  --   returns the witness sets: systems and solutions, their dimensions,
  --   and their symbols for the names of the variables.

  procedure Call_Extrinsic_Diagonal_Homotopies
              ( outfilename : in string; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for two witness sets and runs an extrinsic
  --   cascade of homotopies to compute a witness set representation
  --   of their intersection.

  procedure Call_Intrinsic_Diagonal_Homotopies
               ( outfilename : in string; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the diagonal homotopies in intrinsic coordinates.

  procedure Call_Binomial_Solver
               ( infilename,outfilename : in string;
                 vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls the solver for binomial systems.

  procedure Degrees_of_Monomial_Maps
              ( infilename : in string; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Reads a Laurent system and a list of monomial maps
  --   to compute their degrees.

  procedure Transform_Positive_Corank
               ( infilename,outfilename : in string;
                 vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the rank of the supports of a Laurent system
  --   and if the corank is positive, via a monomial transformation
  --   as many variables as the corank can be eliminated.

  procedure Run_Cascade_Filter
              ( infilename,outfilename : in string; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the cascade filter with verbose level in vrb.

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines phc -c to compute a numerical irreducible decomposition.

  -- ON ENTRY :
  --   nt             the number of tasks, if 0 then no multitasking,
  --                  otherwise nt tasks will be used to track the paths;
  --   infilename     file name for input (a polynomial system),
  --                  if empty, then the user will be prompted to supply
  --                  the name of an input file;
  --   outfilename    name of file for output, if empty, then the user will
  --                  be asked to supply the name of an output file;
  --   verbose        the verbose level.

end Main_Decomposition;
