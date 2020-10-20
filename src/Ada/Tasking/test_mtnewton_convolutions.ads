with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Test_mtNewton_Convolutions is

-- DESCRIPTION :
--   Tests the development of Newton's method on power series
--   with the reverse mode of algorithmic differentation
--   and linearization to solve the matrix series equations,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

  procedure Apply_Fabry
              ( c : in Standard_Complex_VecVecs.VecVec;
                verbose : in boolean := true );
  procedure Apply_Fabry
              ( c : in DoblDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true );
  procedure Apply_Fabry
              ( c : in QuadDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Estimates the radius of convergence of the power series
  --   with the application of the theorem of Fabry,
  --   in double, double double, or quad double precision.

  procedure Standard_Run
              ( nbt,dim,maxit : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true );
  procedure DoblDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true );
  procedure QuadDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs Newton's method with nbt tasks,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   nbt        the number of tasks;
  --   dim        number of equations and variables;
  --   maxit      maximum number of iterations;
  --   s          system of convolution circuits;
  --   scf        its leading coefficients are solution vector;
  --   serelp     the previous elapsed wall clock time of a serial run;
  --   mltelp     the previous elapsed wall clock time of a multitasked run;
  --   speedup    the previous speedup of a multitasked run;
  --   output     if true, then Newton is verbose, else silent;
  --   estco      if true, then the condition number is estimated;
  --   verbose    if true, then timings are shown, otherwise not.

  -- ON RETURN :
  --   serelp     updated elapsed wall clock time of a serial run,
  --              if nbt = 1 and the user did not want multitasking;
  --   mltelp     updated elapsed wall clock time of a multitasked run,
  --              if nbt > 1;
  --   speedup    computed speedup if serelp /= 0.0.
  --   efficiency computed if serelp /= 0.0.


  procedure Standard_Run_Loop
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in Standard_Complex_Vectors.Vector;
	        deg : in integer32 );
  procedure DoblDobl_Run_Loop
              ( p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in DoblDobl_Complex_Vectors.Vector;
	        deg : in integer32 );
  procedure QuadDobl_Run_Loop
              ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in QuadDobl_Complex_Vectors.Vector;
	        deg : in integer32 );

  -- DESCRIPTION :
  --   Runs Newton's method on a solution sol of the system p,
  --   with power series of degree deg, 
  --   in double, double double, or quad double precision.

  procedure Standard_Test;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and tests in double precision.

  procedure DoblDobl_Test;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and tests in double double precision.

  procedure QuadDobl_Test;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and tests in quad double precision.

  procedure Standard_Random_Test ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Tests on a random Newton homotopy in double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure DoblDobl_Random_Test ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Tests on a random Newton homotopy in double double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure QuadDobl_Random_Test ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Tests on a random Newton homotopy in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure Standard_Benchmark
              ( file : in file_type; nbruns,inc,maxit : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.Link_to_VecVec;
                verbose : in boolean := false );
  procedure DoblDobl_Benchmark
              ( file : in file_type; nbruns,inc,maxit : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                verbose : in boolean := false );
  procedure QuadDobl_Benchmark
              ( file : in file_type; nbruns,inc,maxit : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in double, double double,
  --   or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be uses;
  --   inc      increment on the number of tasks, if nbruns /= 0;
  --   maxit    the maximum number of iterations;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output is needed.

  function Prompt_for_Sequence
             ( max : in integer32 )
             return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Prompts the user for a sequence of numbers
  --   and the sequence is returned in a vector.
  --   The length of the sequence may not exceed max.

  procedure Benchmark ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision,
  --   and runs benchmark tests in all three precisions.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure Benchmark
              ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim,deg : in integer32 );

  -- DESCRIPTION :
  --   For the given polynomial system p and some solutions in sols,
  --   prompts the user for the parameters of the benchmark runs
  --   in all three levels of precision.
  
  procedure Prompt_for_Dimensions
              ( dim,deg,nbr,pwr : in out integer32 );

  -- DESCRIPTION :
  --   Prompts the user of the dimensions of the random input data.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the setup of the tests and then runs the tests.

end Test_mtNewton_Convolutions;
