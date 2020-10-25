with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_Matrices;
with PentDobl_Complex_VecMats;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_Matrices;
with OctoDobl_Complex_VecMats;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_Matrices;
with DecaDobl_Complex_VecMats;
with DecaDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Solutions;
with Standard_Complex_Circuits;
with Standard_Coefficient_Circuits;
with DoblDobl_Complex_Circuits;
with TripDobl_Complex_Circuits;
with QuadDobl_Complex_Circuits;
with PentDobl_Complex_Circuits;
with OctoDobl_Complex_Circuits;
with DecaDobl_Complex_Circuits;

package Test_mtHessian_Circuits is

-- DESCRIPTION :
--   Tests the Hessian criterion, for systems of complex circuits,
--   in double, double double, triple double, quad double,
--   penta double, octo double, and deca double arithmetic,
--   with multitasking for shared memory parallel computers.

  procedure Write_Singular_Values 
              ( values : in Standard_Complex_VecVecs.VecVec );
  procedure Write_Singular_Values 
              ( values : in DoblDobl_Complex_VecVecs.VecVec );
  procedure Write_Singular_Values 
              ( values : in TripDobl_Complex_VecVecs.VecVec );
  procedure Write_Singular_Values 
              ( values : in QuadDobl_Complex_VecVecs.VecVec );
  procedure Write_Singular_Values 
              ( values : in PentDobl_Complex_VecVecs.VecVec );
  procedure Write_Singular_Values 
              ( values : in OctoDobl_Complex_VecVecs.VecVec );
  procedure Write_Singular_Values 
              ( values : in DecaDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.

  procedure Standard_Test
              ( s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                static,output : in boolean );

  -- DESCRIPTION :
  --   Runs the test in double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

  procedure Standard_Coefficient_Test
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                output : in boolean );

  -- DESCRIPTION :
  --   Runs the test in double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

  procedure DoblDobl_Test
              ( s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean );

  -- DESCRIPTION :
  --   Runs the test in double double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

  procedure TripDobl_Test
              ( s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean );

  -- DESCRIPTION :
  --   Runs the test in triple double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

  procedure QuadDobl_Test
              ( s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean );

  -- DESCRIPTION :
  --   Runs the test in quad double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

  procedure PentDobl_Test
              ( s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean );

  -- DESCRIPTION :
  --   Runs the test in penta double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

  procedure OctoDobl_Test
              ( s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean );

  -- DESCRIPTION :
  --   Runs the test in octo double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

  procedure DecaDobl_Test
              ( s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean );

  -- DESCRIPTION :
  --   Runs the test in deca double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

  procedure Standard_Random_Test
              ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double precision
  --   and then launches the test.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure DoblDobl_Random_Test
              ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double double precision,
  --   and then launches the test.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure TripDobl_Random_Test
              ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in triple double precision,
  --   and then launches the test.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure QuadDobl_Random_Test
              ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure PentDobl_Random_Test
              ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure OctoDobl_Random_Test
              ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in octo double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure DecaDobl_Random_Test
              ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in deca double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure Standard_User_Test;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions
  --   and launches the test in double precision.

  procedure DoblDobl_User_Test;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions
  --   and launches the test in double double precision.

  procedure TripDobl_User_Test;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions
  --   and launches the test in triple double precision.

  procedure QuadDobl_User_Test;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions
  --   and launches the test in quad double precision.

  procedure PentDobl_User_Test;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions
  --   and launches the test in penta double precision.

  procedure OctoDobl_User_Test;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions
  --   and launches the test in octo double precision.

  procedure DecaDobl_User_Test;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions
  --   and launches the test in deca double precision.

  procedure Standard_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment on the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output is needed.

  procedure Standard_Coefficient_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in double precision,
  --   on coefficient circuits.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment on the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output is needed.

  procedure DoblDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in double double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment on the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of numbers for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output is needed.

  procedure TripDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in triple double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment on the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of numbers for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output is needed.

  procedure QuadDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in quad double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment to the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output needs to be written to screen.

  procedure PentDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in quad double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment to the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output needs to be written to screen.

  procedure OctoDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in octo double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment to the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output needs to be written to screen.

  procedure DecaDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in deca double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment to the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output needs to be written to screen.

  function Prompt_for_Sequence
             ( max : in integer32 )
             return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Prompts the user for a sequence of numbers
  --   and the sequence is returned in a vector.
  --   The length of the sequence may not exceed max.

  procedure Benchmark ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random system of circuits in quad double precision,
  --   and runs benchmark tests in all three precisions.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure Benchmark
              ( p : in DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in DecaDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   For the given polynomial system p with some solutions in sols,
  --   prompts for the parameters and then runs benchmarks.

  function Prompt_for_Precision return character;

  -- DESCRIPTION :
  --   Prompts the user for the work precision
  --   and returns '1', '2', '3', '4', '5', '6', or '7', respectively
  --   for double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Launches the test after prompting for the parameters
  --   to generate a random problem.

end Test_mtHessian_Circuits;
