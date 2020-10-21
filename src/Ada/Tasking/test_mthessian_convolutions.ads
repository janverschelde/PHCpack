with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Test_mtHessian_Convolutions is

-- DESCRIPTION :
--   Tests the development of the Hessian criterion,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

  procedure Write_Singular_Values 
              ( values : in Standard_Complex_VecVecs.VecVec;
                jmvals : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.
  --   The vector jmvals contains the singular values of the Jacobian.

  procedure Write_Singular_Values 
              ( values : in DoblDobl_Complex_VecVecs.VecVec;
                jmvals : in DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.
  --   The vector jmvals contains the singular values of the Jacobian.

  procedure Write_Singular_Values 
              ( values : in QuadDobl_Complex_VecVecs.VecVec;
                jmvals : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.
  --   The vector jmvals contains the singular values of the Jacobian.

  procedure Standard_Test
              ( dim : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.Link_to_VecVec );
  procedure DoblDobl_Test
              ( dim : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.Link_to_VecVec );
  procedure QuadDobl_Test
              ( dim : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Runs the test in double, double double, or quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution.

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double precision
  --   and then launches the test.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double double precision,
  --   and then launches the test.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
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

  procedure QuadDobl_User_Test;

  -- DESCRIPTION :
  --   Prompts for a polynomial system with solutions
  --   and launches the test in quad double precision.

  procedure Standard_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_Vectors.Vector;
                verbose : in boolean := false );
  procedure DoblDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Vector;
                verbose : in boolean := false );
  procedure QuadDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs a benchmark test in double, double double,
  --   or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output is needed.

  procedure Benchmark ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision,
  --   and runs benchmark tests in all three precisions.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  function Prompt_for_Precision return character;

  -- DESCRIPTION :
  --   Prompts the user for the work precision
  --   and returns '0', '1', or '2',
  --   for double, double double, or quad double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Launches the test after prompting for the parameters
  --   to generate a random problem.

end Test_mtHessian_Convolutions;
