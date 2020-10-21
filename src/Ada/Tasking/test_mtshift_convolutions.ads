with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Test_mtShift_Convolutions is

-- DESCRIPTION :
--   Tests the multitasked code to shift coefficients of circuits.

  procedure Standard_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to shift the coefficients in double precision.
  --   Compares with the one task run if output.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

  procedure DoblDobl_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to shift the coefficients in double double precision.
  --   Compares with the one task run if output.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

  procedure QuadDobl_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to shift the coefficients in double double precision.
  --   Compares with the one task run if output.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

  procedure Standard_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in Standard_Speelpenning_Convolutions.Circuits;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs a benchmark in double precision.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

  procedure DoblDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs a benchmark in double double precision.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

  procedure QuadDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs a benchmark in quad double precision.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

  procedure Benchmark ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates random circuits in quad double precision,
  --   and runs benchmark tests in all three precisions.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the parameters of the test
  --   and then launches the test.

end Test_mtShift_Convolutions;
