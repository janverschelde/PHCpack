with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Coefficient_Convolutions;
with DoblDobl_Coefficient_Convolutions;
with QuadDobl_Coefficient_Convolutions;

package Test_mtAlgoDiff_Convolutions is

-- DESCRIPTION :
--   Tests algorithmic differentiation to evaluate and differentiate
--   polynomial systems at truncated power series with multitasking.

  procedure Standard_Test
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to evaluate and differentiate c,
  --   in double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   c        a sequence of convolution circuits;
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

  procedure Standard_Coefficient_Test
              ( c : in Standard_Coefficient_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                static : in boolean; output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to evaluate and differentiate c,
  --   in double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   c        a sequence of coefficient convolution circuits;
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbt      the number of tasks;
  --   static   flag to apply static load balancing;
  --   output   flag for output during multitasking.

  procedure DoblDobl_Test
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to evaluate and differentiate c,
  --   in double double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   c        a sequence of convolution circuits;
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

  procedure DoblDobl_Coefficient_Test
              ( c : in DoblDobl_Coefficient_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                static : in boolean; output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to evaluate and differentiate c,
  --   in double double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   c        a sequence of coefficient convolution circuits;
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbt      the number of tasks;
  --   static   flag to apply static load balancing;
  --   output   flag for output during multitasking.

  procedure QuadDobl_Test
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to evaluate and differentiate c,
  --   in double double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   c        a sequence of convolution circuits;
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

  procedure QuadDobl_Coefficient_Test
              ( c : in QuadDobl_Coefficient_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                static : in boolean; output : in boolean := true );

  -- DESCRIPTION :
  --   Applies multitasking to evaluate and differentiate c,
  --   in quad double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   c        a sequence of coefficient convolution circuits;
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbt      the number of tasks;
  --   static   flag to apply static load balancing;
  --   output   flag for output during multitasking.

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking,
  --   in double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking,
  --   in double double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true );

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking,
  --   in double double precision.  Compares with the one task run.

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
                x : in Standard_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs a benchmark in double precision.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   x        where to evaluate c;
  --   mxe      maximum exponent indices;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

  procedure Standard_Coefficient_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in Standard_Coefficient_Convolutions.Circuits;
                x : in Standard_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs a benchmark in double precision
  --   on coefficient convolution circuits.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   x        where to evaluate c;
  --   mxe      maximum exponent indices;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

  procedure DoblDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs a benchmark in double double precision.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   x        where to evaluate c;
  --   mxe      maximum exponent indices;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

  procedure DoblDobl_Coefficient_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in DoblDobl_Coefficient_Convolutions.Circuits;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs a benchmark in double double precision
  --   on coefficient convolution circuits.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   x        where to evaluate c;
  --   mxe      maximum exponent indices;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

  procedure QuadDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs a benchmark in quad double precision.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   x        where to evaluate c;
  --   mxe      maximum exponent indices;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

  procedure QuadDobl_Coefficient_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in QuadDobl_Coefficient_Convolutions.Circuits;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs a benchmark in quad double precision
  --   on coefficient convolution circuits.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   x        where to evaluate c;
  --   mxe      maximum exponent indices;
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

  function Prompt_for_Precision return character;

  -- DESCRIPTION :
  --   Prompts the user for the precision and returns '0', '1', or '2'
  --   respectively for double, double double, or quad double precision.

  procedure Standard_Test_Problem;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and launches the test in double precision.

  procedure DoblDobl_Test_Problem;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and launches the test in double double precision.

  procedure QuadDobl_Test_Problem;

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and launches the test in quad double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the number of tasks, for the dimensions of the problem,
  --   generates then the data for the problem, and runs tests.

end Test_mtAlgoDiff_Convolutions;
