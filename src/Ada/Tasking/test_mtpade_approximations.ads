with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;

package Test_mtPade_Approximations is

-- DESCRIPTION :
--   Tests multitasked algorithms for rational approximations.

  procedure Standard_Allocate
              ( nbr,numdeg,dendeg : in integer32;
                cff : out Standard_Complex_VecVecs.VecVec;
                numcff1,numcff2 : out Standard_Complex_VecVecs.VecVec;
                dencff1,dencff2 : out Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Generates the input series coefficients and allocates space
  --   for the numerators and the denominators, in double precision.

  -- ON ENTRY :
  --   nbr       the number of components in the vector of approximants;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators.

  -- ON RETURN :
  --   cff       random coefficients for a series;
  --   numcff1   first space allocated for numerator coefficients;
  --   numcff2   second space allocated for numerator coefficients;
  --   dencff1   first space allocated for denominator coefficients;
  --   dencff2   second space allocated for denominator coefficients.

  procedure DoblDobl_Allocate
              ( nbr,numdeg,dendeg : in integer32;
                cff : out DoblDobl_Complex_VecVecs.VecVec;
                numcff1,numcff2 : out DoblDobl_Complex_VecVecs.VecVec;
                dencff1,dencff2 : out DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Generates the input series coefficients and allocates space
  --   for the numerators and the denominators, in double double precision.

  -- ON ENTRY :
  --   nbr       the number of components in the vector of approximants;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators.

  -- ON RETURN :
  --   cff       random coefficients for a series;
  --   numcff1   first space allocated for numerator coefficients;
  --   numcff2   second space allocated for numerator coefficients;
  --   dencff1   first space allocated for denominator coefficients;
  --   dencff2   second space allocated for denominator coefficients.

  procedure QuadDobl_Allocate
              ( nbr,numdeg,dendeg : in integer32;
                cff : out QuadDobl_Complex_VecVecs.VecVec;
                numcff1,numcff2 : out QuadDobl_Complex_VecVecs.VecVec;
                dencff1,dencff2 : out QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Generates the input series coefficients and allocates space
  --   for the numerators and the denominators, in quad double precision.

  -- ON ENTRY :
  --   nbr       the number of components in the vector of approximants;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators.

  -- ON RETURN :
  --   cff       random coefficients for a series;
  --   numcff1   first space allocated for numerator coefficients;
  --   numcff2   second space allocated for numerator coefficients;
  --   dencff1   first space allocated for denominator coefficients;
  --   dencff2   second space allocated for denominator coefficients.

  procedure Standard_Test ( nbr,numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Makes a vector of range 1..nbr of random coefficients to construct
  --   rational approximations of the given degrees, in double precision.

  -- ON ENTRY :
  --   nbr       the number of components in the vector of approximants;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators.

  procedure DoblDobl_Test ( nbr,numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Makes a vector of range 1..nbr of random coefficients to construct
  --   rational approximations of the given degrees,
  --   in double double precision.

  -- ON ENTRY :
  --   nbr       the number of components in the vector of approximants;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators.

  procedure QuadDobl_Test ( nbr,numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Makes a vector of range 1..nbr of random coefficients to construct
  --   rational approximations of the given degrees,
  --   in double double precision.

  -- ON ENTRY :
  --   nbr       the number of components in the vector of approximants;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators.

  procedure Standard_Benchmark
              ( file : in file_type;
                nbr,numdeg,dendeg,nbruns,inc : in integer32 );

  -- DESCRIPTION :
  --   Runs a benchmark test in double precision.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbr      number of components;
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   nbruns   the number of multitasked runs;
  --   inc      increment on the number of tasks;

  procedure DoblDobl_Benchmark
              ( file : in file_type;
                nbr,numdeg,dendeg,nbruns,inc : in integer32 );

  -- DESCRIPTION :
  --   Runs a benchmark test in double double precision.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbr      number of components;
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   nbruns   the number of multitasked runs;
  --   inc      increment on the number of tasks;

  procedure QuadDobl_Benchmark
              ( file : in file_type;
                nbr,numdeg,dendeg,nbruns,inc : in integer32 );

  -- DESCRIPTION :
  --   Runs a benchmark test in quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbr      number of components;
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   nbruns   the number of multitasked runs;
  --   inc      increment on the number of tasks;

  procedure Benchmark ( nbr,numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Generates random problems
  --   and runs benchmark tests in all three precisions.

  -- ON ENTRY :
  --   nbr      number of components;
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension of the problem
  --   and then launches the test.

end Test_mtPade_Approximations;
