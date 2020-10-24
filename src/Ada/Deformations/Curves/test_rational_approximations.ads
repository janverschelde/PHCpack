with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Test_Rational_Approximations is

-- DESCRIPTION :
--    Tests the results of reorganized thread safe construction of a 
--    vector of Pade approximants, comparing with the previous code.

  procedure Standard_Test ( numdeg,dendeg,dim : in integer32 );
  procedure DoblDobl_Test ( numdeg,dendeg,dim : in integer32 );
  procedure TripDobl_Test ( numdeg,dendeg,dim : in integer32 );
  procedure QuadDobl_Test ( numdeg,dendeg,dim : in integer32 );
  procedure PentDobl_Test ( numdeg,dendeg,dim : in integer32 );
  procedure OctoDobl_Test ( numdeg,dendeg,dim : in integer32 );
  procedure DecaDobl_Test ( numdeg,dendeg,dim : in integer32 );

  -- DESCRIPTION :
  --   Runs a test on a random coefficient vector of series,
  --   in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

  procedure Standard_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 );
  procedure DoblDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 );
  procedure TripDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 );
  procedure QuadDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 );
  procedure PentDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 );
  procedure OctoDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 );
  procedure DecaDobl_Test_Vector ( nbr,numdeg,dendeg,dim : in integer32 );

  -- DESCRIPTION :
  --   Runs a test on a vector of range 1..nbr of random coefficient 
  --   vectors of series, in double, double double, triple double,
  --   quad double, penta double, octo double, or deca double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension of the problem and the precision;
  --   then launches the test.

end Test_Rational_Approximations;
