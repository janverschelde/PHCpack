with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;

package Newton_Power_Convolutions is

-- DESCRIPTION :
--   Several steps with Newton's method on power series are performed
--   with convolutions on the coefficients and linearization.

-- NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                nbrit : in integer32; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type; 
                csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                nbrit : in integer32; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf,
  --   in double, double double, or quad double precision,
  --   using LU factorization to solve the linear systems.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   nbrit    maximum number of iterations;
  --   wrk      work space for the matrix series solver;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   info     info from the LU factorization;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

end Newton_Power_Convolutions;
