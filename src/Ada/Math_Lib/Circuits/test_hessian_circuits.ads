with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with PentDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Poly_Systems;
with Standard_Complex_Circuits;
with Standard_Coefficient_Circuits;
with DoblDobl_Complex_Circuits;
with TripDobl_Complex_Circuits;
with QuadDobl_Complex_Circuits;
with PentDobl_Complex_Circuits;
with OctoDobl_Complex_Circuits;
with DecaDobl_Complex_Circuits;

package Test_Hessian_Circuits is

-- DESCRIPTION :
--   Test and development of better algorithms to compute Hessians.

  function Symbolic
             ( c : Standard_Complex_Numbers.Complex_Number;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Computes the Hessian matrix of the product of the variables in x,
  --   multiplied with the coefficient c, symbolically, for testing purposes.

  function Symbolic
             ( c : Standard_Complex_Numbers.Complex_Number;
               idx : Standard_Integer_Vectors.Vector;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Computes the Hessian matrix of the product of the variables in x,
  --   with indices of the participating variables in the product in idx,
  --   multiplied with the coefficient c, symbolically, for testing purposes.

-- code for the Hessian of one product :

  function Algorithmic
             ( c : Standard_Complex_Numbers.Complex_Number;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Applies algorithmic differentiation to compute the Hessian
  --   of the product of the values in x, with coefficient in c,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- REQUIRED : dim >= 2, otherwise zero Hessian.

-- code for the Hession of one indexed product

  function Algorithmic
             ( c : Standard_Complex_Numbers.Complex_Number;
               idx : Standard_Integer_Vectors.Vector;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Applies algorithmic differentiation to compute the Hessian
  --   of the product of the values in x, with coefficient in c,
  --   and with idx the indices of participating variables,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- REQUIRED : dim >= 2, otherwise zero Hessian.

-- updating the Hessian for an indexed product

  procedure Algorithmic
             ( H : in out Standard_Complex_Matrices.Matrix;
               c : in Standard_Complex_Numbers.Complex_Number;
               idx : in Standard_Integer_Vectors.Vector;
               x : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Applies algorithmic differentiation to compute the Hessian
  --   of the product of the values in x, with coefficient in c,
  --   and with idx the indices of participating variables,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- ON ENTRY :
  --   H       the current Hessian matrix,
  --           initialized with zero if called for the first time.
  --   c       coefficient of the term in the circuit;
  --   idx     index of the participating variables;
  --   x       values for all variables.

  -- ON RETURN :
  --   H       updated Hessian matrix, only for upper triangular part.

  -- REQUIRED : idx'last >= 2, otherwise zero Hessian.

-- updates of the Hessian for any monomial

  procedure Algorithmic
              ( H : in out Standard_Complex_Matrices.Matrix;
                c : in Standard_Complex_Numbers.Complex_Number;
                xps : in Standard_Integer_Vectors.Vector;
                idx : in Standard_Integer_Vectors.Vector;
                fac : in Standard_Integer_Vectors.Vector;
                x : in Standard_Complex_Vectors.Vector;
                pwt : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Applies algorithmic differentiation to compute the Hessian
  --   of the product of the values in x, with coefficient in c,
  --   and with idx the indices of participating variables,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- ON ENTRY :
  --   H       the current Hessian matrix,
  --           initialized with zero if called for the first time.
  --   c       coefficient of the term in the circuit;
  --   xps     exponents of all variables in the monomial;
  --   idx     index of the participating variables;
  --   fac     indices to the variables in the common factor;
  --   x       values for all variables;
  --   pwt     values of higher powers of x to evaluate the common factor.

  -- ON RETURN :
  --   H       updated Hessian matrix, only for upper triangular part.

  -- REQUIRED : idx'last >= fac'last >= 1.

-- wrapper functions :

  function Algorithmic
             ( c : Standard_Complex_Circuits.Circuit;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Hessian matrix of the circuit c at x.

  function Algorithmic
             ( c : Standard_Complex_Circuits.Circuit;
               x : Standard_Complex_Vectors.Vector;
               pwt : Standard_Complex_VecVecs.VecVec )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Hessian matrix of the circuit c at x,
  --   with values of higher powers of x in the power table pwt.

  procedure Test_Product ( dim : in integer32; size : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tests the computation of the Hessian of a product of dim variables,
  --   based on the reverse mode to compute the gradient.
  --   The reverse mode computes forward and backward products.
  --   These products then can be used to compute the Hessian.

  -- REQUIRED : dim > 2.

  procedure Test_Circuit ( dim,nbr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random circuit with dim variables and nbr of terms
  --   to test the computation of the Hessian.

  procedure Standard_Test_Power_Circuit ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random circuit with dim variables, nbr of terms,
  --   and highest power pwr, to test the computation of the Hessian,
  --   in standard double precision.

  procedure DoblDobl_Test_Power_Circuit ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random circuit with dim variables, nbr of terms,
  --   and highest power pwr, to test the computation of the Hessian,
  --   in double double precision.

  procedure TripDobl_Test_Power_Circuit ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random circuit with dim variables, nbr of terms,
  --   and highest power pwr, to test the computation of the Hessian,
  --   in triple double precision.

  procedure QuadDobl_Test_Power_Circuit ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random circuit with dim variables, nbr of terms,
  --   and highest power pwr, to test the computation of the Hessian,
  --   in quad double precision.

  procedure PentDobl_Test_Power_Circuit ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random circuit with dim variables, nbr of terms,
  --   and highest power pwr, to test the computation of the Hessian,
  --   in penta double precision.

  procedure OctoDobl_Test_Power_Circuit ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random circuit with dim variables, nbr of terms,
  --   and highest power pwr, to test the computation of the Hessian,
  --   in octo double precision.

  procedure DecaDobl_Test_Power_Circuit ( dim,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random circuit with dim variables, nbr of terms,
  --   and highest power pwr, to test the computation of the Hessian,
  --   in deca double precision.

  procedure Standard_Run_EvalDiff2
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                s : in Standard_Complex_Circuits.Link_to_System;
                cs : in Standard_Coefficient_Circuits.Link_to_System );

  -- DESCRIPTION :
  --   Generates a random point to evaluate the system s of circuits,
  --   and the coefficient system cs, to compute its Jacobian and vector
  --   of Hessians.  The values computed algorithmically are compared with
  --   the symbolic computations on the system p.

  procedure DoblDobl_Run_EvalDiff2
              ( p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                s : in DoblDobl_Complex_Circuits.Link_to_System );

  -- DESCRIPTION :
  --   Generates a random point to evaluate the system s of circuits,
  --   to compute its Jacobian and vector of Hessians.
  --   The values computed algorithmically are compared with
  --   the symbolic computations on the system p.

  procedure TripDobl_Run_EvalDiff2
              ( p : in TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                s : in TripDobl_Complex_Circuits.Link_to_System );

  -- DESCRIPTION :
  --   Generates a random point to evaluate the system s of circuits,
  --   to compute its Jacobian and vector of Hessians.
  --   The values computed algorithmically are compared with
  --   the symbolic computations on the system p.

  procedure QuadDobl_Run_EvalDiff2
              ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                s : in QuadDobl_Complex_Circuits.Link_to_System );

  -- DESCRIPTION :
  --   Generates a random point to evaluate the system s of circuits,
  --   to compute its Jacobian and vector of Hessians.
  --   The values computed algorithmically are compared with
  --   the symbolic computations on the system p.

  procedure PentDobl_Run_EvalDiff2
              ( p : in PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                s : in PentDobl_Complex_Circuits.Link_to_System );

  -- DESCRIPTION :
  --   Generates a random point to evaluate the system s of circuits,
  --   to compute its Jacobian and vector of Hessians.
  --   The values computed algorithmically are compared with
  --   the symbolic computations on the system p.

  procedure OctoDobl_Run_EvalDiff2
              ( p : in OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                s : in OctoDobl_Complex_Circuits.Link_to_System );

  -- DESCRIPTION :
  --   Generates a random point to evaluate the system s of circuits,
  --   to compute its Jacobian and vector of Hessians.
  --   The values computed algorithmically are compared with
  --   the symbolic computations on the system p.

  procedure DecaDobl_Run_EvalDiff2
              ( p : in DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                s : in DecaDobl_Complex_Circuits.Link_to_System );

  -- DESCRIPTION :
  --   Generates a random point to evaluate the system s of circuits,
  --   to compute its Jacobian and vector of Hessians.
  --   The values computed algorithmically are compared with
  --   the symbolic computations on the system p.

  procedure Standard_Test_EvalDiff2;

  -- DESCRIPTION :
  --   Prompts for a polynomial system and tests the correctness
  --   of the computation of the Jacobian and the Hessians,
  --   in standard double precision.

  procedure DoblDobl_Test_EvalDiff2;

  -- DESCRIPTION :
  --   Prompts for a polynomial system and tests the correctness
  --   of the computation of the Jacobian and the Hessians,
  --   in double double precision.

  procedure TripDobl_Test_EvalDiff2;

  -- DESCRIPTION :
  --   Prompts for a polynomial system and tests the correctness
  --   of the computation of the Jacobian and the Hessians,
  --   in triple double precision.

  procedure QuadDobl_Test_EvalDiff2;

  -- DESCRIPTION :
  --   Prompts for a polynomial system and tests the correctness
  --   of the computation of the Jacobian and the Hessians,
  --   in quad double precision.

  procedure PentDobl_Test_EvalDiff2;

  -- DESCRIPTION :
  --   Prompts for a polynomial system and tests the correctness
  --   of the computation of the Jacobian and the Hessians,
  --   in penta double precision.

  procedure OctoDobl_Test_EvalDiff2;

  -- DESCRIPTION :
  --   Prompts for a polynomial system and tests the correctness
  --   of the computation of the Jacobian and the Hessians,
  --   in octo double precision.

  procedure DecaDobl_Test_EvalDiff2;

  -- DESCRIPTION :
  --   Prompts for a polynomial system and tests the correctness
  --   of the computation of the Jacobian and the Hessians,
  --   in deca double precision.

  function Prompt_for_Precision return character;

  -- DESCRIPTION :
  --   Returns '1', '2', '3', '4', '5', '6', or '7', depending
  --   whether double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision is wanted.

  procedure Test_Product_or_Circuit;

  -- DESCRIPTION :
  --   Prompts for the number of terms.  If one, then a product
  --   of variables will be tested, otherwise a number of products
  --   are generated and tested.

  procedure Test_Random_Power_Circuit;

  -- DESCRIPTION :
  --   Prompts for the precision, dimension, number of terms, and 
  --   the highest power to test a randomly generated power circuit.

  procedure Test_EvalDiff2;

  -- DESCRIPTION :
  --   Prompts for the working precision and then tests
  --   for a given polynomial systems.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a dimension and then generates
  --   as many complex random numbers as the dimension.

end Test_Hessian_Circuits;
