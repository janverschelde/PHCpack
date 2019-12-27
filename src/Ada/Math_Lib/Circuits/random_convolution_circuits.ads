with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Random_Convolution_Circuits is

-- DESCRIPTION :
--   The functions in this package return random convolution circuits,
--   mainly for testing and benchmarking purposes.

  function Random_Exponents
             ( dim,nbr : integer32; pwr : integer32 := 1 ) 
             return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Generates as many random vectors of dimension dim
  --   as the value of nbr.  All vectors have a nonzero sum.
  --   The higher power equals pow.

  function Standard_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return Standard_Speelpenning_Convolutions.Convolution_Circuit;
  function DoblDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuit;
  function QuadDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuit;

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and coefficients
  --   in double, double double, or quad double precision,
  --   and returns a convolution circuit.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  function Standard_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return Standard_Speelpenning_Convolutions.Convolution_Circuits;
  function DoblDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuits;
  function QuadDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuits;

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and coefficients
  --   in double, double double, or quad doble precision,
  --   and returns as many convolution circuits as the value of dim.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

end Random_Convolution_Circuits;
