with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer64_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Standard_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random standard numbers.

  function Random_Vector ( first,last,low,upp : integer32 )
                         return Standard_Integer_Vectors.Vector;
  function Random_Vector ( first,last : integer32; low,upp : integer64 )
                         return Standard_Integer64_Vectors.Vector;

  -- DESCRIPTION : 
  --   Returns a vector of range first..last with entries 
  --   generated between low and upp.

  function Random_Vector ( first,last : integer32 )
                         return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random floating numbers
  --   uniformly distributed in the interval [-1,+1].

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random floating numbers
  --   with absolute values in [10^(-m),10^(+m)].

  function Random_Vector ( first,last : integer32 )
                         return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random complex numbers,
  --   on the complex unit circle.

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random complex numbers,
  --   with modulus in [10^(-m),10^(+m)].

end Standard_Random_Vectors;
