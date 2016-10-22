with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer64_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Standard_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random standard numbers.
--   Either the seed remains hidden for the user, or the user can manage
--   the seed for independent and/or reproducible sequences of numbers.

  function Random_Vector ( first,last,low,upp : integer32 )
                         return Standard_Integer_Vectors.Vector;
  function Random_Vector ( first,last : integer32; low,upp : integer64 )
                         return Standard_Integer64_Vectors.Vector;

  -- DESCRIPTION : 
  --   Returns a vector of range first..last with entries 
  --   generated between low and upp.

  procedure Random_Vector
              ( seed : in out integer32; low,upp : in integer32;
                v : out Standard_Integer_Vectors.Vector );
  procedure Random_Vector
              ( seed : in out integer32; low,upp : in integer64;
                v : out Standard_Integer64_Vectors.Vector );

  -- DESCRIPTION :
  --   With the given seed, generates a vector v of random integer
  --   numbers in [lower, upper].  The seed is updated so the next
  --   call to Random_Vector will give a vector of different numbers.

  function Random_Vector ( first,last : integer32 )
                         return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random floating numbers
  --   uniformly distributed in the interval [-1,+1].

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Given a seed, generates a vector v with random doubles in [-1, +1].
  --   The seed is updated so the next call returns a different vector.

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

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Given a seed, generates a vector v with random complex numbers
  --   with modulus one.  The seed is updated for future use.

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random complex numbers,
  --   with modulus in [10^(-m),10^(+m)].

end Standard_Random_Vectors;
