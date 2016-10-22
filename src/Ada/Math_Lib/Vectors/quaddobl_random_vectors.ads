with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;

package QuadDobl_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random quad double numbers.
--   Either the seed remains hidden for the user, or the user can manage
--   the seed for independent and/or reproducible sequences of numbers.

  function Random_Vector ( first,last : integer32 )
                         return Quad_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random quad doubles
  --   in the interval [-1,+1].

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Quad_Double_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random quad doubles,
  --   with each entry in [-1,+1].
  --   The seed is updated so the next call returns a different v.

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return Quad_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random quad doubles
  --   with absolute value in [10^(-m),10^(+m)].  
  
  function Random_Vector ( first,last : integer32 )
                         return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random 
  --   quad double complex numbers on the unit circle.

  procedure Random_Vector
              ( seed : in out integer32;
                v : out QuadDobl_Complex_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random complex numbers
  --   of quad double precision with modulus one.
  --   The seed is updated so the next call returns a different v.
  
  function Random_Vector ( first,last : integer32; m : natural32 )
                         return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random 
  --   quad double complex numbers with modulus in [10^(-m),10^(+m)].

end QuadDobl_Random_Vectors;
