with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;

package QuadDobl_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random double double numbers.

  function Random_Vector ( first,last : integer32 )
                         return Quad_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random quad doubles
  --   in the interval [-1,+1].

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
  
  function Random_Vector ( first,last : integer32; m : natural32 )
                         return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random 
  --   quad double complex numbers with modulus in [10^(-m),10^(+m)].

end QuadDobl_Random_Vectors;
