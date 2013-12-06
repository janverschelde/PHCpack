with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Quad_Double_VecVecs;
with QuadDobl_Complex_VecVecs;

package QuadDobl_Random_VecVecs is

-- DESCRIPTION :
--   Offers routines to generate vecvecs of random quad double numbers.

  function Random_VecVec ( n,m : natural32 )
                         return Quad_Double_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random quad double numbers
  --   in the interval [-1,+1].
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m,g : natural32 )
                         return Quad_Double_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random quad double numbers
  --   with absolute value in the range [10^(-g),10^(+g)].
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m : natural32 )
                         return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random complex quad double
  --   numbers on the unit circle.
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m,g : natural32 )
                         return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random complex quad double
  --   numbers with modulus in the interval [10^(-g),10^(+g)].
  --   The structure on return has range 1..m.

end QuadDobl_Random_VecVecs;
