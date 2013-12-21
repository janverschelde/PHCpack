with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Double_Double_VecVecs;
with DoblDobl_Complex_VecVecs;

package DoblDobl_Random_VecVecs is

-- DESCRIPTION :
--   Offers routines to generate vecvecs of random double double numbers.

  function Random_VecVec ( n,m : natural32 )
                         return Double_Double_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random double double numbers
  --   in the interval [-1,+1].
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m,g : natural32 )
                         return Double_Double_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random double double numbers
  --   with absolute value in the range [10^(-g),10^(+g)].
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m : natural32 )
                         return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random complex double doubles
  --   on the unit circle.
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m,g : natural32 )
                         return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random complex double doubles
  --   with modulus in the range [10^(-g),10^(+g)].
  --   The structure on return has range 1..m.

end DoblDobl_Random_VecVecs;
