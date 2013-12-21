with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_VecVecs;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;

package Standard_Random_VecVecs is

-- DESCRIPTION :
--   Offers routines to generate vecvecs of random standard numbers.

  function Random_VecVec ( n,m : natural32; low,upp : integer32 )
                         return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION : 
  --   Returns m vectors of range 1..n with entries between low and upp.
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m : natural32 )
                         return Standard_Floating_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random floating numbers,
  --   where the numbers in the vectors belong to [-1,+1].
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m,g : natural32 )
                         return Standard_Floating_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random floating numbers,
  --   where the numbers have magnitude in [10^(-g),10^(+g)].
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m : natural32 )
                         return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random complex numbers
  --   on the complex unit circle.
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m,g : natural32 )
                         return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with random complex numbers
  --   of magnitude in the interval [10^(-g),10^(+g)].
  --   The structure on return has range 1..m.

end Standard_Random_VecVecs;
