with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Integer_VecVecs;
with Multprec_Floating_VecVecs;
with Multprec_Complex_VecVecs;

package Multprec_Random_VecVecs is

-- DESCRIPTION :
--   Offers routines to generate vecvecs of random multi-precision numbers.

  function Random_VecVec ( n,m,sz : natural32 )
                         return Multprec_Integer_VecVecs.VecVec;

  -- DESCRIPTION : 
  --   Returns m vectors of range 1..n with entries of size sz.
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m,sz : natural32 )
                         return Multprec_Floating_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with entries of size sz.
  --   The structure on return has range 1..m.

  function Random_VecVec ( n,m,sz : natural32 )
                         return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns m vectors of range 1..n with entries of size sz.
  --   The structure on return has range 1..m.

end Multprec_Random_VecVecs;
