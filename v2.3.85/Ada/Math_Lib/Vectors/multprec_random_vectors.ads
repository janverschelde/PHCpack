with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Vectors;
with Multprec_Floating_Vectors;
with Multprec_Complex_Vectors;

package Multprec_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random multi-precision numbers.

  function Random_Vector ( first,last : integer32; sz : natural32 )
                         return Multprec_Integer_Vectors.Vector;

  -- DESCRIPTION : 
  --   Returns a vector of range first..last with entries of size sz.

  function Random_Vector ( first,last : integer32; sz : natural32 )
                         return Multprec_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with entries of size sz.

  function Random_Vector ( first,last : integer32; sz : natural32 )
                         return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with entries of size sz.

end Multprec_Random_Vectors;
