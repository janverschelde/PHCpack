with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;

package Exponent_Indices is

-- DESCRIPTION :
--   An exponent index stores the indices of the nonzero entries
--   in an exponent vector.  The indices of the higher powers in
--   exponent vectors are stored in the factor index.

  function Index_Size ( v : Standard_Integer_Vectors.Vector ) 
                      return integer32;
  function Index_Size ( v : Standard_Integer_Vectors.Link_to_Vector ) 
                      return integer32;

  -- DESCRIPTION :
  --   Returns the number of positive entries in v.

  function Factor_Size ( v : Standard_Integer_Vectors.Link_to_Vector ) 
                       return integer32;

  -- DESCRIPTION :
  --   Returns the number of entries in v that are larger than one.

  function Exponent_Index
             ( xp : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector;
  function Exponent_Index
             ( xp : Standard_Integer_Vectors.Link_to_Vector )
             return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the vector of positions k in xp for which xp(k) = 1.

  function Exponent_Index
             ( xp : Standard_Integer_VecVecs.VecVec )
             return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the vector of exponent indices for xp.

  function Factor_Index
             ( xp : Standard_Integer_Vectors.Link_to_Vector )
             return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the factor index of the exponent vector.
  --   If none of the entries in xp is larger than one,
  --   then null is returned.
  --   Otherwise the vector on return contains all indices to
  --   the exponets in v that are larger than one.

  function Factor_Index
             ( xp : Standard_Integer_VecVecs.VecVec )
             return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   For each exponent vector, the factor index is computed and returned
  --   in the vector of vectors.  Note that the entries of the vector
  --   on return may be null if there are no powers higher than one.

  function Maxima ( xp : Standard_Integer_VecVecs.VecVec )
                  return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the largest power for each value in the exponents.

  function Polynomial_Degree
             ( xp : Standard_Integer_Vectors.Link_to_Vector )
             return integer32;

  -- DESCRIPTION :
  --   Returns -1 if xp is null,
  --   otherwise returns the sum of all elements in xp.

  function Polynomial_Degree
             ( xp : Standard_Integer_VecVecs.VecVec ) return integer32;
        
  -- DESCRIPTION :
  --   For the exponents in xp, returns the largest degree,
  --   as this is the degree of the polynomial with exponents in xp.

end Exponent_Indices;
