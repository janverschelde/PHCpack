with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;

package Parallel_Directions is

-- DESCRIPTION :
--   The functions in this package are mainly auxiliary to the operations
--   in Parallel_Edges.

  function Pivot ( v : Standard_Integer_Vectors.Vector ) return integer;

  -- DESCRIPTION :
  --   Returns v'first - 1 if v is zero, otherwise it returns the index
  --   in v of the last nonzero entry.

  function Normalize ( v : Standard_Integer_Vectors.Vector; p : integer )
                     return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   The vector on return has 1 at its position p.
  --   Every element in the vector on return is divided by v(p).

  -- REQUIRED : p is in v'range and v(p) /= 0.

  function ">" ( v,w : Standard_Floating_Vectors.Link_to_Vector )
               return boolean;

  -- DESCRIPTION :
  --   Returns true at the first entry of v that is larger than the
  --   corresponding entry in w.  Returns false if v and w are equal
  --   or if at the first different position v is smaller than w.
  --   By different, the difference must be larger than 1.0E-8.

  function Equal ( v,w : Standard_Floating_Vectors.Link_to_Vector )
                 return boolean;

  -- DESCRIPTION :
  --   Two directions are equal if all their components differ
  --   by less than 1.0E-8.

  procedure Swap ( v : in out Standard_Natural_Vectors.Vector;
                   i,j : in integer );

  -- DESCRIPTION :
  --   Swaps the entries i and j in the vector v.

  function Min ( v : Standard_Floating_VecVecs.Array_of_VecVecs;
                 sv : Standard_Natural_VecVecs.VecVec;
                 ind : Standard_Natural_Vectors.Vector ) return integer;

  -- DESCRIPTION :
  --   Returns the index of the current minimum vector in v, sorted as sv,
  --   with ind pointing at the current position in v.

  -- ON ENTRY :
  --   v         normalized direction vectors for tuple of supports;
  --   sv        labels to sort the vectors increasingly lexicographically;
  --   ind       index vector to current position in v.

  -- ON RETURN :
  --   Returns -1 if ind(i) for all i is larger than v(i)'last,
  --   otherwise, returns the index of the minimal vector.

  function Last_Index ( v : Standard_Floating_VecVecs.Array_of_VecVecs;
                        ind : Standard_Natural_Vectors.Vector )
                      return integer;

  -- DESCRIPTION :
  --   Returns -1 if there are two or more indices k in ind for which
  --   ind(k) <= v(k)'last.  Otherwise the index k of the last entry 
  --   in ind for which ind(k) <= v(k)'last still holds is returned.
  --   If for supports ind(k) > v(k)'last, then 0 is returned.

end Parallel_Directions;
