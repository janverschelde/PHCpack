with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Lexicographical_Supports is

-- DESCRIPTION :
--   Ordering supports in lexicographical increasing order
--   is convenient for evaluation.

  function LexLess ( a,b : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if a is lexicographically less than b.
  --   Returns false otherwise.

  -- REQUIRED : a'range = b'range.

  function DegLess ( a,b : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if a has lower degree than b, using lexicographical
  --   order as tie breaker.  Returns false otherwise.

  -- REQUIRED : a'range = b'range.

  procedure Swap ( a,b : in out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Swaps the elements in a with b.

  -- REQUIRED : a'range = b'range.

  function Sort ( s : List ) return List;

  -- DESCRIPTION :
  --   The list on return has the same content as s,
  --   sorted in degree lexicographically increasing order.

  function Index ( s : Standard_Integer_VecVecs.VecVec;
                   v : Standard_Integer_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Applies binary search to return the index of a vector
  --   in the degree lexicographically sorted vector of vectors s.
  --   Returns either 0 if v does not belong to s
  --   or returns k with Equal(s(k).all,v).

-- FACTORED REPRESENTATION :

  function First_Positive ( s : List ) return natural32;

  -- DESCRIPTION :
  --   Returns the index of the first vector in s that has
  --   no negative components and at least one component positive.
  --   Returns 0 if s is empty or there is no such first vector.
  --   What is returned serves as "start" in Factor below.

  procedure Factor ( v : in out Standard_Integer_VecVecs.VecVec;
                     org : in out Standard_Integer_VecVecs.VecVec;
                     start,k : in integer32 );

  -- DESCRIPTION :
  --   Updates the factored representation of v at the k-th vector in v,
  --   searching for the lexicographically largest factor in start..(k-1)
  --   that divides v(k).  The original vectors are in the org parameter.

  function Factor ( s : List ) return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Given on input a lexicographically sorted list of vectors of range 1..n
  --   returns the list in factored form, of vectors of range 0..n.
  --   For a vector v in the vector of vectors on return, we have:
  --     v(0) = 0 : no vector divides v and v(1..n) contains the original
  --     v(0) = k : then the k-th vector divides v and v(1..n) contains
  --                the other factor that when multiplied gives the original.

  -- REQUIRED :
  --   The list s is the output of the Sort from above
  --   and all components of all elements of s are positive.

-- COMPRESSION :

  function Compress ( v : Standard_Integer_Vectors.Vector )
                    return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   A compressed monomial vector is a vector of pairs:
  --   odd entries contain the index of the variable,
  --   followed by the power of the variable.
  --   Because of sparsity, the vector on return is often much shorter.

  function Compress ( v : Standard_Integer_VecVecs.VecVec )
                    return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Applies the function Compress from above to all vectors in v.

-- PRE-FACTORIZATON :

  function Decrement ( v : Standard_Integer_Vectors.Vector ) 
                     return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a decremented vector v, subtracting from every component one
  --   and setting the component to zero if the result would be negative.
  --   This decremented vector potentially becomes an auxiliary node in
  --   an arithmetic network for efficient evaluation in factored form.

  function Nodes ( s : Standard_Integer_VecVecs.VecVec ) return List;

  -- DESCRIPTION :
  --   Given a sorted support in s, returns a list of nodes not already in s.
  --   Note that the list on return may be empty.
  --   A node for a vector v is the vector with all components decremented
  --   by one or zero.

end Lexicographical_Supports;
