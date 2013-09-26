with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;

package Parallel_Edges is

-- DESCRIPTION :
--   This package provides routines to compute parallel edges
--   of a tuple of edges of polytopes.

  function Edge_Direction
              ( n : natural; f : Face )
              return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the normalized direction vector defined by the vertices
  --   that span the edge f in n-space.  The normalized direction vector
  --   on return has its last nonzero component equal to one.

  function Edge_Directions
              ( n : natural; f : Faces )
              return Standard_Floating_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns an array of normalized directions in n-space, 
  --   defined by the faces in the list f. 
  --   The array on return is of range 1..Length_Of(f).

  function Edge_Directions
              ( n : natural; f : Array_of_Faces )
              return Standard_Floating_VecVecs.Array_of_VecVecs;

  -- DESCRIPTION :
  --   Returns an array of normalized directions in n-space, 
  --   defined by the faces in the lists f. 
  --   The array on return has the same range as f.

  function Sort ( v : Standard_Floating_VecVecs.VecVec )
                return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of indices to v that sorts the vectors
  --   increasingly in lexicographical order.
  --   The range of the vector on return is the same as the range of v.

  function Sort ( v : Standard_Floating_VecVecs.Array_of_VecVecs )
                return Standard_Natural_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Sorts all vectors in the components of v separately
  --   and returns for each component the index vector that puts the
  --   vectors in increasing lexicographic order.
  --   The vector on return has the same range as v.

  procedure Merge_Sort
              ( v,w : in Standard_Floating_VecVecs.VecVec;
                sv,sw : in Standard_Natural_Vectors.Vector;
                vnb,wnb : in natural;
                snb,dst : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Applies merge sort to the directions in v and w.

  -- REQUIRED :
  --   snb'range = dst'range = 1..v'length+w'length

  -- ON ENTRY :
  --   v        normalized direction vectors to edges of support vnb;
  --   w        normalized direction vectors to edges of support wnb;
  --   sv       labels to sort the directions in v lexicographically;
  --   sw       labels to sort the directions in w lexicographically;
  --   vnb      number of the support for the directions of v;
  --   wnb      number of the support for the directions of w.

  -- ON RETURN :
  --   snb      support labels to sort the directions, snb(i) is either
  --            vnb or wnb: if i-th direction comes from v or w;
  --   dst      direction labels to sort, dst(i) is index to a direction
  --            in v or w depending on what the corresponding snb(i) is.

  function Sum_of_Lengths
              ( v : in Standard_Floating_VecVecs.Array_of_VecVecs ) 
              return natural;

  -- DESCRIPTION :
  --   Returns the sum of all lengths of the vectors in v.

  procedure Merge_Sort
              ( v : in Standard_Floating_VecVecs.Array_of_VecVecs;
                sv : in Standard_Natural_VecVecs.VecVec;
                snb,dst : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Applies merge sort to the directions in v.

  -- REQUIRED :
  --   snb'range = dst'range = 1..Sum_of_Lengths(v).

  -- ON ENTRY :
  --   v        normalized direction vectors for a tuple of edges;
  --   sv       labels to sort the directions in each v lexicographically,.
  --            obtained via the sort above, i.e: sv = Sort(v).

  -- ON RETURN :
  --   snb      sorted support numbers, snb(i) is in v'range, 
  --            indicating from which support the i-th direction comes;
  --   dst      labels to sort the directions, dst(i) is an index to
  --            a direction in v(snb(i)).

end Parallel_Edges;
