with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;

package Tropisms_of_Supports is

-- DESCRIPTION :
--   This package offers tools to compute tropisms for inputs
--   given by supports.

  function Edges ( n : natural; s : Array_of_Lists ) return Array_of_Faces;

  -- DESCRIPTION :
  --   The array on returns contains the edges of the supports s.
  --   The ambient dimension of the supports equals n.

  procedure Show_Parallel_Edges
              ( n : in natural; s : in Array_of_Lists;
                e : in Array_of_Faces );

  -- DESCRIPTION :
  --   Enumerates all pairs of edges and writes to screen those
  --   that are parallel to each other.
  --   This is a first ineffecient algorithm.

  procedure Sort_Parallel_Edges
              ( n : in natural; s : in Array_of_Lists;
                e : in Array_of_Faces );

  -- DESCRIPTION :
  --   Computes parallel edges by sorting the direction vectors.

end Tropisms_of_Supports;
