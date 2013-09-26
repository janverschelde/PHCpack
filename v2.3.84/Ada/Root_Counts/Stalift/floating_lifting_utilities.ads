with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;

package Floating_Lifting_Utilities is

-- DESCRIPTION :
--   This package provides some utilities for dealing with lifting functions.

  function Adaptive_Lifting ( L : Array_of_Lists ) return Vector;

  -- DESCRIPTION :
  --   Returns upper bounds for a random lifting, depending on the lengths
  --   of the lists in L.

  procedure Search_Lifting ( L : in List; pt : in Vector;
                             found : out boolean; lif : out double_float );

  -- DESCRIPTION :
  --   Searches the lifting of the point in the lifted list L.
  --   If found, then lif equals the lifting, otherwise lif is meaningless.

  function Search_and_Lift ( L : List; pt : Vector ) return Vector;

  -- DESCRIPTION :
  --   Given a lifted list of points and a unlifted vector, the function
  --   either returns the corresponding lifted vector from the list, or
  --   the same point, when there is no lifted point in l whose projection
  --   equals the given point pt.

  function Occurred_Lifting ( mixsub : Mixed_Subdivision; k : integer32;
                              pt : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the point pt augmented with its lifting value when it
  --   occurs in the kth component of some cell in the subdivision.
  --   Otherwise the point itself is returned.

  function Occurred_Lifting
               ( n : integer32; mix : Standard_Integer_Vectors.Vector;
                 points : Array_of_Lists; mixsub : Mixed_Subdivision )
               return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns the lifted points for those points that belong to the
  --   mixed cells in the subdivision.

  function Induced_Lifting ( mixsub : Mixed_Subdivision; k : integer32;
                             pt : Vector ) return Vector;
  function Induced_Lifting
               ( n : integer32; mix : Standard_Integer_Vectors.Vector;
                 points : Array_of_Lists; mixsub : Mixed_Subdivision )
               return Array_of_Lists;

  -- DESCRIPTION :
  --   Given a mixed subdivision for a tuple of supports,
  --   then the lifted points will be returned as induced by the
  --   subdivision.   When points do not occur in the mixed subdivision,
  --   they will be lifted conservatively.

  function Conservative_Lifting 
               ( mic : Mixed_Cell; k : integer32; point : Vector )
               return double_float;
  function Conservative_Lifting 
               ( mixsub : Mixed_Subdivision; k : integer32; point : Vector )
               return double_float;
  
  -- DESCRIPTION :
  --   Returns the value of the conservative lifting function of the point
  --   to be considered w.r.t. the kth polytope.

  -- REQUIRED :
  --   The given point must already be in the lifted space and its last
  --   coordinate must contain already a lower bound for the lifting value.

  function Lifted_Supports ( r : integer32; mixsub : Mixed_Subdivision )
                           return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns an array of lists of range 1..r with the lifted supports
  --   for all points in the given mixed subdivision.

  -- REQUIRED : r = c.pts'last for all cells c in mixsub.

end Floating_Lifting_Utilities;
