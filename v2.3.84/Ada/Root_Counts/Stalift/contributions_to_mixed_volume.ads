with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;

package Contributions_to_Mixed_Volume is

-- DESCRIPTION :
--   This package provides utilities to determine whether a point
--   contributes to the mixed volume of a tuple of point sets, i.e.: whether 
--   it can be removed of the tuple without decreasing the mixed volume.

-- THE CRITERION : (simple/exhaustive)

  function Simple_Zero_Contribution 
               ( pts : Array_of_Lists;
                 x : Vector; i : integer32 ) return boolean;

  function Simple_Zero_Contribution
               ( pts : Array_of_Lists; ifacets : Faces;
                 x : Vector; i : integer32 ) return boolean;

  function Exhaustive_Zero_Contribution
               ( pts : Array_of_Lists;
                 x : Vector; i : integer32 ) return boolean;

  function Exhaustive_Zero_Contribution
               ( pts : Array_of_Lists; ifacets : Faces;
                 x : Vector; i : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the point x does not influence the mixed volume of
  --   the tuple pts, when considered w.r.t. the ith component.
  --   There is a simple and an exhaustive criterion to check the contribution
  --   of the point.  When the simple criterion returns false, the point does
  --   not necessarily make a contribution to the mixed volume, whereas when
  --   the exhaustive criterion returns false, we can be sure that the point
  --   contributes to the mixed volume.

  -- ON ENTRY :
  --   pts       tuple of points sets, with dimension = x'length;
  --   ifacets   all facets of the polytope spanned by pts(i), that contain x;
  --   x         a point, considered w.r.t. pts(i);
  --   i         component of pts, i is in pts'range.

  -- ON RETURN :
  --   true      if Mixed_Volume(pts) = Mixed_Volume(pts, pts(i) contains x);
  --   false     when the criterion fails to decide.

-- SWEEPING THROUGH THE POINT LISTS : (with simple/exhaustive criterion)

  function Simple_Sweep ( pts : Array_of_Lists ) return Array_of_Lists;

  function Simple_Sweep ( pts : Array_of_Lists; facets : Array_of_Faces )
                        return Array_of_Lists;

  function Exhaustive_Sweep ( pts : Array_of_Lists ) return Array_of_Lists;

  function Exhaustive_Sweep ( pts : Array_of_Lists; facets : Array_of_Faces )
                            return Array_of_Lists;

  -- DESCRIPTION :
  --   Applies the simple/exhaustive criterion to the points.
  --   Returns those points for which the criterion has failed to 
  --   demonstrate that the contribution to the mixed volume was zero.

  -- ON ENTRY :
  --   pts       tuple of point sets;
  --   facets    tuple of facets of the polytopes spanned by pts,
  --             facets'range = pts'range.

  -- ON RETURN :
  --   A tuple of point sets where the points that do not contribute to the 
  --   mixed volume have been eliminated, according to the simple or exhaustive
  --   criterion that has been used.
  --   Note that essential sets are not necessarily unique, because they depend
  --   on the order in which points with zero contribution are eliminated.

end Contributions_to_Mixed_Volume;
