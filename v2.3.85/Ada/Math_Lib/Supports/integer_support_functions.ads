with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Integer_Support_Functions is

-- DESCRIPTION :
--   This package provides support functions for polytopes spanned by
--   lists of integer vectors.

  function Maximal_Support ( L : List; v : Vector ) return integer32;
  function Minimal_Support ( L : List; v : Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns respectively max <d,v> or min <d,v>, for all d in L.

  procedure Min_Max ( L : in List; k : in integer32;
                      min,max : in out integer32 );

  -- DESCRIPTION :
  --   Computes the minimum and maximum of the k-th entry of all
  --   points belonging to L.

  function Graded_Max ( L : List ) return Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the maximal element in the list w.r.t. the graded
  --   lexicographical ordening.

  -- REQUIRED : List is not empty.

  function Face ( L : List; v : Vector; m : integer32 ) return List;

  -- DESCRIPTION :
  --   Returns a list of vectors d for which <d,v> = m.
  --   For m = Maximal_Support(L,v), we get the face w.r.t. the outer normal v;
  --   for m = Minimal_Support(L,v), we get the face w.r.t. the inner normal v.

  function Inner_Face ( L : List; v : Vector ) return List;
  function Outer_Face ( L : List; v : Vector ) return List;

  -- DESCRIPTION :
  --   Returns the face of the list L, where v is inner or outer normal.

end Integer_Support_Functions;
