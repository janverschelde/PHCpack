with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;

package Arrays_of_Lists_Utilities is

-- DESCRIPTION :
--   This package offers some utilities for working with
--   arrays of lists of integer vectors.

  function All_Equal ( al : Array_of_Lists ) return boolean;

  -- DESCRIPTION :
  --   returns true if all lists in al are the same.

  function Interchange2 ( al : Array_of_Lists ) return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns a permuted array such that if the original array contains
  --   a list of length less than or equal to 2, the permuted array has
  --   this list on its first position.

  function Index2 ( al : Array_of_Lists ) return integer32;

  -- DESCRIPTION :
  --   Returns the index of the list in the array which has a length less
  --   than or equal to 2.  If no such list can be found, then al'first
  --   will be returned.

  procedure Mixture ( al : in Array_of_Lists; perm,mix : out Link_to_Vector );

  -- DESCRIPTION :
  --   Computes the mixture of the sets in al: mix(i) indicates how many
  --   times the ith list, given by al(perm(i)), occurs in the array.
  --   The list with minimal occurence stands in front.

  function Permute ( perm : Vector; al : in Array_of_Lists )
                   return Array_of_Lists;

  -- DESCRIPTION :
  --   The returned array is permuted according to the permutation vector.

  function Different_Points ( al : Array_of_Lists ) return List;

  -- DESCRIPTION :
  --   Returns a list of all different points from the lists
  --   of al, the first list will not be considered.

  function Different_Points ( al : Array_of_Lists ) return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns lists of all different points from the lists
  --   of al, the first list will not be considered.

  procedure Remove_Duplicates ( al : in out Array_of_Lists );

  -- DESCRIPTION :
  --   Except for the first list, all duplicates will be removed from
  --   the lists in al.

  procedure Shift ( al : in out Array_of_Lists; shiftvecs : in VecVec );
  function  Shift ( al : Array_of_Lists; shiftvecs : VecVec )
                  return Array_of_Lists;

  -- DESCRIPTION :
  --   All lists in al will be shifted along the degrees in shiftvecs.

  -- REQUIRED :
  --   The ranges of shiftvecs and al are the same.

  procedure Projection ( al : in Array_of_Lists; v : in Vector;
                         ind : integer32; res : in out Array_of_Lists;
                         degenerate : out boolean );
  -- DESCRIPTION :
  --   The first list in al will not be considered.
  --   after termination:
  --   1. degenerate = false:
  --    res is an array of degrees lists where for all points x in
  --    the list res(i) the following holds: < x , v > = pv, where
  --    pv equals the support of the list res(i) in the direction v.
  --    The ind-entry of each point in the lists has been deleted.
  --    Duplicates have been removed from the lists.
  --   2. degenerate = true:
  --    if there exists an entry i for which Length_Of(res(i)) = 1.

  -- REQUIRED : al'first = res'first and al'last = res'last + 1

end Arrays_of_Lists_Utilities;
