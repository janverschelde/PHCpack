with Symmetry_Group;                     use Symmetry_Group;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Linear_Symmetric_Reduction is

-- DESCRIPTION :
--   This package contains two routines that enable the faster
--   solution of a symmetric product system, by extracting the
--   generating positions.

  function Linear_Symmetric_Reduce ( sign : boolean ) return List;
  function Linear_Symmetric_Reduce ( v,w : List_of_Permutations ) return List;

  -- DESCRIPTION :
  --   Returns the generating list of positions in the random product system.

  -- REQUIRED : data in package Random_Product_System is not empty.

  procedure Linear_Symmetric_Reduce ( lp : in out List; sign : in boolean );
  procedure Linear_Symmetric_Reduce 
               ( v,w : in List_of_Permutations; lp : in out List );

  -- DESCRIPTION :
  --   Given a (G,V,W)-symmetric product system, a list of positions
  --   that indicate the generating subsystems will be returned.

  -- REQUIRED : data in package Random_Product_System is not empty.

  -- ON ENTRY :
  --   v,w       group representations,
  --             if not provided, then the full permutation group is assumed;
  --   sign      if true, then there is also sign symmetry to consider;
  --   lp        a list of positions, indicating the useful
  --             linear systems in the product system.

  -- ON RETURN :
  --   lp        a sublist of the former list, 
  --             contains only the generating linear systems.

end Linear_Symmetric_Reduction;
