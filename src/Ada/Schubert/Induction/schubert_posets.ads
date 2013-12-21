with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Multprec_Natural_Numbers;          use Multprec_Natural_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Bracket_Monomials;                 use Bracket_Monomials;
with Checker_Posets;                    use Checker_Posets;
with Intersection_Posets;               use Intersection_Posets;

package Schubert_Posets is

-- DESCRIPTION :
--   A Schubert poset looks at an intersection poset starting at the leaves
--   to implement the generalization of the flags to solve a Schubert problem.
--   The Specialize procedures below extend the Intersect procedures of
--   the package Intersection_Posets.

  procedure Initialize_Leaves_Coefficients ( ps : in out Poset );

  -- DESCRIPTION :
  --   Initializes the Littlewood-Richardson coefficient at all leaves
  --   of the checker poset ps to one.

  procedure Initialize_Leaves_Coefficients
               ( ps : in out Poset; pl : in Poset_List );

  -- DESCRIPTION :
  --   Initializes the Littlewood-Richardson coefficients at all leaves
  --   of the checkers poset ps with the coefficients at the roots of
  --   those posets of pl that can be considered as children.

  procedure Specialize_Leaves ( ps : in out Poset; w : in Vector );

  -- DESCRIPTION :
  --   Adjusts root counts stored in the Littlewood-Richardson coefficients
  --   in the leaves of the poset ps, taking into account the condition in w.

  procedure Count_Roots ( ps : in out Poset );

  -- DESCRIPTION :
  --   Assigns the Littlewood-Richardson coefficients of all nodes in the
  --   poset ps, except for the coefficients at the leaves.  On return,
  --   the coefficient at the root node of ps equals the root count.

  procedure Specialize ( ips : in out Intersection_Poset;
                         pnd : in Link_to_Poset_Node;
                         w : in Vector; fail : out boolean );

  -- DESCRIPTION :
  --   Starting at the leaves of the current poset pnd of ips, the
  --   intersection poset ips is updated with the intersection condition w
  --   and new posets are added, after checking for duplicates.
  --   On return fail is true, if w made unhappy combinations with
  --   all intersection conditions at the leaves of the poset pnd.

  -- REQUIRED : ips.level < ips.m.

  procedure Specialize ( ips : in out Intersection_Poset;
                         w : in Vector; fail : out boolean );

  -- DESCRIPTION :
  --   Adds one intersection condition w to the Schubert problem with
  --   its partial resolution already stored in ips.  If the updated
  --   Schubert problem admits no solution, then fail is true on return.

  -- REQUIRED : ips.level < ips.m.

  function Specialize ( n,k : integer32; b : Bracket_Monomial )
                      return Intersection_Poset;

  -- DESCRIPTION :
  --   Returns an intersection posets based on the intersection conditions
  --   on k-planes in n-space, stored in the brackets of b.
  --   If ips.level on return equals ips.m, then the Schubert problem admits
  --   solutions and the root count is in the leaves.

  -- REQUIRED : number of brackets in b must be at least 3.

  function Root_Count_at_Leaves
              ( ips : Intersection_Poset ) return Natural_Number;

  -- DESCRIPTION :
  --   Returns the number of solution collected at the leaves of
  --   all posets at the leaf levels of ips.

  procedure Count_Roots ( ips : in Intersection_Poset );

  -- DESCRIPTION :
  --   Applies the Count_Roots procedure to all posets in ips,
  --   starting at the leaves.  If all goes well, the root count
  --   is then stored at the root of the poset at the root of ips.

end Schubert_Posets;
