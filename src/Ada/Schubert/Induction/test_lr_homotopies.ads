with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Complex_Matrices;
with Bracket_Monomials;                  use Bracket_Monomials;
with Intersection_Posets;

package Test_LR_Homotopies is

-- DESCRIPTION :
--   Interactive development of the application of Littlewood-Richardson
--   homotopies to solve general Schubert problems.

  procedure Read_Brackets ( bm : out Bracket_Monomial );

  -- DECRIPTION :
  --   Prompts the user for a bracket monomial.

  procedure Walk_from_Root_to_Leaves
              ( n,k : in integer32; bm : in Bracket_Monomial );

  -- DESCRIPTION :
  --   Prompts the user for m intersection conditions on k-planes in n-space,
  --   and writes the expansion to resolve the conditions from the top down,
  --   walking the intersection poset from the root to the leaves.

  -- ON ENTRY :
  --   n        ambient space
  --   k        dimension of the solution planes;
  --   bm       product of k-brackets, with conditions on the k-planes.

  procedure Connect_Checker_Posets_to_Count
              ( pl : in Intersection_Posets.Poset_List;
                nd : in Intersection_Posets.Poset_Node );

  -- DESCRIPTION :
  --   Connects the root counts at the root of the child poset with
  --   those leaves of the parent poset for which the conditions match.
  --   The main purpose of the connection is to count the roots
  --   bottom up, from the leaves to the root.

  -- ON ENTRY :
  --   pl       list of checker posets at some level of the parent nodes
  --            to the node nd in the intersection poset;
  --   nd       poset of the child.

  procedure Walk_from_Leaves_to_Root
              ( n,k : in integer32; bm : in Bracket_Monomial );

  -- DESCRIPTION :
  --   Prompts the user for m intersection conditions on k-planes in n-space,
  --   and writes the evolution of the root count from the leaves to the root.

  -- ON ENTRY :
  --   n        ambient space
  --   k        dimension of the solution planes;
  --   bm       product of k-brackets, with conditions on the k-planes.

  procedure Stubbing_Littlewood_Richardson
              ( n,k : in integer32; bm : in Bracket_Monomial );

  -- DESCRIPTION :
  --   This procedure shows the progress of the application of Littlewood-
  --   Richardson homotopies to solve general Schubert problems.
  --   Prompts the user for m intersection conditions on k-planes in n-space,
  --   and writes the evolution of the root count from the leaves to the root.

  -- ON ENTRY :
  --   n        ambient space
  --   k        dimension of the solution planes;
  --   bm       product of k-brackets, with conditions on the k-planes.

  procedure Make_Solution_Poset
              ( n,k : in integer32; bm : in Bracket_Monomial );

  -- DESCRIPTION :
  --   Creates the solution poset corresponding to the intersection
  --   poset for the given Schubert conditions.

  -- ON ENTRY :
  --   n        ambient space
  --   k        dimension of the solution planes;
  --   bm       product of k-brackets, with conditions on the k-planes.

  procedure Resolve_Schubert_Problem
              ( n,k : in integer32; bm : in Bracket_Monomial );

  -- DESCRIPTION :
  --   Prompts the user for m intersection conditions on k-planes in n-space,
  --   and writes the evolution of the root count from the leaves to the root.

  -- ON ENTRY :
  --   n        ambient space
  --   k        dimension of the solution planes;
  --   bm       product of k-brackets, with conditions on the k-planes.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu a prompts for a test.

end Test_LR_Homotopies;
