with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Intersection_Posets;               use Intersection_Posets;

package Intersection_Posets_io is

-- DESCRIPTION :
--   This package provides output routines for intersection posets
--   to resolve a general Schubert problem with checker games.

  procedure Write_Parents ( pl : in Poset_List; nd : in Poset_Node );

  -- DESCRIPTION :
  --   Enumerates all parents from the list pl of the node nd
  --   and writes for each parent one line: white checkers at
  --   the root and white checkers at the corresponding leaf.

  procedure Write_Formal_Equations
               ( ips : in Intersection_Poset; k : in integer32 );
  procedure Write_Formal_Equations
               ( file : in file_type;
                 ips : in Intersection_Poset; k : in integer32 );

  -- DESCRIPTION :
  --   Writes the formal equations for level k in the poset.

  -- REQUIRED : 1 <= k <= ips.level.

  procedure Write_Formal_Equations ( ips : in Intersection_Poset );
  procedure Write_Formal_Equations
              ( file : in file_type; ips : in Intersection_Poset );

  -- DESCRIPTION :
  --   Writes a formal equation for each level in the poset,
  --   either to standard output or to file.

  procedure Write_Lefthand_Product
              ( ips : in Intersection_Poset; k : in integer32 );
  procedure Write_Lefthand_Product
              ( file : in file_type;
                ips : in Intersection_Poset; k : in integer32 );

  -- DESCRIPTION :
  --   Writes the lefthand product of intersection conditions resolved
  --   by the intersection poset, starting at level k.
  --   If k equals 1, then the writing starts at the root of the poset.

  procedure Write_Lefthand_Product ( ips : in Intersection_Poset );
  procedure Write_Lefthand_Product
              ( file : in file_type; ips : in Intersection_Poset );

  -- DESCRIPTION :
  --   Writes the complete list of intersection conditions as a product
  --   of brackets as defined by the first nodes at each level.

  procedure Write_Final_Sum ( pl : in Poset_List );
  procedure Write_Final_Sum ( file : in file_type; pl : in Poset_List );

  -- DESCRIPTION :
  --   Writes the formal sum of intersection conditions at the leaves
  --   of all the posets in the list of poset in pl.

  procedure Write_Expansion ( ips : in Intersection_Poset );
  procedure Write_Expansion
              ( file : in file_type; ips : in Intersection_Poset );

  -- DESCRIPTION :
  --   Writes the expansion of the intersection conditions as computed
  --   in the intersection poset, to standard output or to file.

end Intersection_Posets_io;
