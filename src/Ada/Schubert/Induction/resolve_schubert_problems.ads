with text_io;                            use text_io;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Standard_Complex_VecMats;
with Intersection_Posets;                use Intersection_Posets;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Resolve_Schubert_Problems is

-- DESCRIPTION :
--   We resolve general Schubert problems tracking solution paths as defined
--   by Littewood-Richardson homotopies.

  procedure Initialize_Leaves ( pl : in out Poset_List );

  -- DESCRIPTION :
  --   Initializes the root count at the leaves to one and
  --   at all other nodes in the poset list to zero.

  procedure Initialize_Nodes ( pl : in out Poset_List );

  -- DESCRIPTION :
  --   Initializes the coefficient at every node in pl to zero.

  procedure Connect_Checker_Posets_to_Count
              ( file : in file_type;
                pl : in Poset_List; nd : in Poset_Node );

  -- DESCRIPTION :
  --   Connects the root counts at the root of the child poset with
  --   those leaves of the parent poset for which the conditions match.
  --   The main purpose of the connection is to stub the homotopies
  --   and show the progress of tracking of the paths.

  -- ON ENTRY :
  --   file     for intermediate output;
  --   pl       list of checker posets at some level of the parent nodes
  --            to the node nd in the intersection poset;
  --   nd       poset of the child.

  procedure Count_Roots
              ( file : in file_type; ips : in out Intersection_Poset;
                roco : out Natural_Number );

  -- DESCRIPTION :
  --   Performs a bottom up root count, in the same order as the
  --   resolution with Littlewood-Richardson homotopies proceeds.
  --   This root counter can be viewed as a stub for the resolution
  --   of the intersection conditions by Littlewood-Richardson homotopies.

  -- REQUIRED :
  --   The intersection conditions are processed into an intersection poset,
  --   given in ips on entry.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   ips      an intersection poset built to resolve Schubert conditions.

  -- ON RETURN :
  --   ips      intersection poset with Littlewood-Richardson coefficients,
  --            computed from the bottom leaves to the top root;
  --   roco     the formal root count.

  procedure Resolve
              ( file : in file_type; ips : in out Intersection_Poset;
                flags : in out Standard_Complex_VecMats.VecMat;
                sols : out Solution_List );

  -- DESCRIPTION :
  --   Performs a bottom up root count, in the same order as the
  --   resolution with Littlewood-Richardson homotopies proceeds.
  --   This root counter can be viewed as a stub for the resolution
  --   of the intersection conditions by Littlewood-Richardson homotopies.

  -- REQUIRED :
  --   The intersection conditions are processed into an intersection poset,
  --   given in ips on entry.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   ips      an intersection poset built to resolve Schubert conditions;
  --   flags    generic complex matrices that represented nested linear
  --            space for use in the homotopies.

  -- ON RETURN :
  --   ips      intersection poset with Littlewood-Richardson coefficients,
  --            computed from the bottom leaves to the top root;
  --   flags    to glue the intersection conditions, the linear spaces
  --            have been multiplied with invertible matrices;
  --   sols     solutions to the Schubert problem, the length of
  --            this list must equal the formal root count.
 
end Resolve_Schubert_Problems;
