with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Solutions;
with Standard_Solution_Posets;
with DoblDobl_Complex_Solutions;
with DoblDobl_Solution_Posets;
with QuadDobl_Complex_Solutions;
with QuadDobl_Solution_Posets;
with Intersection_Posets;                use Intersection_Posets;

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

  procedure Initialize_Solution_Nodes
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                nodes : in out Standard_Solution_Posets.Solnode_List;
                res : out double_float );
  procedure Initialize_Solution_Nodes
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                nodes : in out DoblDobl_Solution_Posets.Solnode_List;
                res : out double_double );
  procedure Initialize_Solution_Nodes
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                nodes : in out QuadDobl_Solution_Posets.Solnode_List;
                res : out quad_double );

  -- DESCRIPTION :
  --   Initializes the nodes with the start solutions, computed
  --   in standard double, double double, or quad double precision,
  --   with diagnostic output written to file.

  -- ON ENTRY :
  --   file     for diagnostic output;
  --   n        ambient dimension;
  --   k        dimension of the solution planes;
  --   conds    intersection conditions;
  --   flags    generic flags in n-space;
  --   nodes    list of solution nodes.

  -- ON RETURN :
  --   nodes    initialized with start solutions;
  --   res      sum of all residuals.

  procedure Initialize_Solution_Nodes
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                nodes : in out Standard_Solution_Posets.Solnode_List;
                res : out double_float );
  procedure Initialize_Solution_Nodes
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                nodes : in out DoblDobl_Solution_Posets.Solnode_List;
                res : out double_double );
  procedure Initialize_Solution_Nodes
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                nodes : in out QuadDobl_Solution_Posets.Solnode_List;
                res : out quad_double );

  -- DESCRIPTION :
  --   Initializes the nodes with the start solutions, computed
  --   in standard double, double double, or quad double precision,
  --   without diagnostic output written to file.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution planes;
  --   conds    intersection conditions;
  --   flags    generic flags in n-space;
  --   nodes    list of solution nodes.

  -- ON RETURN :
  --   nodes    initialized with start solutions;
  --   res      sum of all residuals.

  procedure Start_Solution 
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                snd : in Standard_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out double_float );
  procedure Start_Solution 
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                snd : in DoblDobl_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out double_double );
  procedure Start_Solution 
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                snd : in QuadDobl_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out quad_double );

  -- DESCRIPTION :
  --   Computes the start solution at a solution node,
  --   in standard double, double double, or quad double precision,
  --   positioned at the leaves of the intersection poset,
  --   as the analogue to the Initialize_Leaves from above.

  -- REQUIRED :
  --   The solution node comes with a valid checker poset.
  --   Assumed is that there conditions turn into a linear system
  --   and there is at most only one solution.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        the ambient dimension;
  --   k        dimension of the solution plane;
  --   conds    conditions on the solution planes meeting the flags in vf;
  --   flags    fixed flags for setting up the intersection conditions;
  --   snd      a solution node.

  -- ON RETURN :
  --   snd      if not fail, the solution node contains a start solution;
  --   fail     true if the residual is higher than the 1.0e-8 threshold;
  --   x        the computed solution vector;
  --   res      the residual as the two norm of the solution evaluated
  --            at the polynomial equations that express the intersection
  --            conditions imposed by the brackets in the poset,
  --            the brackets in cond and the fixed flags.

  procedure Start_Solution 
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                snd : in Standard_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out double_float );
  procedure Start_Solution 
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                snd : in DoblDobl_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out double_double );
  procedure Start_Solution 
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                snd : in QuadDobl_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out quad_double );

  -- DESCRIPTION :
  --   Computes the start solution at a solution node,
  --   in standard double, double double, or quad double precision,
  --   positioned at the leaves of the intersection poset,
  --   as the analogue to the Initialize_Leaves from above.
  --   These versions are silent: they produce no output.

  -- REQUIRED :
  --   The solution node comes with a valid checker poset.
  --   Assumed is that there conditions turn into a linear system
  --   and there is at most only one solution.

  -- ON ENTRY :
  --   n        the ambient dimension;
  --   k        dimension of the solution plane;
  --   conds    conditions on the solution planes meeting the flags in vf;
  --   flags    fixed flags for setting up the intersection conditions;
  --   snd      a solution node.

  -- ON RETURN :
  --   snd      if not fail, the solution node contains a start solution;
  --   fail     true if the residual is higher than the 1.0e-8 threshold;
  --   x        the computed solution vector;
  --   res      the residual as the two norm of the solution evaluated
  --            at the polynomial equations that express the intersection
  --            conditions imposed by the brackets in the poset,
  --            the brackets in cond and the fixed flags.

  procedure Transform_Start_Solutions
              ( n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Transform_Start_Solutions
              ( n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Transform_Start_Solutions
              ( n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Applies the tm transformation to all solutions and column reduces
  --   the result to fit in the pattern defined by rows and cols,
  --   in standard double, double double, or quad double precision,
  --   without intermediate output to file.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane
  --   r_src    rows of the conditions at the source: incoming solutions;
  --   c_src    columns of the conditions at the source;
  --   r_src    rows of the conditions at the target: transformed solutions;
  --   c_src    columns of the conditions at the target;
  --   tm       transformation matrix, equals T1*M;
  --   sols     solutions that fit the pattern as defined by r_src, c_src.

  -- ON RETURN :
  --   sols     transformed solutions that fit the pattern prescribed by
  --            r_tgt and c_tgt.

  procedure Transform_Start_Solutions
              ( file : in file_type; n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Transform_Start_Solutions
              ( file : in file_type; n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Transform_Start_Solutions
              ( file : in file_type; n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Applies the tm transformation to all solutions and column reduces
  --   the result to fit in the pattern defined by rows and cols,
  --   in standard double, double double, or quad double precision,
  --   with diagnostic intermediate output written to file.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics,
  --            if omitted, then there is no intermediate output;
  --   n        ambient dimension;
  --   k        dimension of the solution plane
  --   r_src    rows of the conditions at the source: incoming solutions;
  --   c_src    columns of the conditions at the source;
  --   r_src    rows of the conditions at the target: transformed solutions;
  --   c_src    columns of the conditions at the target;
  --   tm       transformation matrix, equals T1*M;
  --   sols     solutions that fit the pattern as defined by r_src, c_src.

  -- ON RETURN :
  --   sols     transformed solutions that fit the pattern prescribed by
  --            r_tgt and c_tgt.

  procedure Connect_Checker_Posets_to_Count
              ( file : in file_type;
                pl : in Poset_List; nd : in Poset_Node;
                vrblvl : in integer32 := 0 );

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

  procedure Connect_Checker_Posets_to_Track
              ( n,k,level,nt : in integer32; tol : in double_float;
                pl : in Poset_List;
                snd : in Standard_Solution_Posets.Link_to_Solution_Node;
                tmfo : in Standard_Complex_Matrices.Link_to_Matrix;
                sps : in out Standard_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Connect_Checker_Posets_to_Track
              ( n,k,level,nt : in integer32; tol : in double_float;
                pl : in Poset_List;
                snd : in DoblDobl_Solution_Posets.Link_to_Solution_Node;
                tmfo : in DoblDobl_Complex_Matrices.Link_to_Matrix;
                sps : in out DoblDobl_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                vrblvl : in integer32 := 0 );
  procedure Connect_Checker_Posets_to_Track
              ( n,k,level,nt : in integer32; tol : in double_float;
                pl : in Poset_List;
                snd : in QuadDobl_Solution_Posets.Link_to_Solution_Node;
                tmfo : in QuadDobl_Complex_Matrices.Link_to_Matrix;
                sps : in out QuadDobl_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Extension of the Connect_Checker_Posets_to_Count to track the
  --   paths in standard double, double double, or quad double precision.
  --   This version is silent, does not write diagnostic output to file.

  -- ON ENTRY :
  --   n        the ambient dimension;
  --   k        dimension of the solution planes;
  --   level    level of the parent nodes in the list pl;
  --   nt       number of tasks, if zero, then no multitasking;
  --   pl       list of checker posets at some level of the parent nodes
  --            to the node nd in the intersection poset;
  --   snd      solution node that contains the poset of the child;
  --   tmfo     transformation for use at start solution if not null;
  --   sps      solution poset constructed up to the proper level;
  --   minrep   to use a more efficient problem formulation;
  --   tosqr    to square overdetermined homotopies;
  --   conds    conditions on the current fixed flags;
  --   flags    current fixed flags;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sps      updated solution poset.

  procedure Connect_Checker_Posets_to_Track
              ( file : in file_type;
                n,k,level,nt : in integer32; tol : in double_float;
                pl : in Poset_List;
                snd : in Standard_Solution_Posets.Link_to_Solution_Node;
                tmfo : in Standard_Complex_Matrices.Link_to_Matrix;
                sps : in out Standard_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Connect_Checker_Posets_to_Track
              ( file : in file_type;
                n,k,level,nt : in integer32; tol : in double_float;
                pl : in Poset_List;
                snd : in DoblDobl_Solution_Posets.Link_to_Solution_Node;
                tmfo : in DoblDobl_Complex_Matrices.Link_to_Matrix;
                sps : in out DoblDobl_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                vrblvl : in integer32 := 0 );
  procedure Connect_Checker_Posets_to_Track
              ( file : in file_type;
                n,k,level,nt : in integer32; tol : in double_float;
                pl : in Poset_List;
                snd : in QuadDobl_Solution_Posets.Link_to_Solution_Node;
                tmfo : in QuadDobl_Complex_Matrices.Link_to_Matrix;
                sps : in out QuadDobl_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Extension of the Connect_Checker_Posets_to_Count to track the
  --   paths in standard double, double double, or quad double precision.
  --   This version writes diagnostic output to file.

  -- ON ENTRY :
  --   file     for intermediate output, if provided, otherwise silent;
  --   n        the ambient dimension;
  --   k        dimension of the solution planes;
  --   level    level of the parent nodes in the list pl;
  --   nt       number of tasks, if zero, then no multitasking;
  --   pl       list of checker posets at some level of the parent nodes
  --            to the node nd in the intersection poset;
  --   snd      solution node that contains the poset of the child;
  --   tmfo     transformation for use at start solution if not null;
  --   sps      solution poset constructed up to the proper level;
  --   verify   flag to indicate whether diagnostic verification is needed;
  --   minrep   to use a more efficient problem formulation;
  --   tosqr    to square overdetermined homotopies;
  --   conds    conditions on the current fixed flags;
  --   flags    current fixed flags;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sps      updated solution poset.

  procedure Count_Roots
              ( ips : in out Intersection_Poset;
                roco : out Natural_Number; vrblvl : in integer32 := 0 );
  procedure Count_Roots
              ( file : in file_type; ips : in out Intersection_Poset;
                roco : out Natural_Number; vrblvl : in integer32 := 0 );

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
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ips      intersection poset with Littlewood-Richardson coefficients,
  --            computed from the bottom leaves to the top root;
  --   roco     the formal root count.

  procedure Resolve
              ( file : in file_type; extopt,repcon : in boolean;
                n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out Standard_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                sols : out Standard_Complex_Solutions.Solution_List;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Resolve
              ( file : in file_type; extopt,repcon : in boolean;
                n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out DoblDobl_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );
  procedure Resolve
              ( file : in file_type; extopt,repcon : in boolean;
                n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out QuadDobl_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies the Littlewood-Richardson homotopies running in
  --   standard double, double double, or quad double precision,
  --   to resolve a sequence of intersection conditions,
  --   with diagnostic output written to file.

  -- REQUIRED :
  --   The intersection conditions are processed into an intersection poset,
  --   given in ips on entry.  Moreover, conds'range = flags'range.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   extopt   extra output about monitoring the Littlewood-Richardson
  --            homotopies in each and every checker game;
  --   repcon   true for the path trackers to run in reporting version,
  --            false if the path trackers have to stay mute;
  --   n        the ambient dimension;
  --   k        dimension of the solution plane;
  --   nt       number of tasks, if zero, then no multitasking;
  --   tol      tolerance on residual to decide failure in checker games;
  --   ips      an intersection poset built to resolve Schubert conditions;
  --   sps      an initialized solution poset corresponding to ips;
  --   verify   flag to indicate whether diagnostic verification is needed;
  --   minrep   to use a more efficient problem formulation;
  --   tosqr    to square overdetermined homotopies;
  --   conds    intersection conditions on the fixed flags;
  --   flags    generic complex matrices that represented nested linear
  --            space for use in the homotopies;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ips      intersection poset with Littlewood-Richardson coefficients,
  --            computed from the bottom leaves to the top root;
  --   sps      solution poset with at each level the corresponding solutions;
  --   sols     solutions to the Schubert problem, the length of
  --            this list must equal the formal root count.

  procedure Resolve
              ( n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out Standard_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                sols : out Standard_Complex_Solutions.Solution_List;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Resolve
              ( n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out DoblDobl_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );
  procedure Resolve
              ( n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out QuadDobl_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 );
 
  -- DESCRIPTION :
  --   Applies the Littlewood-Richardson homotopies running in
  --   standard double, double double, or quad double precision,
  --   to resolve a sequence of intersection conditions,
  --   without writing any diagnostic output to file.

  -- REQUIRED :
  --   The intersection conditions are processed into an intersection poset,
  --   given in ips on entry.  Moreover, conds'range = flags'range.

  -- ON ENTRY :
  --   extopt   extra output about monitoring the Littlewood-Richardson
  --            homotopies in each and every checker game;
  --   repcon   true for the path trackers to run in reporting version,
  --            false if the path trackers have to stay mute;
  --   n        the ambient dimension;
  --   k        dimension of the solution plane;
  --   nt       number of tasks, if zero, then no multitasking;
  --   tol      tolerance on residual to decide failure in checker games;
  --   ips      an intersection poset built to resolve Schubert conditions;
  --   sps      an initialized solution poset corresponding to ips;
  --   minrep   to use a more efficient problem formulation;
  --   tosqr    to square overdetermined homotopies;
  --   conds    intersection conditions on the fixed flags;
  --   flags    generic complex matrices that represented nested linear
  --            space for use in the homotopies;
  --   rpt      flag to run the robust path tracker;
  --   verblvl  is the verbose level.

  -- ON RETURN :
  --   ips      intersection poset with Littlewood-Richardson coefficients,
  --            computed from the bottom leaves to the top root;
  --   sps      solution poset with at each level the corresponding solutions;
  --   sols     solutions to the Schubert problem, the length of
  --            this list must equal the formal root count.

end Resolve_Schubert_Problems;
