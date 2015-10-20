with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Complex_Solutions;
with Checker_Posets;                      use Checker_Posets;
with Intersection_Solution_Posets;        use Intersection_Solution_Posets;

package Checker_Poset_Deformations is

-- DESCRIPTION :
--   The procedures in this package track paths defined by checker games.

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in out Standard_Complex_Matrices.Matrix;
                ls : in out Standard_Complex_Solutions.Link_to_Solution;
                tol : in double_float; unhappy : out boolean );

  -- DESCRIPTION :
  --   Tracks one path in the poset, given as an array of nodes.
  --   Start solutions are computed.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ps       checker poset for one game;
  --   path     path of nodes in the poset;
  --   count    number of the path;
  --   minrep   to use a more efficient problem formulation;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   mf       coordinates of the moving flag,
  --            should be equal to the identity matrix at the start;
  --   tol      tolerance on the residual to decide failure;

  -- ON RETURN :
  --   mf       moving flag at the end of the path;
  --   ls       solution at the end of the path;
  --   unhappy  true if the configuration of checkers is unhappy and
  --            gives no solution, also true if tolerance is not met.

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in out Standard_Complex_Matrices.Matrix;
                start : in Standard_Complex_Solutions.Solution_List;
                sols : out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; unhappy : out boolean );

  -- DESCRIPTION :
  --   Tracks one path in the poset, given as an array of nodes.
  --   Start solutions are provided in the nodes.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ps       checker poset for one game;
  --   path     path of nodes in the poset;
  --   count    number of the path;
  --   minrep   to use a more efficient problem formulation;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   mf       coordinates of the moving flag,
  --            should be equal to the identity matrix at the start;
  --   start    the start solutions from the nodes at the previous level
  --            in the intersection poset;
  --   tol      tolerance on residuals to decide failure.

  -- ON RETURN :
  --   mf       moving flag at the end of the path;
  --   sols     solutions at the end of the path in the poset;
  --   unhappy  true if the configuration of checkers is unhappy
  --            and gives no solution, true also if tolerance is not met.

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                tol : in double_float;
                sols : out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Tracks paths for one entire checker game in n-space,
  --   computing all start solutions.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ps       checker poset for one game;
  --   minrep   to use a more efficient problem formulation;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   tol      tolerance on residuals to decide failure.

  -- ON RETURN :
  --   sols     all solutions at the end of the paths.

  procedure Track_All_Paths_in_Poset
              ( n,k : in integer32; ps : in Poset;
                child : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                tol : in double_float;
                start : in Standard_Complex_Solutions.Solution_List;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                child : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                tol : in double_float;
                start : in Standard_Complex_Solutions.Solution_List;
                sols : out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Tracks paths for one entire checker game in n-space,
  --   with start solutions provided in the nodes,
  --   but only those paths that start at the matching child condition.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics,
  --            if omitted, then there is no intermediate output;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ps       checker poset for one game;
  --   child    conditions on the child for which the start solutions
  --            are provided and which should match leaves of ps;
  --   minrep   to use a more efficient problem formulation;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   start    contains solutions of the previous level, transformed
  --            to serve as the start solutions for the current level;
  --   tol      tolerance on residuals to decide failure.

  -- ON RETURN :
  --   sols     all solutions at the end of the paths.

end Checker_Poset_Deformations;
