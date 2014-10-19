with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Complex_Poly_Systems;       use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;          use Standard_Complex_Solutions;
with Checker_Posets;                      use Checker_Posets;
with Intersection_Solution_Posets;        use Intersection_Solution_Posets;

package Moving_Flag_Continuation is

-- DESCRIPTION :
--   One checker game defines a flag moving from special to general position.
--   Following the moving flag, path trackers compute solutions to a fixed
--   general flag and a special flag.

  procedure Set_Parameters ( file : in file_type; report : out boolean );

  -- DESCRIPTION :
  --   Interactive determination of the continuation and output parameters.

  procedure Call_Path_Trackers
              ( n : in integer32; h : in Poly_Sys;
                xt : in out Standard_Complex_Vectors.Vector;
                sol : out Link_to_Solution ); 
  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32; h : in Poly_Sys;
                xt : in out Standard_Complex_Vectors.Vector;
                sol : out Link_to_Solution ); 

  -- DESCRIPTION :
  --   Tracks one path starting at the solution in xt,
  --   as defined by the homotopy h.

  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics,
  --            if omitted, then there is no intermediate output;
  --   n        number of variables in the ambient space;
  --   h        homotopy in n+1 variables;
  --   xt       start solution with its last component equal to zero,
  --            satisfies the homotopy h (upto tolerance).

  -- ON RETURN :
  --   xt       solution at the end of the path, tracked to the
  --            last component of xt to be equal to one;
  --   sol      standard representation of the solution.

  procedure Track_First_Move
              ( file : in file_type; n : in integer32; h : in Poly_Sys;
                sol : in out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   Given a homotopy with last variable (with index n+1) the 
  --   continuation parameter, the start solution is computed and
  --   if that did not fail, a path tracker is launched.

  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics;
  --   n        number of variables in the ambient space;
  --   homtp    type of homotopy, if 0, then only start solution
  --            is computed and verified;
  --   h        homotopy in n+1 variables;
  --   sol      if not null, then the start solution for the homotopy,
  --            otherwise, the start solution will be computed.

  -- ON RETURN :
  --   sol      solution at the end of the path if not fail;
  --   fail     true if no start solution could be computed,
  --            or if the path tracker failed to reach a solution.

  procedure Track_Next_Move
              ( file : in file_type; n : in integer32; h : in Poly_Sys;
                sol : in out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   Tracks a path for the next move in the checker poset,
  --   given a homotopy with last variable (with index n+1) the 
  --   continuation parameter and a start solution.

  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics;
  --   n        number of variables in the ambient space;
  --   h        homotopy in n+1 variables.

  -- ON RETURN :
  --   sol      solution at the end of the path if not fail;
  --   fail     true if there was no solution,
  --            or if the path tracker failed to reach a solution.

  procedure Generalizing_Homotopy
              ( file : in file_type; n,k : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,nf : in Standard_Complex_Matrices.Matrix;
                h : out Link_to_Poly_Sys; dim : out integer32 );

  -- DESCRIPTION :
  --   A generalizing homotopy to move the black checkers from p to q.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   q        parent permutation in the specializing poset;
  --   p        current permutation determines shape of moving flag;
  --   rows     row positions of the white checkers at root of poset;
  --   cols     column position of the white checkers at root of poset;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   mf       coordinates of the moving flag;
  --   nf       current accumulated coordinates of the moving flag.

  -- ON RETURN :
  --   h        polynomial equations in the generalizing homotopy;
  --   dim      number of variables (excluding the continuation parameter)
  --            in the homotopy h.

  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                x : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Verifies the intersection conditions imposed by the input flags vf
  --   on a k-plane in n-space with attidutes of the intersection in cond.

  -- ON ENTRY :
  --   file     for intermedidate output and diagnostics;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   q        parent permutation in the poset, used for the localization;
  --   rows     position of the rows of the white checkers;
  --   cols     position of the columns of the white checkers;
  --   cond     specifications for the intersection conditions;
  --   mf       coordinates for the moving flag;
  --   vf       coordinates for the input flags;
  --   x        current solution plane.

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                ls : in out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   In the trivial stay case instead of a homotopy,
  --   only a coordinate transformation is needed.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ctr      index of the critical row;
  --   ind      index of the move, if 0 then ls will be a first solution;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   cond     intersection conditions for the general fixed flags;
  --   mf       coordinates of the moving flag;
  --   vf       coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   ls       solution in the proper coordinates;
  --   fail     true if no longer a solution, false otherwise.

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                ls : in out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   In a stay homotopy, the white checkers stay in position.
  --   After tracking a path, the intersection conditions are verified.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ctr      index of the critical row;
  --   ind      index of the move, if 0 then ls will be a first solution;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   mf       the new moving flag at the target;
  --   start_mf is the moving flag at the start of the homotopy.

  -- ON RETURN :
  --   ls       solution in the proper coordinates;
  --   fail     true if no longer a solution, false otherwise.

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                ls : in out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   In a swap homotopy, two white checkers get swapped.
  --   After tracking a path, the intersection conditions are verified.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ctr      index of the critical row;
  --   ind      index of the move, if 0 then ls will be a first solution;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   cond     intersection conditions for the general fixed flags;
   --  mf       coordinates of the moving flag;
  --   start_mf is the moving flag at the start of the homotopy;
  --   vf       coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   ls       solution in the proper coordinates;
  --   fail     true if no longer a solution, false otherwise.

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in out Standard_Complex_Matrices.Matrix;
                ls : in out Link_to_Solution; unhappy : out boolean );

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
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   mf       coordinates of the moving flag,
  --            should be equal to the identity matrix at the start.

  -- ON RETURN :
  --   mf       moving flag at the end of the path;
  --   ls       solution at the end of the path;
  --   unhappy  true if the configuration of checkers is unhappy
  --            and gives no solution.

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in out Standard_Complex_Matrices.Matrix;
                snd : in Link_to_Solution_Node; ls : in out Link_to_Solution;
                unhappy : out boolean );

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
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   mf       coordinates of the moving flag,
  --            should be equal to the identity matrix at the start;
  --   snd      start solutions from the nodes at the previous level
  --            in the intersection poset.

  -- ON RETURN :
  --   mf       moving flag at the end of the path;
  --   ls       solution at the end of the path;
  --   unhappy  true if the configuration of checkers is unhappy
  --            and gives no solution.

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : out Solution_List );

  -- DESCRIPTION :
  --   Tracks paths for one entire checker game in n-space,
  --   computing all start solutions.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ps       checker poset for one game;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   sols     all solutions at the end of the paths.

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                snd : in Link_to_Solution_Node; sols : out Solution_List );

  -- DESCRIPTION :
  --   Tracks paths for one entire checker game in n-space,
  --   with start solutions provided in the nodes.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ps       checker poset for one game;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   snd      solution nodes at the previous level with the start
  --            solutions for the current level.

  -- ON RETURN :
  --   sols     all solutions at the end of the paths.

end Moving_Flag_Continuation;
