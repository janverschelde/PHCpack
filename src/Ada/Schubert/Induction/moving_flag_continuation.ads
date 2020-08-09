with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Natural_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Matrices;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Moving_Flag_Continuation is

-- DESCRIPTION :
--   One checker game defines a flag moving from special to general position.
--   Following the moving flag, path trackers compute solutions to a fixed
--   general flag and a special flag.

  procedure Track_First_Move
              ( file : in file_type; n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out Standard_Complex_Solutions.Link_to_Solution;
                fail : out boolean; rpt : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure Track_First_Move
              ( file : in file_type; n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 );
  procedure Track_First_Move
              ( file : in file_type; n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given a homotopy with last variable (with index n+1) the 
  --   continuation parameter, the start solution is computed and
  --   if that did not fail, the path tracker is launched,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics;
  --   n        number of variables in the ambient space;
  --   h        homotopy in n+1 variables;
  --   tosqr    true if the overdetermined homotopy will be squared up, or
  --            false if the Gauss-Newton path trackers will be applied;
  --   tol      tolerance on the residual to decide failure;
  --   sol      if not null, then the start solution for the homotopy,
  --            otherwise, the start solution will be computed;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sol      solution at the end of the path if not fail;
  --   fail     true if no start solution could be computed,
  --            or if the path tracker failed to reach a solution.

  procedure Track_Next_Move
              ( file : in file_type; nv : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out Standard_Complex_Solutions.Link_to_Solution;
                fail : out boolean; rpt : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure Track_Next_Move
              ( file : in file_type; nv : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 );
  procedure Track_Next_Move
              ( file : in file_type; nv : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks a path for the next move in the checker poset,
  --   given a homotopy with last variable (with index nv+1) the 
  --   continuation parameter and one start solution,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics;
  --   nv       number of variables in the ambient space;
  --   h        homotopy in nv+1 variables;
  --   tosqr    true if the overdetermined homotopy will be squared up, or
  --            false if the Gauss-Newton path trackers will be applied;
  --   tol      tolerance on the residual to decide failure;
  --   rpt      flag to run robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sol      solution at the end of the path if not fail;
  --   fail     true if there was no solution,
  --            or if the path tracker failed to reach a solution.

  procedure Track_Next_Move
              ( nv,nt : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; rpt : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure Track_Next_Move
              ( nv,nt : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 );
  procedure Track_Next_Move
              ( nv,nt : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks a path for the next move in the checker poset,
  --   given a homotopy with last variable (with index nv+1) the 
  --   continuation parameter and a start solution, without output,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   nv       number of variables in the ambient space;
  --   nt       number of tasks, if zero, then no multitasking;
  --   h        homotopy in nv+1 variables;
  --   tosqr    true if the overdetermined homotopy will be squared up, or
  --            false if the Gauss-Newton path trackers will be applied;
  --   tol      tolerance on the residual to decide failure;
  --   sols     start solutions for the homotopy;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sols     solutions at the end of the path if not fail;
  --   fail     true if there was no solution,
  --            or if the path tracker failed to reach a solution.

  procedure Track_Next_Move
              ( file : in file_type; nv,nt : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; rpt : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure Track_Next_Move
              ( file : in file_type; nv,nt : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 );
  procedure Track_Next_Move
              ( file : in file_type; nv,nt : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks a path for the next move in the checker poset,
  --   given a homotopy with last variable (with index n+1) the 
  --   continuation parameter and a start solution, with output to file,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics,
  --            if omitted, then no output is written to file;
  --   nv       number of variables in the ambient space;
  --   nt       number of tasks, if zero, then no multitasking;
  --   h        homotopy in nv+1 variables;
  --   tosqr    true if the overdetermined homotopy will be squared up, or
  --            false if the Gauss-Newton path trackers will be applied;
  --   tol      tolerance on the residual to decide failure;
  --   sols     start solutions for the homotopy;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sols     solutions at the end of the path if not fail;
  --   fail     true if there was no solution,
  --            or if the path tracker failed to reach a solution.

  procedure Initialize_Symbol_Table
              ( n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                dim : out integer32 );

  -- DESCRIPTION :
  --   Uses the localization pattern to initialize the symbol table.
  --   On return in dim is the number of free variables.

  procedure Generalizing_Homotopy
              ( n,k : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,nf : in Standard_Complex_Matrices.Matrix;
                h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                dim : out integer32 );
  procedure Generalizing_Homotopy
              ( file : in file_type; n,k : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,nf : in Standard_Complex_Matrices.Matrix;
                h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                dim : out integer32 );

  -- DESCRIPTION :
  --   A generalizing homotopy to move the black checkers from p to q.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics,
  --            if omitted, then there is no output.
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
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                x : in Standard_Complex_Vectors.Vector );
  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                x : in DoblDobl_Complex_Vectors.Vector );
  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                x : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Verifies the intersection conditions imposed by the input flags vf
  --   on a k-plane in n-space with attidutes of the intersection in cond,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     for intermedidate output and diagnostics;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   q        parent permutation in the poset, used for the localization;
  --   rows     position of the rows of the white checkers;
  --   cols     position of the columns of the white checkers;
  --   minrep   use a more efficient representation for the Schubert problem;
  --   cond     specifications for the intersection conditions;
  --   mf       coordinates for the moving flag;
  --   vf       coordinates for the input flags;
  --   x        current solution plane.

  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean );
  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean );
  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean );

  -- DESCRIPTION :
  --   This verification applies to a list of solutions.
  --   Verifies the intersection conditions imposed by the input flags vf
  --   on k-planes in n-space with attidutes of the intersection in cond,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     for intermedidate output and diagnostics;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   q        parent permutation in the poset, used for the localization;
  --   rows     position of the rows of the white checkers;
  --   cols     position of the columns of the white checkers;
  --   minrep   if true, then use the more efficient representation
  --            for the Schubert problems;
  --   cond     specifications for the intersection conditions;
  --   mf       coordinates for the moving flag;
  --   vf       coordinates for the input flags;
  --   sols     current solution planes;
  --   tol      tolerance on the residual to decide failure.

  -- ON RETURN :
  --   fail     true if all solutions evaluate to a residual that is
  --            higher than the fixed threshold of tol.

  procedure Copy_Flags ( src : in Standard_Complex_VecMats.VecMat;
                         dst : in out Standard_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   Copies the flags from source src to destination dst.

  -- REQUIRED : dst'range = src'range.

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                ls : in out Standard_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 );
  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 );
  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   In the trivial stay case instead of a homotopy,
  --   only a coordinate transformation is needed.  Computations happen
  --   in standard double, double double, or quad double precision.
  --   This version of Trivial_Stays deals with only one solution.

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
  --   verify   if verification is needed after coordinate transformation;
  --   minrep   to use the more efficient representation for the problem;
  --   cond     intersection conditions for the general fixed flags;
  --   mf       coordinates of the moving flag;
  --   vf       coordinates of general flags to keep fixed;
  --   ls       current solution if ind > 0;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ls       solution in the proper coordinates;
  --   fail     true if no longer a solution, false otherwise.

  procedure Trivial_Stay
              ( n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 );
  procedure Trivial_Stay
              ( n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 );
  procedure Trivial_Stay
              ( n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   In the trivial stay case instead of a homotopy,
  --   only a coordinate transformation is needed.  Computations happen
  --   in standard double, double double, or quad double precision.
  --   This is the silent version on a list of solutions.

  -- ON ENTRY :
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ctr      index of the critical row;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   minrep   to use a more efficient representation of the problem;
  --   sols     current list of solutions;
  --   tol      tolerance on the residual to decide failure;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sols     solutions in the proper coordinates;
  --   fail     true if no longer a solution, false otherwise.

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   In the trivial stay case instead of a homotopy,
  --   only a coordinate transformation is needed.  Computations happen
  --   in standard double, double double, or quad double precision.
  --   This version writes diagnostics to file.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics,
  --            if omitted, then there is no output and no checks either;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ctr      index of the critical row;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   verify   flag to indicate if diagnostic verification is needed;
  --   minrep   to use a more efficient representation of the problem;
  --   cond     intersection conditions for the general fixed flags;
  --   mf       coordinates of the moving flag;
  --   vf       coordinates of general flags to keep fixed;
  --   sols     current list of solutions;
  --   tol      tolerance on the residual to decide failure;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sols     solutions in the proper coordinates;
  --   fail     true if no longer a solution, false otherwise.

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                ls : in out Standard_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   In a stay homotopy, the white checkers stay in position.
  --   After tracking a path, the intersection conditions are verified,
  --   in standard double, double double, or quad double precision.
  --   This version is for one single solution.

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
  --   verify   flag to indicate if diagnostic verification is needed;
  --   minrep   to use a more efficient representation for the problem;
  --   tosqr    flag to indicate squaring of overdetermined homotopies;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   mf       the new moving flag at the target;
  --   start_mf is the moving flag at the start of the homotopy;
  --   ls       link to the current solution;
  --   tol      tolerance on the residual to decide failure;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ls       solution in the proper coordinates;
  --   fail     true if no longer a solution, false otherwise.

  procedure Stay_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Stay_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
  procedure Stay_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   In a stay homotopy, the white checkers stay in position.
  --   After tracking a path, the intersection conditions are verified,
  --   in standard double, double double, or quad double precision.
  --   This version is silent for a list of solutions.

  -- ON ENTRY :
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ctr      index of the critical row;
  --   nt       number of tasks, if zero, then no multitasking;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   minrep   to use a more efficient representation for the problem,
  --            only for checking with intermediate output;
  --   tosqr    flag to square overdetermined homotopies;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   mf       the new moving flag at the target;
  --   start_mf is the moving flag at the start of the homotopy;
  --   sols     the start solutions for the homotopy;
  --   tol      tolerance on the residual to decide failure;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sols     solution in the proper coordinates;
  --   fail     true if sols contains no longer a solution, false otherwise.

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   In a stay homotopy, the white checkers stay in position.
  --   After tracking a path, the intersection conditions are verified,
  --   in standard double, double double, or quad double precision.
  --   Diagnostics are written to file, for a list of solutions.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ctr      index of the critical row;
  --   nt       number of tasks, if zero, then no multitasking;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   verify   flag to indicate if diagnostic verification is needed;
  --   minrep   to use a more efficient representation for the problem,
  --            only for checking with intermediate output;
  --   tosqr    flag to indicate squaring of overdetermined homotopies;
  --   cond     intersection conditions for the general fixed flags;
  --   vf       coordinates of general flags to keep fixed;
  --   mf       the new moving flag at the target;
  --   start_mf is the moving flag at the start of the homotopy;
  --   sols     the start solutions for the homotopy;
  --   tol      tolerance on the residual to decide failure;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sols     solution in the proper coordinates;
  --   fail     true if sols contains no longer a solution, false otherwise.

  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                ls : in out Standard_Complex_Solutions.Link_to_Solution;
                rlq : out Standard_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 );
  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                rlq : out DoblDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 );
  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                rlq : out QuadDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   The variable k-plane x may have no pivot in column s in case of
  --   a swap homotopy which may cause that column to slide to zero.
  --   This procedure reconditions the swap homotopy for one start solution.
  --   Output is written to file.

  -- ON ENTRY :
  --   file     must be opened for writing output;
  --   dim      number of free variables in the localization map, locmap;
  --   r        index of the critical column;
  --   s        the swap column;
  --   locmap   localization pattern for x;
  --   x        variable coordinates for the solution k-plane;
  --   ls       start solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   x        matrix x with one scaling variable added if pividx /= 0;
  --   ls       reconditioned start solution if pividx /= 0;
  --   rlq      moving linear equation added for the extra scaling variable,
  --            is not null if pividx /= 0;
  --   pividx   index of the pivot in the swap column, which equals the
  --            index of x(r+1,s+1); if zero, then no recondition happened.

  procedure Recondition_Swap_Homotopy
              ( dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                rlq : out Standard_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 );
  procedure Recondition_Swap_Homotopy
              ( dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                rlq : out DoblDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 );
  procedure Recondition_Swap_Homotopy
              ( dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                rlq : out QuadDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   The variable k-plane x may have no pivot in column s in case of
  --   a swap homotopy which may cause that column to slide to zero.
  --   This procedure reconditions the swap homotopy for many start solutions.
  --   No output is written.

  -- ON ENTRY :
  --   dim      number of free variables in the localization map, locmap;
  --   r        index of the critical column;
  --   s        the swap column;
  --   locmap   localization pattern for x;
  --   x        variable coordinates for the solution k-plane;
  --   sols     start solutions;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   x        matrix x with one scaling variable added if pividx /= 0;
  --   ls       reconditioned start solution if pividx /= 0;
  --   rlq      moving linear equation added for the extra scaling variable,
  --            is not null if pividx /= 0;
  --   pividx   index of the pivot in the swap column, which equals the
  --            index of x(r+1,s+1); if zero, then no recondition happened.

  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                rlq : out Standard_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 );
  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                rlq : out DoblDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 );
  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                rlq : out QuadDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   The variable k-plane x may have no pivot in column s in case of
  --   a swap homotopy which may cause that column to slide to zero.
  --   This procedure reconditions the swap homotopy for many start solutions.
  --   Output is written to file.

  -- ON ENTRY :
  --   file     must be opened for writing output;
  --   dim      number of free variables in the localization map, locmap;
  --   r        index of the critical column;
  --   s        the swap column;
  --   locmap   localization pattern for x;
  --   x        variable coordinates for the solution k-plane;
  --   sols     start solutions;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   x        matrix x with one scaling variable added if pividx /= 0;
  --   ls       reconditioned start solution if pividx /= 0;
  --   rlq      moving linear equation added for the extra scaling variable,
  --            is not null if pividx /= 0;
  --   pividx   index of the pivot in the swap column, which equals the
  --            index of x(r+1,s+1); if zero, then no recondition happened.

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                ls : in out Standard_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   In a swap homotopy, two white checkers get swapped.
  --   After tracking a path, the intersection conditions are verified,
  --   in standard double, double double, or quad double precision.
  --   This version is for one solution.

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
  --   verify   flag to indicate if diagnostic verification is needed;
  --   minrep   to use a more efficient representation for the problem;
  --   tosqr    flag to indicate squaring of overdetermined homotopies;
  --   cond     intersection conditions for the general fixed flags;
  --   mf       coordinates of the moving flag;
  --   start_mf is the moving flag at the start of the homotopy;
  --   vf       coordinates of general flags to keep fixed;
  --   ls       link to the current solution;
  --   tol      tolerance on the residual to decide failure;
  --   rpt      flag to run robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   ls       solution in the proper coordinates;
  --   fail     true if no longer a solution, false otherwise.

  procedure Swap_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Swap_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
  procedure Swap_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );


  -- DESCRIPTION :
  --   In a swap homotopy, two white checkers get swapped.
  --   After tracking a path, the intersection conditions are verified,
  --   in standard double, double double, or quad double precision.
  --   This silent version works on a list of solutions.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics,
  --            if omitted, then there is no intermediate output;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ctr      index of the critical row;
  --   nt       number of tasks, if zero, then no multitasking;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   minrep   to use a more efficient representation for the problem;
  --   tosqr    to square overdetermined homotopies;
  --   cond     intersection conditions for the general fixed flags;
  --   mf       coordinates of the moving flag;
  --   start_mf is the moving flag at the start of the homotopy;
  --   vf       coordinates of general flags to keep fixed;
  --   sols     start solutions for the homotopy;
  --   tol      tolerance on the residual to decide failure;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sols     solutions in the proper coordinates;
  --   fail     true if sols contains no longer a solution, false otherwise.

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   In a swap homotopy, two white checkers get swapped.
  --   After tracking a path, the intersection conditions are verified,
  --   in standard double, double double, or quad double precision.
  --   Diagnostics are written on file, for a list of solutions.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics,
  --            if omitted, then there is no intermediate output;
  --   n        dimension of the ambient space, number of black checkers;
  --   k        dimension of the plane, number of white checkers;
  --   ctr      index of the critical row;
  --   ind      index of the move, if 0 then ls will be a first solution;
  --   nt       number of tasks, if zero, then no multitasking;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   verify   flag to indicate if diagnostic verification is needed;
  --   minrep   to use a more efficient representation for the problem;
  --   tosqr    to square the overdetermined homotopies;
  --   cond     intersection conditions for the general fixed flags;
  --   mf       coordinates of the moving flag;
  --   start_mf is the moving flag at the start of the homotopy;
  --   vf       coordinates of general flags to keep fixed;
  --   sols     start solutions for the homotopy;
  --   tol      tolerance on the residual to decide failure;
  --   rpt      flag to run the robust path tracker;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   sols     solutions in the proper coordinates;
  --   fail     true if sols contains no longer a solution, false otherwise.

end Moving_Flag_Continuation;
