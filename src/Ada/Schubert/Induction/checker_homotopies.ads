with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Matrices;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Poly_Matrices;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Checker_Homotopies is

-- DESCRIPTION :
--   Captures the definitions of the Littlewood-Richardson homotopies,
--   as defined by the relative positions of white and black checkers.
--   This library consists in three parts:
--   (1) local classification of Littlewood-Richardson homotopies;
--   (2) coordinate transformations on input flags and solution plane;
--   (3) coordinate definitions for the two types of swap homotopies.

-- PART I : classification of Littlewood-Richardson homotopies

  procedure Define_Specializing_Homotopy
              ( n : in integer32;
                p,row,col : in Standard_Natural_Vectors.Vector;
                homtp,ctr : out integer32 );
  procedure Define_Specializing_Homotopy
              ( file : in file_type; n : in integer32;
                p,row,col : in Standard_Natural_Vectors.Vector;
                homtp,ctr : out integer32 );

  -- DESCRIPTION :
  --   To deform an k-plane in n-space, information about a homotopy
  --   is printed to screen, based on a checker configuration defined
  --   by a permutation in p and row and column indices.
  --   This routine is called when traversing the poset from the root
  --   to the leaves, so there is choice at stay/swap case.

  -- ON RETURN :
  --   homtp    type of homotopy, 0 if none, 1 if stay, 2 if swap;
  --   ctr      index of the critical row.

  procedure Define_Generalizing_Homotopy 
              ( n : in integer32;
                p,row,col : in Standard_Natural_Vectors.Vector;
                stay_child : in boolean; homtp,ctr : out integer32 );
  procedure Define_Generalizing_Homotopy 
              ( file : in file_type; n : in integer32;
                p,row,col : in Standard_Natural_Vectors.Vector;
                stay_child : in boolean; homtp,ctr : out integer32 );

  -- DESCRIPTION :
  --   To deform an k-plane in n-space, information about a homotopy
  --   is printed to screen, based on a checker configuration defined
  --   by a permutation in p and row and column indices.
  --   This routine is called when traversing the poset from the leaves,
  --   in which case type of homotopy depends on whether the child is
  --   a "stay child" or a "swap child", as indicated by the boolean
  --   variable stay_child.

  -- ON RETURN :
  --   homtp    type of homotopy, 0 if none, 1 if stay, 2 if swap;
  --   ctr      index of the critical row.

-- PART II : coordinate transformations on input flags and solution plane

  procedure Inverse_Row_Transformation
              ( r : in integer32;
                x : in out Standard_Complex_Matrices.Matrix );
  procedure Inverse_Row_Transformation
              ( r : in integer32;
                x : in out DoblDobl_Complex_Matrices.Matrix );
  procedure Inverse_Row_Transformation
              ( r : in integer32;
                x : in out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the inverse coordinate transformation to a solution plane
  --   in case no homotopies are needed.  The transformation happens in
  --   standard double, double double, or quad double precision.

  -- REQUIRED : r < x'last(1).

  -- ON ENTRY :
  --   r        the index of the critical row;
  --   x        matrix representation of a solution k-plane,
  --            as an n-by-k matrix.

  -- ON RETURN :
  --   x        the inverse coordinate transformation is applied to x.

  procedure Inverse_Row_Transformation
              ( mf : in Standard_Complex_Matrices.Matrix;
                x : in out Standard_Complex_Matrices.Matrix );
  procedure Inverse_Row_Transformation
              ( mf : in DoblDobl_Complex_Matrices.Matrix;
                x : in out DoblDobl_Complex_Matrices.Matrix );
  procedure Inverse_Row_Transformation
              ( mf : in QuadDobl_Complex_Matrices.Matrix;
                x : in out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Applies the inverse coordinate transformation to a solution plane
  --   in the case homotopies are needed.  The transformation happens in
  --   standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   mf       the new current moving flag;
  --   x        matrix representation of a solution k-plane,
  --            as an n-by-k matrix.

  -- ON RETURN :
  --   x        the inverse coordinate transformation is applied to x.

  procedure Normalize_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Matrices.Matrix );
  procedure Normalize_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Matrices.Matrix );
  procedure Normalize_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Normalizes and column reduces the solution in x to fit 
  --   the localization pattern, after application of the inverse
  --   row transformation on the solution plane x.  The normalization is
  --   done in standard double, double double, or quad double precision.

  -- REQUIRED : The rows where ones are expected according to the 
  --   localization pattern in  the matrix x must be different from zero.

  -- ON ENTRY :
  --   pattern  a matrix of 0, 1, and 2 entries represents a localization
  --            pattern for a k-plane: 2 for arbitrary numbers,
  --            and 0 and 1 at position where zero and ones are expected;
  --   x        matrix representation of a solution k-plane,
  --            as an n-by-k matrix, after Inverse_Row_Transformation(r,x),
  --            where r is the critical row.

  -- ON RETURN :
  --   x        normalized so ones appear at the expected places.

  procedure Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Matrices.Matrix );
  procedure Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Matrices.Matrix );
  procedure Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Reduces the columns of the matrix x according to the zeres in
  --   the localization pattern.  The reduction happens in
  --   standard double, double double, or quad double precision.

  -- REQUIRED :
  --   The matrix x has been normalized so the ones appear on the
  --   places prescribed by the localization pattern.

  -- ON ENTRY :
  --   pattern  a matrix of 0, 1, and 2 entries represents a localization
  --            pattern for a k-plane: 2 for arbitrary numbers,
  --            and 0 and 1 at position where zero and ones are expected;
  --   x        matrix representation of a solution k-plane,
  --            as an n-by-k matrix, after Inverse_Row_Transformation(r,x),
  --            where r is the critical row, and after Normalize_to_Fit
  --            has been applied to have ones in the proper places.

  -- ON RETURN :
  --   x        normalized so ones appear at the expected places,
  --            and column reduced so zeros are where expected.

  procedure Normalize_and_Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Matrices.Matrix );
  procedure Normalize_and_Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Matrices.Matrix );
  procedure Normalize_and_Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Normalizes and column reduces the solution in x to fit 
  --   the localization pattern, after application of the inverse
  --   coordinate row transformation on the solution plane x.
  --   The normalization and the reduction to fit a pattern happen
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   pattern  a matrix of 0, 1, and 2 entries represents a localization
  --            pattern for a k-plane: 2 for arbitrary numbers,
  --            and 0 and 1 at position where zero and ones are expected;
  --   x        matrix representation of a solution k-plane,
  --            as an n-by-k matrix, after Inverse_Row_Transformation(r,x),
  --            where r is the critical row.

  -- ON RETURN :
  --   x        normalized so ones appear at the expected places,
  --            and column reduced so zeros are where expected.

  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out Standard_Complex_Vectors.Vector );
  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out DoblDobl_Complex_Vectors.Vector );
  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs the change of variables on the solution x,
  --   in standard double, double double, or quad double precision,
  --   in case of a trivial stay case when no homotopy is needed.
  --   There is no intermediate output.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   x        current solution.

  -- ON RETURN :
  --   x        solution with transformed coordinates.

  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out Standard_Complex_Vectors.Vector );
  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out DoblDobl_Complex_Vectors.Vector );
  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs the change of variables on the solution x,
  --   in standard double, double double, or quad double precision,
  --   in case of a trivial stay case when no homotopy is needed.
  --   There is intermediate output to file during computations.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent version;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   x        current solution.

  -- ON RETURN :
  --   x        solution with transformed coordinates.

  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs the change of variables on the solutions in sols,
  --   in standard double, double double, or quad double precision,
  --   in case of a trivial stay case when no homotopy is needed.
  --   There is no intermediate output during computations.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   sols     list of current solutions.

  -- ON RETURN :
  --   sols     solutions with transformed coordinates.

  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs the change of variables on the solutions in sols,
  --   in standard double, double double, or quad double precision,
  --   in case of a trivial stay case when no homotopy is needed.
  --   There is intermediate output during computations.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent version;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   sols     list of current solutions.

  -- ON RETURN :
  --   sols     solutions with transformed coordinates.

  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector );
  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector );
  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs the change of variables on the solution x,
  --   in standard double, double double, or quad double precision,
  --   after a homotopy in the stay case, without intermediate output.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   p        permutation indicates location of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   x        current solution.

  -- ON RETURN :
  --   x        solution with transformed coordinates.

  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector );
  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector );
  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs the change of variables on the solution x,
  --   in standard double, double double, or quad double precision,
  --   after a homotopy in the stay case, with output to file.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent version;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   p        permutation indicates location of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   x        current solution.

  -- ON RETURN :
  --   x        solution with transformed coordinates.

  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs the change of variables on the solutions in the list sols,
  --   in standard double, double double, or quad double precision.
  --   after a homotopy in the stay case.  There is no output.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   p        permutation indicates location of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   sols     current solutions.

  -- ON RETURN :
  --   sols     list of solutions with transformed coordinates.

  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs the change of variables on the solutions in the list sols,
  --   in standard double, double double, or quad double precision,
  --   after a homotopy in the stay case, with output to file.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent version;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   p        permutation indicates location of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   sols     current solutions.

  -- ON RETURN :
  --   sols     list of solutions with transformed coordinates.

  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector );
  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector );
  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs the change of variables on the solution x,
  --   in standard double, double double, or quad double precision.
  --   after a swap homotopy of the first type.
  --   This version of the change of variables is silent.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   big_r    row of the swapped white checker, R > r + 1;
  --   dc       index of the descending black checker in specializing poset;
  --   s        columns s and s+1 in x are swapped;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   x        current solution.

  -- ON RETURN :
  --   x        solution with transformed coordinates.

  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector );
  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector );
  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs the change of variables on the solution x,
  --   in standard double, double double, or quad double precision,
  --   after a swap homotopy of the first type,
  --   with intermediate output written to file.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent version;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   big_r    row of the swapped white checker, R > r + 1;
  --   dc       index of the descending black checker in specializing poset;
  --   s        columns s and s+1 in x are swapped;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   x        current solution.

  -- ON RETURN :
  --   x        solution with transformed coordinates.

  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs the change of variables on the solutions in the list sols
  --   after a swap homotopy of the first type, without output,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   big_r    row of the swapped white checker, R > r + 1;
  --   dc       index of the descending black checker in specializing poset;
  --   s        columns s and s+1 in x are swapped;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   sols     current solutions.

  -- ON RETURN :
  --   sols     list of solutions with transformed coordinates.

  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs the change of variables on the solutions in the list sols,
  --   in standard double, double double, or quad double precision,
  --   after a swap homotopy of the first type, 
  --   writing intermediate output to file.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent version;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   big_r    row of the swapped white checker, R > r + 1;
  --   dc       index of the descending black checker in specializing poset;
  --   s        columns s and s+1 in x are swapped;
  --   q        parent permutation in the poset used for target;
  --   p        current permutation in the poset used for start;
  --   qr       position of the rows of the white checkers with q;
  --   qc       position of the columns of the white checkers with q;
  --   pr       position of the rows of the white checkers with p;
  --   pc       position of the columns of the white checkers with p;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   sols     current solutions.

  -- ON RETURN :
  --   sols     list of solutions with transformed coordinates.

  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector );
  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector );
  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs the change of variables on the solution x,
  --   in standard double, double double, or quad double precision,
  --   after a swap homotopy of the second type, without output.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   s        columns s and s+1 in x are swapped;
  --   p        permutation indicates location of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   x        current solution.

  -- ON RETURN :
  --   x        solution with transformed coordinates.

  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector );
  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector );
  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs the change of variables on the solution x,
  --   in standard double, double double, or quad double precision,
  --   after a swap homotopy of the second type,
  --   writing intermediate output to file.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent version;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   s        columns s and s+1 in x are swapped;
  --   p        permutation indicates location of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   x        current solution.

  -- ON RETURN :
  --   x        solution with transformed coordinates.

  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs the change of variables on the solutions in the list sols,
  --   in standard double, double double, or quad double precision,
  --   after a swap homotopy of the second type, without output.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   s        columns s and s+1 in x are swapped;
  --   p        permutation indicates location of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   sols     current solutions.

  -- ON RETURN :
  --   sols     list of solutions with transformed coordinates.

  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs the change of variables on the solutions in the list sols,
  --   in standard double, double double, or quad double precision,
  --   after a swap homotopy of the second type,
  --   with intermediate output written to file.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent version;
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   r        the critical row;
  --   s        columns s and s+1 in x are swapped;
  --   p        permutation indicates location of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers;
  --   mf       the coordinates of the new moving flag;
  --   xtm      localization pattern extended with t;
  --   sols     current solutions.

  -- ON RETURN :
  --   sols     list of solutions with transformed coordinates.

-- PART III : coordinate definitions for the stay and swap homotopies

  function Swap_Column
              ( r : integer32; m : Standard_Natural_Matrices.Matrix )
              return integer32;

  -- DESCRIPTION :
  --   Given in r the critical row and in m a column pattern for a k-plane,
  --   the index of the column for which a one occurs in row r is returned.
  --   If the number on return equals zero, then no swap is possible.

  function Swap_Column
              ( r : integer32; rows : Standard_Natural_Vectors.Vector )
              return integer32;

  -- DESCRIPTION :
  --   Given in r the critical row and in rows the row indices of the
  --   white checkers, the index of the column to be swapped is returned.
  --   If the number on return equals zero, then no swap is possible.

  function Swap_Checker
              ( q,rows,cols : Standard_Natural_Vectors.Vector )
              return integer32;

  -- DESCRIPTION :
  --  The big R identifies the position of the second white checker
  --  swapped with the first one.  If R = r+1, then we have type II
  --  of the swap homotopy, type I corresponds to R > r+1.
  --  This function returns this big R.

  procedure Initialize_Moving_Plane
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix );
  procedure Initialize_Moving_Plane
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix );
  procedure Initialize_Moving_Plane
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix );

  -- DESCRIPTION :
  --   Initializes the coordinates of the moving plane x, using the
  --   localization pattern in m, with x matrix with coefficients in
  --   standard double, double double, or quad double precision.

  -- REQUIRED :
  --   x'range(1) = 1..n = m'range(1) and x'range(2) = 1..k = m'range(2).

  -- ON ENTRY :
  --   x        matrix of n rows and k rows;
  --   m        column localization pattern.

  -- ON RETURN :
  --   x        initialized coordinates.

  procedure Initialize_Moving_Plane
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix; s : in integer32 );
  procedure Initialize_Moving_Plane
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix; s : in integer32 );
  procedure Initialize_Moving_Plane
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix; s : in integer32 );

  -- DESCRIPTION :
  --   Initializes the coordinates of the moving plane x, using the
  --   localization pattern in m, for all columns except s and s + 1,
  --   if s is provided as an input parameter.  The matrix x is in
  --   standard double, double double, or quad double precision.

  -- REQUIRED :
  --   x'range(1) = 1..n = m'range(1) and x'range(2) = 1..k = m'range(2).

  -- ON ENTRY :
  --   x        matrix of n rows and k rows;
  --   m        column localization pattern;
  --   s        columns s and s+1 will be swapped.

  -- ON RETURN :
  --   x        initialized coordinates except columns s and s+1.

  function Is_Zone_A_Empty
             ( locmap : Standard_Natural_Matrices.Matrix;
               p : Standard_Natural_Vectors.Vector;
               r,s,dc : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if zone A is empty, returns false otherwise.

  -- ON ENTRY :
  --   locmap   Checker_Localization_Patterns.Column_Pattern(p,rows,cols);
  --   p        permutation defines current position of black checkers;
  --   r        critical row, row of the descending black checker;
  --   s        columns s and s+1 in x are swapped;
  --   dc       index of the descending black checker in specializing poset.

  procedure First_Swap_Plane
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );
  procedure First_Swap_Plane
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );
  procedure First_Swap_Plane
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );

  -- DESCRIPTION :
  --   In the first type of swap homotopy r+1 < R = big_r, where
  --   x is represented as a matrix of polynomials with coefficients
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   x        matrix of n rows and k rows;
  --   r        critical row, row of the descending black checker;
  --   big_r    row of the swapped white checker, R > r + 1;
  --   dc       index of the descending black checker in specializing poset;
  --   s        columns s and s+1 in x are swapped;
  --   p        permutation defines current position of black checkers;
  --   locmap   Checker_Localization_Patterns.Column_Pattern(p,rows,cols).

  -- ON RETURN :
  --   x        coordinates for swap homotopy of first type.

  procedure First_Swap_Plane
              ( file : in file_type;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );
  procedure First_Swap_Plane
              ( file : in file_type;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );
  procedure First_Swap_Plane
              ( file : in file_type;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );

  -- DESCRIPTION :
  --   In the first type of swap homotopy r+1 < R = big_r, where
  --   x is represented as a matrix of polynomials with coefficients
  --   in standard double, double double, or quad double precision.
  --   This version writes intermediate output to file.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent;
  --   x        matrix of n rows and k rows;
  --   r        critical row, row of the descending black checker;
  --   big_r    row of the swapped white checker, R > r + 1;
  --   dc       index of the descending black checker in specializing poset;
  --   s        columns s and s+1 in x are swapped;
  --   p        permutation defines current position of black checkers;
  --   locmap   Checker_Localization_Patterns.Column_Pattern(p,rows,cols).

  -- ON RETURN :
  --   x        coordinates for swap homotopy of first type.

  procedure Second_Swap_Plane
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );
  procedure Second_Swap_Plane
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );
  procedure Second_Swap_Plane
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );

  -- DESCRIPTION :
  --   In the second type of swap homotopy r+1 = R, where
  --   x is represented as a matrix of polynomials with coefficients
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   x        matrix of n rows and k rows;
  --   r        critical row, row of the descending black checker;
  --   dc       index of the descending black checker in specializing poset;
  --   s        columns s and s+1 in x are swapped;
  --   p        permutation defines current position of black checkers;
  --   locmap   Checker_Localization_Patterns.Column_Pattern(p,rows,cols).

  -- ON RETURN :
  --   x        coordinates for swap homotopy of second type.

  procedure Second_Swap_Plane
              ( file : in file_type;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );
  procedure Second_Swap_Plane
              ( file : in file_type;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );
  procedure Second_Swap_Plane
              ( file : in file_type;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix );

  -- DESCRIPTION :
  --   In the second type of swap homotopy r+1 = R, where
  --   x is represented as a matrix of polynomials with coefficients
  --   in standard double, double double, or quad double precision.
  --   Intermediate output is written to file.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted, then silent;
  --   x        matrix of n rows and k rows;
  --   r        critical row, row of the descending black checker;
  --   dc       index of the descending black checker in specializing poset;
  --   s        columns s and s+1 in x are swapped;
  --   p        permutation defines current position of black checkers;
  --   locmap   Checker_Localization_Patterns.Column_Pattern(p,rows,cols).

  -- ON RETURN :
  --   x        coordinates for swap homotopy of second type.

  function Stay_Moving_Plane
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector )
              return Standard_Complex_Poly_Matrices.Matrix;
  function Stay_Moving_Plane
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector )
              return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Stay_Moving_Plane
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector )
              return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a moving k-plane in n-space for use in a stay homotopy,
  --   represented as a matrix of polynomials with coefficients in
  --   standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   n        number of rows of the matrix on return;
  --   k        number of columns of the matrix on return;
  --   r        critical row;
  --   p        permutation defines current position of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers.

  function Swap_Moving_Plane
              ( n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector )
              return Standard_Complex_Poly_Matrices.Matrix;
  function Swap_Moving_Plane
              ( n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector )
              return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Swap_Moving_Plane
              ( n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector )
              return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a moving k-plane in n-space for use in a swap homotopy,
  --   represented as a matrix of polynomials with coefficients in
  --   standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   n        number of rows of the matrix on return;
  --   k        number of columns of the matrix on return;
  --   r        critical row;
  --   s        column to be swapped;
  --   q        parent permutation in specializing poset;
  --   p        permutation defines current position of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers.

  function Swap_Moving_Plane
              ( file : in file_type; n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector )
              return Standard_Complex_Poly_Matrices.Matrix;
  function Swap_Moving_Plane
              ( file : in file_type; n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector )
              return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Swap_Moving_Plane
              ( file : in file_type; n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector )
              return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a moving k-plane in n-space for use in a swap homotopy,
  --   represented as a matrix of polynomials with coefficients in
  --   standard double, double double, or quad double precision.
  --   This version writes intermediate output to file.

  -- ON ENTRY :
  --   file     for intermediate diagnostics, if omitted, then silent;
  --   n        number of rows of the matrix on return;
  --   k        number of columns of the matrix on return;
  --   r        critical row;
  --   s        column to be swapped;
  --   q        parent permutation in specializing poset;
  --   p        permutation defines current position of black checkers;
  --   rows     row indices for the location of the white checkers;
  --   cols     column indices for the location of the white checkers.

end Checker_Homotopies;
