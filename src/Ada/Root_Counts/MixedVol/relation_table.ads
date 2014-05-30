with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Zero_Index_tree;                   use Zero_Index_tree;
with Index_Tree_LP;                     use Index_Tree_LP;

package Relation_Table is

-- DESCRIPTION :
--   This package offers routines to create a "relation table" which
--   records all edges on the lower hulls of all pairwise Minkowski sums.

  type Boolean_Matrix is
    array ( integer32 range <> , integer32 range <> ) of boolean;

  procedure RelTable
               ( nVar,nSpt : in integer32;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 RelTab : out Boolean_Matrix;
                 L0 : in out Link_to_L0_IML );

  -- DESCRIPTION :
  --   Sets up the relation table for the given lifted supports.
  --   This is the main routine which calls the other functions
  --   with prototypes listed below.
  --
  -- ON ENTRY :
  --   nVar      ambient dimension;
  --   nSpt      number of different supports;
  --   Spt       supports, in array of range 0..nSpt-1;
  --   SptIdx    the i-th supports starts at SptIdx[i] in Spt,
  --             SptIdx[nSpt] equals the total #points in Spt;
  --   lft       lifting for each point in Spt;
  --   RelTab    matrix of dimension equal to the total #points;
  --   L0        initialized level zero tree.
  --
  -- ON RETURN :
  --   RelTab    RelTab[i][j] is true if points i and j span an edge
  --             on the lower hull of the polytope (if i and j are
  --             from the same support set) or on the lower hull
  --             of the Minkowski sum of the polytopes spanned by
  --             the support set which contains point i and the
  --             other support set which contains point j,
  --             RelTab[i][j] is false otherwise,
  --             note that this matrix is symmetric;
  --   L0        updated tree with stored LP information.

  procedure RlTbLP2_a
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 info : out integer32 );

  -- DESCRIPTION :
  --   Applies LP for nondegenerate solutions to perform a 1-point test.
  --
  -- ON ENTRY :
  --   ma        number of involved constraints;
  --   na        number of variables (= nVar) plus one;
  --   a         matrix of inequality constraints,
  --             as many rows as total #points plus one,
  --             as many columns as nVar + 2;
  --   nv1       #variables plus one;
  --   c         vector of range 0..nVar+1,
  --             used to determine outgoing constraint;
  --   LPidx     vector of dimension total #points plus one,
  --             index to the involved constaints;
  --   Bidx      index vector of range 0..nVar+1;
  --   x         solution vector of range 0..nVar+1;
  --   Binv      square matrix of dimension nVar+2.
  --
  -- ON RETURN :
  --   Bidx      updated index vector;
  --   x         updated solution vector;
  --   Binv      updated matrix;
  --   info      information about the rank.

  procedure RlTbLP2_e
               ( ma,na,NumCol : in integer32;
                 a : in out Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32; 
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 info : out integer32 );

  -- DESCRIPTION :
  --   Applies LP to extend the point using a 1-point test,
  --   returns a nondegenerate solution for 1-point tests.
  --
  -- ON ENTRY :
  --   ma        number of involved constraints;
  --   na        number of variables (= nVar), original dimension;
  --   NumCol    total number of points in the support;
  --   a         matrix of inequality constraints,
  --             as many rows as total #points plus one,
  --             as many columns as nVar + 2;
  --   nv1       #variables plus one;
  --   LPidx     vector of dimension total #points plus one,
  --             index to the involved constraints;
  --   Bidx      index vector of range 0..nVar+1;
  --   x         solution vector of range 0..nVar+1;
  --   Binv      square matrix of dimension nVar+2.
  --
  -- ON RETURN :
  --   a         extended matrix with the "epsilon variable";
  --   Bidx      updated index vector;
  --   x         updated solution vector;
  --   Binv      updated matrix;
  --   info      information about the rank.

  procedure dlp2_1pt_i
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 FixPt,TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 RelTab : in out Boolean_Matrix );

  -- DESCRIPTION :
  --   Applies LP to perform a 1-point test in case the current support
  --   is not the first one, for the variable info < 0.
  --
  -- ON ENTRY :
  --   ma        number of involved constraints;
  --   na        number of variables (= nVar);
  --   a         matrix of inequality constraints,
  --             as many rows as total #points plus one,
  --             as many columns as nVar + 2;
  --   nv1       #variables plus one;
  --   c         vector of range 0..nVar+1,
  --             used to determine outgoing constraint;
  --   LPidx     vector of dimension total #points plus one;
  --   FixPt     index to a point in the supports;
  --   TstPt     index to an involved constraint;
  --   Bidx      index vector of range 0..nVar+1;
  --   x         solution vector of range 0..nVar+1;
  --   Binv      square matrix of dimension nVar+2;
  --   RelTab    matrix of dimension equal to the total #points.
  --
  -- ON RETURN :
  --   Bidx      updated index vector;
  --   x         updated solution vector;
  --   Binv      updated matrix;
  --   RelTab    updated relation table.

  procedure dlp2_1pt_s
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 FixPt,TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 RelTab : in out Boolean_Matrix;
                 L0 : in out Link_to_L0_IML );

  -- DESCRIPTION :
  --   Applies LP to perform a 1-point test for the first support,
  --   in case the variable info < 0.
  --
  -- ON ENTRY :
  --   ma        number of involved constraints;
  --   na        number of variables (= nVar);
  --   a         matrix of inequality constraints,
  --             as many rows as total #points plus one,
  --             as many columns as nVar + 2;
  --   nv1       #variables plus one;
  --   c         vector of range 0..nVar+1,
  --             used to determine outgoing constraint;
  --   LPidx     vector of dimension total #points plus one,
  --             index to the involved constraints;
  --   FixPt     index to a point in the supports;
  --   TstPt     index to an involved constraint;
  --   Bidx      index vector of range 0..nVar+1;
  --   x         solution vector of range 0..nVar+1;
  --   Binv      square matrix of dimension nVar+2;
  --   RelTab    matrix of dimension equal to the total #points;
  --   L0        current level zero tree.
  --
  -- ON RETURN :
  --   Bidx      updated index vector;
  --   x         updated solution vector;
  --   Binv      updated matrix;
  --   RelTab    updated relation table;
  --   L0        LP data with Bidx, x, and Binv.
   
  procedure dlp1_1pt_i
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 FixPt,TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 RelTab : in out Boolean_Matrix );

  -- DESCRIPTION :
  --   Applies LP to perform a 1-point test when the current support
  --   is not the first one, and for variable info >= 0.
  --
  -- ON ENTRY :
  --   ma        number of involved constraints;
  --   na        number of variables;
  --   a         matrix of inequality constraints,
  --             as many rows as total #points plus one,
  --             as many columns as nVar + 2;
  --   nv1       #variables plus one;
  --   c         vector of range 0..nVar+1,
  --             used to determine outgoing constraint;
  --   LPidx     vector of dimension total #points plus one,
  --             index to the involved constraints;
  --   FixPt     index to a point in the supports;
  --   TstPt     index to an involved constraint;
  --   Bidx      index vector of range 0..nVar+1;
  --   x         solution vector of range 0..nVar+1;
  --   Binv      square matrix of dimension nVar+2;
  --   RelTab    matrix of dimension equal to the total #points.
  --
  -- ON RETURN :
  --   Bidx      updated index vector;
  --   x         updated solution vector;
  --   Binv      updated matrix;
  --   RelTab    updated relation table.

  procedure dlp1_1pt_s
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nv1 : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 LPidx : in Standard_Integer_Vectors.Link_to_Vector;
                 FixPt,TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 RelTab : in out Boolean_Matrix;
                 L0 : in out Link_to_L0_IML );

  -- DESCRIPTION :
  --   Applies LP to perform a 1-point test for the first support
  --   and for variable info >= 0.
  --
  -- ON ENTRY :
  --   ma        number of involved constraints;
  --   na        number of variables;
  --   a         matrix of inequality constraints,
  --             as many rows as total #points plus one,
  --             as many columns as nVar + 2;
  --   nv1       #variables plus one;
  --   c         vector of range 0..nVar+1,
  --             used to determine outgoing constraint;
  --   LPidx     vector of dimension total #points plus one,
  --             index to the involved constraints;
  --   FixPt     index to a point in the supports;
  --   TstPt     index to an involved constraint;
  --   Bidx      index vector of range 0..nVar+1;
  --   x         solution vector of range 0..nVar+1;
  --   Binv      square matrix of dimension nVar+2;
  --   RelTab    matrix of dimension equal to the total #points;
  --   L0        current level zero tree.
  --
  -- ON RETURN :
  --   x         updated vector;
  --   Binv      updated matrix;
  --   RelTab    updated relation table;
  --   L0        saved LP data Bidx, x, and Binv.

end Relation_Table;
