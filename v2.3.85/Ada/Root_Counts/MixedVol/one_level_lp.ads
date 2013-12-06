with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Index_Tree_LP;                     use Index_tree_LP;

package One_Level_LP is

-- DESCRIPTION :
--    This packages offers routines to perform tests at one level,
--    as called by the main program "mixed_volume".
 
  procedure one_level_LP 
               ( Strt1Pt,End1Pt : in integer32;
                 PtIn : in out Standard_Integer_Vectors.Link_to_Vector;
                 LPdim : in integer32;
                 A : in Standard_Floating_Matrices.Link_to_Matrix;
                 nVar : in integer32;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 ItLp : in out Link_to_IT_LP );

  -- DESCRIPTION :
  --   This routine is called in the main program in "mixed_volume"
  --   and uses the other functions below in its implementation.
  --
  -- ON ENTRY :
  --   Strt1Pt   index where to start the 1-point tests;
  --   End1Pt    index where to end the 1-point tests;
  --   PtIn      one row of a matrix with as many columns as the
  --             total number of points in the supports;
  --   LPDim     current dimension of the LP problem,
  --             its value is used as "na" in the routines below;
  --   A         matrix of constraints at the current level;
  --   nVar      number of variables;
  --   x         solution vector of dimension nVar;
  --   Binv      basis inverse, matrix of dimension nVar;
  --   Bidx      index vector of dimension nVar;
  --   ItLp      index tree with data for LP problems.
  --
  -- ON RETURN :
  --   PtIn      updated row of indices;
  --   x         updated solution vector;
  --   Binv      updated basis inverse;
  --   Bidx      updated index vector;
  --   ItLp      updated index tree with data for LP problems.

  procedure dnulp2_a
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nVar : in integer32; 
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 info : out integer32 );

  -- DESCRIPTION :
  --   Auxiliary routine called by one_level_LP,
  --   invoked the first time when performing a 1-point test.
  --
  -- ON ENTRY :
  --   ma        #constraints involved in the LP problem;
  --   na        dimension of the LP problem;
  --   a         matrix with coefficients of the inequality constraints;
  --   nVar      number of variables;
  --   c         vector of dimension nVar+2,
  --             used to determine the outgoing constraint;
  --   Bidx      index vector of dimension nVar;
  --   x         solution vector of dimension nVar;
  --   Binv      basis inverse, matrix of dimension nVar.
  --
  -- ON RETURN :
  --   Bidx      updated index vector;
  --   x         updated solution vector;
  --   Binv      updated basis inverse;
  --   info      information about the dimension of the problem.

  procedure dlp2_1pts
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nVar : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 PtIn : in out Standard_Integer_Vectors.Link_to_Vector;
                 ItLp  : in out Link_to_IT_LP );

  -- DESCRIPTION :
  --   Auxiliary routine called by one_level_LP,
  --   invoked after the first time, when dnulp2_a returns info >= 0.
  --
  -- ON ENTRY :
  --   ma        #constraints involved in the LP problem;
  --   na        dimension of the LP problem;
  --   a         matrix with coefficients of the inequality constraints;
  --   nVar      number of variables;
  --   c         vector of dimension nVar+2,
  --             used to determine the outgoing constraint;
  --   TstPt     index of a point ranging between Strt1Pt+1 and End1Pt
  --             of the routine one_level_LP from above;
  --   Bidx      index vector of dimension nVar;
  --   x         solution vector of dimension nVar;
  --   Binv      basis inverse, matrix of dimension nVar;
  --   PtIn      one row of a matrix with as many columns as the
  --             total number of points in the supports;
  --   ItLp      index tree with data for LP problems.
  --
  -- ON RETURN :
  --   Bidx      updated index vector;
  --   x         updated solution vector;
  --   Binv      updated basis inverse;
  --   PtIn      vector with one updated entry;
  --   ItLp      updated index tree with data for LP problems.

  procedure dlp1_1pts
               ( ma,na : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix;
                 nVar : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector; 
                 TstPt : in integer32;
                 Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 x : in out Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 PtIn : in out Standard_Integer_Vectors.Link_to_Vector;
                 ItLp : in out Link_to_IT_LP );

  -- DESCRIPTION :
  --   Auxiliary routine called by one_level_LP,
  --   invoked after the first time, when dnulp2_a returns info < 0.
  --
  -- ON ENTRY :
  --   ma        #constraints involved in the LP problem;
  --   na        dimension of the LP problem;
  --   a         matrix with coefficients of the inequality constraints;
  --   nVar      number of variables;
  --   c         vector of dimension nVar+2,
  --             used to determine the outgoing constraint;
  --   TstPt     index of a point ranging between Strt1Pt+1 and End1Pt
  --             of the routine one_level_LP from above;
  --   Bidx      index vector of dimension nVar;
  --   x         solution vector of dimension nVar;
  --   Binv      basis inverse, matrix of dimension nVar;
  --   PtIn      one row of a matrix with as many columns as the
  --             total number of points in the supports;
  --   ItLp      index tree with data for LP problems.
  --
  -- ON RETURN :
  --   Bidx      updated index vector;
  --   x         updated solution vector;
  --   Binv      updated basis inverse;
  --   PtIn      vector with one updated entry;
  --   ItLp      updated index tree with data for LP problems.

  procedure Sort ( n : in integer32; 
                   a : in out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Sorts the first n numbers in the array a in ascending order.
  --
  -- ON ENTRY :
  --   n         dimension of the array a;
  --   a         sequence of n integer numbers.
  --
  -- ON RETURN :
  --   a         sequence sorted in ascending order.

end One_Level_LP;
