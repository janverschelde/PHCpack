with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;

package Simplex_Pivoting is

-- DESCRIPTION :
--   This package collects the common pivoting routines used in the
--   simplex algorithms to solve linear programming problems.

-- COMMON CONSTANTS :

  eps : constant double_float := 1.0E-6;
  DBL_MAX : constant double_float := 1.0E+20;

  singular_base : exception; -- raised when base matrix is singular
  unbounded_LP : exception;  -- raised when LP problem is unbounded

-- PIVOTING OPERATIONS :

  procedure Search_Outgoing
               ( na : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 Binv : in Standard_Floating_Matrices.Link_to_Matrix;
                 vmax : out double_float; k : out integer32 );

  -- DESCRIPTION :
  --   Searches the constraint which will leave the base.

  -- ON ENTRY :
  --   na        dimension of the LP problem;
  --   c         vector used to determine the index;
  --   Binv      current inverse of the base matrix;

  -- ON RETURN :
  --   vmax      value computed to determine the index k;
  --   k         index of the constraint which will leave the base.

  procedure Search_Outgoing
               ( na : in integer32;
                 c : in Standard_Floating_Vectors.Link_to_Vector;
                 Bidx : in Standard_Integer_Vectors.Link_to_Vector;
                 Binv : in Standard_Floating_Matrices.Link_to_Matrix;
                 vmax : out double_float; k : out integer32 );

  -- DESCRIPTION :
  --   Searches the constraint which will leave the base,
  --   skipping those indices i for which Bidx(i) = -1.

  -- ON ENTRY :
  --   na        dimension of the LP problem;
  --   c         vector used to determine the index;
  --   Bidx      index vector to constraints currently involved;
  --   Binv      current inverse of the base matrix;

  -- ON RETURN :
  --   vmax      value computed to determin the index k;
  --   k         index of the constraint which will leave the base.

  procedure Search_Incoming
               ( ma,na,nv1,k : in integer32;
                 LPidx,Bidx : in Standard_Integer_Vectors.Link_to_Vector;
                 x : in Standard_Floating_Vectors.Link_to_Vector;
                 a,Binv : in Standard_Floating_Matrices.Link_to_Matrix;
                 sigj : out double_float; ell : out integer32 );

  -- DESCRIPTION :
  --   Searches for the new constraint coming in to the base.
  --   Raises the exception "unbounded_LP" when the LP problem is unbounded.

  -- ON ENTRY :
  --   ma        number of involved constraints;
  --   na        number of variables;
  --   nv1       number of variables plus one;
  --   k         index to the outgoing constraint;
  --   LPidx     vector of dimension equal to the total #points plus one;
  --   Bidx      index vector to the involved constraints;
  --   x         current solution vector;
  --   a         coefficients of the constraints in the LP problem;
  --   Binv      inverse of the base matrix.

  -- ON RETURN :
  --   sigj      value to be used to update the solution vector;
  --   ell       index of the incoming constraint.

  procedure Search_Incoming
               ( ma,na,nv1,k : in integer32;
                 LPidx,Bidx : in Standard_Integer_Vectors.Link_to_Vector;
                 x : in Standard_Floating_Vectors.Link_to_Vector;
                 v : in Standard_Floating_Vectors.Vector;
                 a,Binv : in Standard_Floating_Matrices.Link_to_Matrix;
                 sigj : out double_float; ell : out integer32 );

  -- DESCRIPTION :
  --   Searches for the new constraint coming in to the base.
  --   Uses v to determine the index and raises no exceptions,
  --   even when the LP problem might be unbounded.

  -- ON ENTRY :
  --   ma        number of involved constraints;
  --   na        number of variables;
  --   nv1       number of variables plus one;
  --   k         index to the outgoing constraint;
  --   LPidx     vector of dimension equal to the total #points plus one;
  --   Bidx      index vector to the involved constraints;
  --   x         current solution vector;
  --   v         used to the determine the index ell;
  --   a         coefficients of the constraints in the LP problem;
  --   Binv      inverse of the base matrix.

  -- ON RETURN :
  --   sigj      value to be used to update the solution vector;
  --   ell       index of the incoming constraint.

  procedure Search_Incoming
               ( ma,na,nVar,k : in integer32;
                 Bidx : in Standard_Integer_Vectors.Link_to_Vector;
                 x : in Standard_Floating_Vectors.Link_to_Vector;
                 a,Binv : in Standard_Floating_Matrices.Link_to_Matrix;
		 sigj : out double_float; ell : out integer32 );

  -- DESCRIPTION :
  --   Searches for the new constraint coming in to the base.
  --   Raises an exception when the LP problem is unbounded.

  -- ON ENTRY :
  --   ma        number of involved constraints;
  --   na        number of variables;
  --   nVar      number of variables;
  --   k         index to the outgoing constraint;
  --   Bidx      index vector to the involved constraints;
  --   x         current solution vector;
  --   a         coefficients of the constraints in the LP problem;
  --   Binv      inverse of the base matrix.

  -- ON RETURN :
  --   sigj      value to be used to update the solution vector;
  --   ell       index of the incoming constraint.

  procedure Update_Base
               ( Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                 Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 na,k,ell : in integer32;
                 a : in Standard_Floating_Matrices.Link_to_Matrix );

  -- DESCRIPTION :
  --   Updates the base matrix of the LP problem.
  --   Raises the exception "singular_base" when the base is singular.

  -- ON ENTRY :
  --   Bidx      index vector to the constraints involved;
  --   Binv      inverse of the base matrix;
  --   na        dimension of the LP problem;
  --   k         index of the outgoing constraint;
  --   ell       index of the incoming constraint;
  --   a         coefficient matrix of the LP problem.

  -- ON RETURN :
  --   Bidx      updated index vector;
  --   Binv      updated basis inverse.

end Simplex_Pivoting;
