with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Integer_Vectors;

package Linear_Programming is

-- DESCRIPTION :
--   This package contains some routines for solving the standard primal and
--   dual linear-optimization problems, by means of the simplex algorithm.

  procedure Primal_Simplex
                ( dic : in out Matrix; eps : in double_float;
                  in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                  nit : in out natural32; unbounded : out boolean );

  generic

    with procedure Report ( dic : in Matrix;
                            in_bas,out_bas : Standard_Integer_Vectors.Vector );

  procedure Generic_Primal_Simplex
                ( dic : in out Matrix; eps : in double_float;
                  in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                  nit : in out natural32; unbounded : out boolean );

  -- DESCRIPTION :
  --   This is a very simple implementation of the simplex procedure
  --   for solving the following problem:
  --
  --      max <c,x>
  --          <a,x> <= 0
  --   
  --   where x = (1,x1,x2,..,xn).
  --   The generic procedure allows to report on intermediate dictionaries.

  -- REQUIRED : The primal dictionary is already initialized.

  -- ON ENTRY :
  --   dic        the matrix for the initial dictionary;
  --   eps        constant float to determine wether number is zero;
  --   in_bas     initial unknowns in the basis;
  --   out_bas    initial unknowns out the basis;
  --   nit        counter for number of iterations.

  -- ON RETURN :
  --   dic        the modified matrix for the dictionary;
  --   in_bas     unknowns in the basis;
  --   out_bas    unknowns out the basis;
  --   nit        augmented with number of iterations made;
  --   unbounded  true when it is detected that the solution is unbounded.

  procedure Dual_Simplex 
                ( dic : in out Matrix; eps : in double_float;
                  in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                  nit : in out natural32; feasible : out boolean );

  generic

    with procedure Report ( dic : in Matrix;
                            in_bas,out_bas : Standard_Integer_Vectors.Vector );

  procedure Generic_Dual_Simplex 
                ( dic : in out Matrix; eps : in double_float;
                  in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                  nit : in out natural32; feasible : out boolean );

  -- DESCRIPTION :
  --   This is a very simple implementation of the simplex procedure
  --   for solving the following problem:
  --
  --      min <c,x>
  --          <a,x> >= 0
  --   
  --   where x = (1,x1,x2,..,xn).
  --   The generic procedure allows to report on intermediate dictionaries.

  -- REQUIRED : The dual dictionary is already initialized.

  -- ON ENTRY :
  --   dic        the matrix for the initial dictionary;
  --   eps        constant float to determine wether number is zero;
  --   in_bas     initial unknowns in the basis;
  --   out_bas    initial unknowns out the basis;
  --   nit        counter for number of iterations.

  -- ON RETURN :
  --   dic        the modified matrix for the dictionary; 
  --   in_bas     unknowns in the basis;
  --   out_bas    unknowns out the basis;
  --   nit        augmented with number of iterations made;
  --   feasible   is true when the problem is detected to be infeasible.

end Linear_Programming;
