with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Integer_Matrices;
with Standard_Complex_VecVecVecs;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Double_Lseries_Polynomials;        use Double_Lseries_Polynomials;

package Double_Lseries_Newton_Steps is

-- DESCRIPTION :
--   Runs Newton's method on Laurent series in double precision.

  procedure Make_Series
              ( sol : in Standard_Complex_Vectors.Vector;
                deg : in integer32;
                lead : out Standard_Integer_Vectors.Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Returns a regular power series of degree deg,
  --   with as leading coefficients the coordinates in sol.
  --   Allocates all space for cffs.

  -- REQUIRED : 
  --   lead'range = sol'range.

  procedure Make_Series
              ( sol : in Laur_Sys; deg : in integer32;
                lead : out Standard_Integer_Vectors.Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Returns a regular power series of degree deg,
  --   with as leading terms the coordinates in sol.
  --   Allocates all space for cffs.

  -- REQUIRED : 
  --   lead'range = sol'range.

  procedure Set_Leading_Exponents 
              ( lead : in out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Interactive setup of the leading exponents.
  --   Prompts for every element of lead.

  function Max_Norm ( v : Standard_Complex_VecVecs.Link_to_VecVec )
                    return double_float;

  -- DESCRIPTION :
  --   Returns the maximum of all max norms of the vectors in v.

  procedure Newton_Step
              ( deg : in integer32;
                p : in Table_Vector; jp : in Table_Vector_Array;
                xlead : in out Standard_Integer_Vectors.Vector;
                xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                ylead : in out Standard_Integer_Vectors.Vector;
                ycffs : in out Standard_Complex_VecVecs.Link_to_VecVec;
                Alead : in out Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Blead : in out Standard_Integer_Matrices.Matrix;
                Bcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                dxlead : in out Standard_Integer_Vectors.Vector;
                dxcffs : in out Standard_Complex_VecVecs.Link_to_VecVec;
                rlead : in out Standard_Integer_Vectors.Vector;
                rcffs : in out Standard_Complex_VecVecs.Link_to_VecVec;
                dxnrm,pxnrm : out double_float;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Does one step with Newton's method.

  -- REQUIRED : ycffs, dxcffs, Acffs, and Bcffs are allocated.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   p        table representation of a series system;
  --   jp       table representation of the Jacobian matrix;
  --   xlead    leading exponents of the current solution series;
  --   xcffs    coefficient vector of the current solution series;
  --   verbose  flag for the verbosity.

  -- ON RETURN :
  --   dxlead   leading exponents of the update to the solution series;
  --   dxcffs   coefficient vector of the update to the solution series;
  --   Alead    leading exponents of the LU factorization
  --   Acffs    factors of the Jacobian matrix evaluated at (xlead, xcffs);
  --   Blead    leading exponents of the evaluated Jacobian matrix;
  --   Bcffs    Jacobian matrix evaluated at (xlead, xcffs).
  --   dxnrm    maximum norm of the update, forward error estimate;
  --   pxnrm    maximum norm of the residual, backward error estimate.

  procedure Run_Newton_Steps
              ( deg : in integer32;
                tv : in Table_Vector; tva : in Table_Vector_Array;
                xlead : in out Standard_Integer_Vectors.Vector;
                xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                dxnrm,pxnrm : out double_float;
                numit : out integer32; maxit : in integer32 := 4;
                dxtol : in double_float := 1.0E-8;
                pxtol : in double_float := 1.0E-8;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Wraps the allocation of the coefficients
  --   and runs a number of Newton steps.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   p        table representation of a series system;
  --   jp       table representation of the Jacobian matrix;
  --   xlead    leading exponents of the current solution series;
  --   xcffs    coefficient vector of the current solution series;
  --   maxit    maximum number of steps;
  --   dxtol    tolerance on the forward error, on the update;
  --   pxtol    tolerance on the backward error, on the residual;
  --   verbose  flag for the verbosity.

  -- ON RETURN :
  --   xlead    leading exponents of the current solution series;
  --   xcffs    coefficient vector of the current solution series;
  --   dxnrm    maximum of the max norms of the update vectors;
  --   pxnrm    maximum of the residual vectors;
  --   numit    number of steps done.

end Double_Lseries_Newton_Steps;
