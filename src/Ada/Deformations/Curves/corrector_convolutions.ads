with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;

package Corrector_Convolutions is

-- DESCRIPTION :
--   The procedures in this package run Newton's method at a point,
--   on the leading coefficients of systems of convolution circuits,
--   to correct a solution after a predictor step.

-- MANAGEMENT OF LEADING COEFFICIENTS :

  procedure Allocate_Leading_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                lead : out Standard_Complex_VecVecs.Link_to_VecVec );
  procedure Allocate_Leading_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                lead : out DoblDobl_Complex_VecVecs.Link_to_VecVec );
  procedure Allocate_Leading_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                lead : out QuadDobl_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Allocates space for the leading coefficients in c.

  procedure Store_Leading_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                lead : in Standard_Complex_Vectors.Link_to_Vector );
  procedure Store_Leading_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                lead : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure Store_Leading_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                lead : in QuadDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Stores the leading coefficients of c into the vector lead.

  -- REQUIRED :
  --   c /= null and lead /= 0,
  --   lead'range = 0..c.nbr, where lead(0) is zero if c.cst is null.

  procedure Restore_Leading_Coefficients
              ( lead : in Standard_Complex_Vectors.Link_to_Vector;
                c : in Standard_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Restore_Leading_Coefficients
              ( lead : in DoblDobl_Complex_Vectors.Link_to_Vector;
                c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Restore_Leading_Coefficients
              ( lead : in QuadDobl_Complex_Vectors.Link_to_Vector;
                c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit );

  -- DESCRIPTION :
  --   Restores the values of the leading coefficients of c,
  --   using the values in the vector lead.

  -- REQUIRED :
  --   lead /= null and c /= null;
  --   lead'range = 0..c.nbr, if c.cst is zero, then it will
  --   not be assigned, regardless of the value of lead(0) is zero.

  procedure Store_Leading_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                lead : in Standard_Complex_VecVecs.Link_to_VecVec );
  procedure Store_Leading_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                lead : in DoblDobl_Complex_VecVecs.Link_to_VecVec );
  procedure Store_Leading_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                lead : in QuadDobl_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Stores the leading coefficients of the circuits in c
  --   into the vector of lead.

  -- REQUIRED :
  --   Allocate_Leading_Coefficients(c,lead) was executed.

  procedure Restore_Leading_Coefficients
              ( lead : in Standard_Complex_VecVecs.Link_to_VecVec;
                c : in Standard_Speelpenning_Convolutions.Circuits );
  procedure Restore_Leading_Coefficients
              ( lead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits );
  procedure Restore_Leading_Coefficients
              ( lead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits );

  -- DESCRIPTION :
  --   Using the leading coefficients in lead, restores the leading
  --   coefficients in the circuits c.

  -- REQUIRED :
  --   Store_Leading_Coefficients(c,lead) was executed.

-- MANAGEMENT OF ALL COEFFICIENTS :

  procedure Allocate_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                cff : out Standard_Complex_VecVecs.Link_to_VecVec );
  procedure Allocate_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                cff : out DoblDobl_Complex_VecVecs.Link_to_VecVec );
  procedure Allocate_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                cff : out QuadDobl_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Allocates space for all coefficients in the circuit.
  --   The constant coefficient is stored at index 0 in cff.

  procedure Allocate_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
              cff : out Standard_Speelpenning_Convolutions.Link_to_VecVecVec );
  procedure Allocate_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
              cff : out DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec );
  procedure Allocate_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
              cff : out QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Allocates space for all coefficients in the circuits c.

  procedure Store_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                cff : in Standard_Complex_VecVecs.Link_to_VecVec );
  procedure Store_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                cff : in DoblDobl_Complex_VecVecs.Link_to_VecVec );
  procedure Store_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                cff : in QuadDobl_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Stores the coefficients of c in the vector of vectors cff.

  -- REQUIRED : Allocated_Coefficients(c,cff) was executed.

  procedure Restore_Coefficients
              ( cff : in Standard_Complex_VecVecs.Link_to_VecVec;
                c : in Standard_Speelpenning_Convolutions.
                       Link_to_Circuit );
  procedure Restore_Coefficients
              ( cff : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                c : in DoblDobl_Speelpenning_Convolutions.
                       Link_to_Circuit );
  procedure Restore_Coefficients
              ( cff : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                c : in QuadDobl_Speelpenning_Convolutions.
                       Link_to_Circuit );

  -- DESCRIPTION :
  --   Restores the coefficients of c using the coefficients in cff.

  -- REQUIRED : Store_Coefficients(cff,c) was executed.

  procedure Store_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                cff : in Standard_Speelpenning_Convolutions.
                         Link_to_VecVecVec );
  procedure Store_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                cff : in DoblDobl_Speelpenning_Convolutions.
                         Link_to_VecVecVec );
  procedure Store_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                cff : in QuadDobl_Speelpenning_Convolutions.
                         Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Stores the coefficients in the circuits c in the vector of
  --   vectors of vectors cff.

  -- REQUIRED : Allocated_Coefficients(c,cff)) was executed.

  procedure Restore_Coefficients
              ( cff : in Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
                c : in Standard_Speelpenning_Convolutions.Circuits );
  procedure Restore_Coefficients
              ( cff : in DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits );
  procedure Restore_Coefficients
              ( cff : in QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits );

  -- DESCRIPTION :
  --   Restores the coefficients of c using the coefficients in cff.

  -- REQUIRED : Store_Coefficients(cff,c) was executed.

-- STEP COEFFICIENTS :
--   The correct is applied to check whether a step size is good.
--   The Step_Coefficient functions compute the constant coefficient
--   of the power series in the convolution circuits, given a value
--   for the step size t, which is evaluated at the power series.

  function Step_Coefficient
              ( c : Standard_Complex_Vectors.Vector; t : double_float )
              return Standard_Complex_Numbers.Complex_Number;
  function Step_Coefficient
              ( c : DoblDobl_Complex_Vectors.Vector; t : double_double )
              return DoblDobl_Complex_Numbers.Complex_Number;
  function Step_Coefficient
              ( c : QuadDobl_Complex_Vectors.Vector; t : quad_double )
              return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Given in c is the coefficient vector of a series
  --   and in t is some step size, t > 0.
  --   Returns the value of the series c evaluated at t,
  --   computed in double, double double, or quad double precision.

  procedure Step_Coefficient
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                t : in double_float );
  procedure Step_Coefficient
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                t : in double_double );
  procedure Step_Coefficient
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                t : in quad_double );

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t,
  --   in double, double double, or quad double precision.

  procedure Step_Coefficient
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                t : in double_float );
  procedure Step_Coefficient
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                t : in double_double );
  procedure Step_Coefficient
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                t : in quad_double );

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t,
  --   in double, double double, and quad double precision.

-- ONE STEP OF NEWTON'S METHOD WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );
  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true );
  procedure LU_Newton_Step
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );
  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true );
  procedure LU_Newton_Step
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );
  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true );

  -- DESCRIPTION :
  --   Does one Newton step with LU factorization,
  --   without an estimate for the inverse condition number,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write extra output to if verbose (optional);
  --   hom      convolution system for a homotopy;
  --   sol      an initial value for a solution at t = 0;
  --   verbose  flag to indicate if vectors need to be written,
  --            if a file is given on input.

  -- ON RETURN :
  --   sol      the updated solution;
  --   dx       the update vector applied to the solution;
  --   ipvt     pivoting information for the LU factorization;
  --   info     zero is all went well, if nonzero,
  --            then the matrix in hom.vm(0) may be singular.

-- ONE STEP OF NEWTON'S METHOD WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float );
  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := true );
  procedure LU_Newton_Step
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double );
  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double; verbose : in boolean := true );
  procedure LU_Newton_Step
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double );
  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double; verbose : in boolean := true );

  -- DESCRIPTION :
  --   Does one Newton step with LU factorization,
  --   with an estimate for the inverse condition number,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write extra output to if verbose (optional);
  --   hom      convolution system for a homotopy;
  --   sol      an initial value for a solution at t = 0;
  --   verbose  flag to indicate if vectors need to be written,
  --            if a file is given on input.

  -- ON RETURN :
  --   sol      the updated solution;
  --   dx       the update vector applied to the solution;
  --   ipvt     pivoting information for the LU factorization;
  --   rcond    estimate for the inverse of the condition number,
  --            if 1.0 + rcond = 1.0, then the Jacobian matrix is singular.

-- MANY STEPs OF NEWTON'S METHOD WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true );
  procedure LU_Newton_Steps
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true );
  procedure LU_Newton_Steps
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs several Newton steps using LU factorization,
  --   without an estimate for the inverse condition number,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write extra output to if verbose (optional);
  --   hom      convolution system for a homotopy;
  --   abh      homotopy with radii for coefficients for mixed residual;
  --   psv      psv.sol contains an initial value for a solution at t = 0;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on the update dx;
  --   extra    number of some extra corrector iterations;
  --   verbose  flag to indicate if vectors need to be written,
  --            if a file is given on input.

  -- ON RETURN :
  --   psv.sol  the updated solution;
  --   maxdx    magnitude of the largest number in ix;
  --   mixres   the mixed residual;
  --   dx       the update vector applied to the solution;
  --   ipvt     pivoting information for the LU factorization;
  --   info     zero is all went well, if nonzero,
  --            then the matrix in hom.vm(0) may be singular;
  --   fail     true if tolres > tol after maxit iterations,
  --            false otherwise.

-- MANY STEPs OF NEWTON'S METHOD WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float; fail : out boolean;
                extra : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true );
  procedure LU_Newton_Steps
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double; fail : out boolean;
                extra : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true );
  procedure LU_Newton_Steps
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double; fail : out boolean;
                extra : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs several Newton steps using LU factorization,
  --   with an estimate for the inverse condition number,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write extra output to if verbose (optional);
  --   hom      convolution system for a homotopy;
  --   abh      homotopy with radii for coefficients for mixed residual;
  --   psv      psv.sol contains an initial value for a solution at t = 0;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on the update dx;
  --   extra    number of some extra corrector iterations;
  --   verbose  flag to indicate if vectors need to be written,
  --            if a file is given on input.

  -- ON RETURN :
  --   psv.sol  the updated solution;
  --   maxdx    magnitude of the largest number in dx;
  --   mixres   the mixed residual;
  --   dx       the update vector applied to the solution;
  --   ipvt     pivoting information for the LU factorization;
  --   rcond    estimate for the inverse of the condition number,
  --            if 1.0 + rcond = 1.0, then the Jacobian matrix is singular;
  --   fail     true if tolres > tol after maxit iterations,
  --            false otherwise.

end Corrector_Convolutions;
