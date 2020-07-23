with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Speelpenning_Convolutions; 
with Standard_Coefficient_Circuits;
with Standard_Coefficient_Convolutions; 

package Standard_Predictor_Convolutions is

-- DESCRIPTION :
--   Exports data structures and procedures to construct a rational prediction
--   for the next solution on a path defined by a polynomial homotopy,
--   represented by convolution circuits, in double precision.

-- DATA STRUCTURE FOR LU :
--   The predictor type stores a solution in a vector of dimension dim,
--   of power series of degree deg, for a rational approximation of degrees
--   numdeg and dendeg, respectively of numerator and denominator.
--   The work space for the pivoting information of the linear system
--   solving of the series computation and the rational approximation.
--   The work space (wrk below) for Newton's method on power series
--   has a dimension equal to the number of equations in the homotopy.

-- REQUIRED : deg >= numdeg + dendeg + 2.

  type LU_Predictor ( dim,deg,numdeg,dendeg : integer32 ) is record
    sol : Standard_Complex_VecVecs.VecVec(1..dim);        -- degree deg series
    wrk : Standard_Complex_Vectors.Link_to_Vector;        -- Newton work space
    numcff : Standard_Complex_VecVecs.VecVec(1..dim);     -- numerator coeffs
    dencff : Standard_Complex_VecVecs.VecVec(1..dim);     -- denominator coeffs
    newtpiv : Standard_Integer_Vectors.Vector(1..dim);    -- pivots for Newton
    padepiv : Standard_Integer_Vectors.Vector(1..dendeg); -- pivots for Pade
    mat : Matrix(1..dendeg,1..dendeg); -- matrix for the rational approximation
    rhs : Standard_Complex_Vectors.Vector(1..dendeg); -- right hand side
   -- work space for the inlined Newton on linear systems of series
    rc,ic : Standard_Floating_VecVecs.Link_to_VecVec;       -- cols of A(0)
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec; -- row space
    rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;       -- rhs
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;       -- work vectors
  end record;

  type Link_to_LU_Predictor is access LU_Predictor;

  type LU_Predictor_Array is
    array ( integer32 range <> ) of Link_to_LU_Predictor;

-- DATA STRUCTURE FOR SVD :
--   The predictor with the Singular Value Decomposition works for systems
--   that are overdetermined, with more equations than variables.
--   The number of equations is stored in the neq attribute,
--   the number of variables is stored in the dim attribute.
--   The number of variables plus one is the dim1 attribute.
--   The update in Newton's method is stored twice, once as dx,
--   as a vector of power series, and once as xd, in linearized format,
--   as a series with vector coefficients. 

  type SVD_Predictor ( neq,dim,dim1,deg,numdeg,dendeg : integer32 ) is record
    sol : Standard_Complex_VecVecs.VecVec(1..dim);        -- degree deg series
    wrk : Standard_Complex_Vectors.Link_to_Vector;        -- of range 1..neq
    dx : Standard_Complex_VecVecs.VecVec(1..dim);         -- vector of series
    xd : Standard_Complex_VecVecs.VecVec(0..deg);         -- series of vectors
    svl : Standard_Complex_Vectors.Vector(1..dim1);       -- singular values
    U : Standard_Complex_Matrices.Matrix(1..neq,1..neq);  -- U in SVD
    V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);  -- V in SVD
    ewrk : Standard_Complex_Vectors.Link_to_Vector;       -- of range 1..dim
    numcff : Standard_Complex_VecVecs.VecVec(1..dim);     -- numerator coeffs
    dencff : Standard_Complex_VecVecs.VecVec(1..dim);     -- denominator coeffs
    padepiv : Standard_Integer_Vectors.Vector(1..dendeg); -- pivots for Pade
    mat : Matrix(1..dendeg,1..dendeg); -- matrix for the rational approximation
    rhs : Standard_Complex_Vectors.Vector(1..dendeg); -- right hand side
  end record;

  type Link_to_SVD_Predictor is access SVD_Predictor;

  type SVD_Predictor_Array is
    array ( integer32 range <> ) of Link_to_SVD_Predictor;

  type Predictor_Type is (LU, SVD);

  type Predictor ( kind : Predictor_Type := SVD ) is record
    case kind is
      when LU  => ludata : Link_to_LU_Predictor;
      when SVD => svdata : Link_to_SVD_Predictor;
    end case;
  end record;

  type Predictor_Array is array ( integer32 range <> ) of Predictor;

  type Predictor_Vectors ( dim,neq : integer32 ) is record
    sol : Standard_Complex_Vectors.Vector(1..dim);    -- solution work space
    radsol : Standard_Complex_Vectors.Vector(1..dim); -- radii of sol(k)
    res : Standard_Complex_Vectors.Vector(1..neq);    -- residual work space
    radres : Standard_Complex_Vectors.Vector(1..neq); -- evaluated at radsol
  end record;

  type Link_to_Predictor_Vectors is access Predictor_Vectors;

  type Predictor_Vectors_Array is
    array ( integer32 range <> ) of Link_to_Predictor_Vectors;

-- DATA STRUCTURE FOR CURVATURE :
--   The singular values of the Hessians are computed to estimate the
--   distance from the current solution to the nearest path.
--   With the data structure we store the work space and the largest
--   singular values of the Hessians.  At spot zero, we keep the smallest
--   singular value of the Jacobian matrix at the solution.

  type SVD_Hessians ( dim,dim1 : integer32 ) is record   -- dim1 = dim+1
    H : Standard_Complex_Matrices.Matrix(1..dim,1..dim); -- work Hessian
    U : Standard_Complex_Matrices.Matrix(1..dim,1..dim); -- work U of SVD
    V : Standard_Complex_Matrices.Matrix(1..dim,1..dim); -- work V of SVD
    svl : Standard_Complex_Vectors.Vector(1..dim1);      -- work values
    ewrk : Standard_Complex_Vectors.Vector(1..dim);      -- work values
    vals : Standard_Complex_Vectors.Vector(0..dim);      -- singular values
    work : Standard_Complex_Vectors.Vector(1..dim);      -- work space for SVD
    first : boolean;                                     -- first time flag
   -- work space for the inlined singular value computation
    sr,si : Standard_Floating_Vectors.Link_to_Vector;   -- singular values
    er,ei : Standard_Floating_Vectors.Link_to_Vector;   -- real and imag error
    xrv,xiv : Standard_Floating_VecVecs.Link_to_VecVec; -- columns of H
    urv,uiv : Standard_Floating_VecVecs.Link_to_VecVec; -- columns of U
    vrv,viv : Standard_Floating_VecVecs.Link_to_VecVec; -- columns of V
    wr,wi : Standard_Floating_Vectors.Link_to_Vector;   -- work space
  end record;

  type Link_to_SVD_Hessians is access SVD_Hessians;

  type SVD_Hessians_Array is
    array ( integer32 range <> ) of Link_to_SVD_Hessians;

-- CONSTRUCTORS :

  function Create ( sol : Standard_Complex_Vectors.Vector;
	            neq,deg,numdeg,dendeg : integer32;
                    inlined : boolean := true ) return LU_Predictor;
  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32;
                    inlined : boolean := true )
                  return Link_to_LU_Predictor;
  function Create ( sol : Standard_Complex_Vectors.Vector;
	            neq,deg,numdeg,dendeg : integer32 ) return SVD_Predictor;
  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 )
                  return Link_to_SVD_Predictor;

  -- DESCRIPTION :
  --   Given a solution vector and dimensions,
  --   makes the power series for sol and allocate work space.

  -- ON ENTRY :
  --   sol      a solution vector;
  --   neq      number of equations in the homotopy;
  --   deg      degree of the power series;
  --   numdeg   degree of the numerator of the rational approximation;
  --   dendeg   degree of the denominator of the rational approximation;
  --   inlined  by default, space will be allocated for the inlined
  --            versions of the linearization.

  -- ON RETURN :
  --   data allocated to run the Newton-Fabry-Pade predictor.

  function Create ( p : Link_to_LU_Predictor ) return Predictor;
  function Create ( p : Link_to_SVD_Predictor ) return Predictor;

  -- DESCRIPTION :
  --   Returns an instance of the predictor of the corresponding kind.

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32;
                    kind : Predictor_Type; 
                    inlined : boolean := true ) return Predictor;

  -- DESCRIPTION :
  --   Given solution vector, dimensions, and the kind,
  --   returns the corresponding predictor.

  function Create ( nbr : integer32;
                    sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32;
                    kind : Predictor_Type;
                    inlined : boolean := true ) return Predictor_Array;

  -- DESCRIPTION :
  --   Given solution vector, dimensions, and the kind,
  --   returns the corresponding predictor array of range 1..nbr.

  function Create ( nbr,dim,neq : integer32 ) return Predictor_Vectors_Array;

  -- DESCRIPTION :
  --   Returns a vector of range 1..nbr with Predictor_Vectors allocated
  --   for dimension dim and number of equations neq.

  function Create ( dim : integer32 ) return Link_to_SVD_Hessians;

  -- DESCRIPTION :
  --   Returns the allocated and initialized data structures for
  --   the computation of the curvature step size.

  function Create ( nbr,dim : integer32 ) return SVD_Hessians_Array;

  -- DESCRIPTION :
  --   Returns an array of range 1..nbr of data structures of dimension dim
  --   for the computation of the curvature step size.

-- AUXILIARY PREDICTOR PROCEDURES FOR SETUP :

  procedure Set_Lead_Coefficients
              ( p : in Predictor;
                s : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sets the leading coefficients of the data in the predictor p
  --   with the values in the vector s. 
  --   Resets all higher order coefficients to zero.

  procedure Newton_Fabry_Report 
              ( file : in file_type;
                nbrit : in integer32; absdx : in double_float;
                fail : in boolean;
                z : in Standard_Complex_Numbers.Complex_Number;
                r,err,step : in double_float;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                output : in boolean );

  -- DESCRIPTION :
  --   Writes the results of the Newton-Fabry-Pade predictor,
  --   in double precision.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbrit    number of iterations done with Newton's method;
  --   absdx    absolute value of last corrector;
  --   fail     true if required tolerance was not reached;
  --   z        estimate for the closest singularity;
  --   r        radius of z;
  --   err      error on the estimate z;
  --   step     the pole step, equals beta1*r;
  --   numcff   coefficients of the numerator of the Pade approximants;
  --   dencff   coefficients of the denominator of the Pade approximants;
  --   output   if true, then numcff and dencff are written to screen.

-- NEWTON-FABRY ON COEFFICIENT CONVOLUTIONS :

  procedure Newton_Fabry
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; vrblvl : in integer32 := 0 );
  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with LU,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision,
  --   on a system of coefficient convolution circuits.

  -- ON ENTRY :
  --   file     for writing data to if output (optional);
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   output   flag to indicate data output during computations,
  --            if a file is provided;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the last corrector step;
  --   fail     indicates failure status;
  --   z        closest singularity estimated by Fabry's theorem;
  --   rad      estimates radius of convergence of the series;
  --   err      error estimate on the location of z.

  procedure Newton_Fabry
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; vrblvl : in integer32 := 0 );
  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   file     for writing data to if output (optional);
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   output   flag to indicate data output during computations,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the last corrector step;
  --   rcond    estimate for the inverse of the condition number;
  --   fail     indicates failure status;
  --   z        closest singularity estimated by Fabry's theorem;
  --   rad      estimates radius of convergence of the series;
  --   err      error estimate on the location of z.

-- NEWTON-FABRY ON COMPLEX CONVOLUTION CIRCUITS :

  procedure Newton_Fabry
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; vrblvl : in integer32 := 0 );
  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with LU,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   file     for writing data to if output (optional);
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   output   flag to indicate data output during computations,
  --            if a file is provided;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the last corrector step;
  --   fail     indicates failure status;
  --   z        closest singularity estimated by Fabry's theorem;
  --   rad      estimates radius of convergence of the series;
  --   err      error estimate on the location of z.

  procedure Newton_Fabry
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; vrblvl : in integer32 := 0 );
  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   file     for writing data to if output (optional);
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   output   flag to indicate data output during computations,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the last corrector step;
  --   rcond    estimate for the inverse of the condition number;
  --   fail     indicates failure status;
  --   z        closest singularity estimated by Fabry's theorem;
  --   rad      estimates radius of convergence of the series;
  --   err      error estimate on the location of z.

  function Distance ( svh : in Link_to_SVD_Hessians ) return double_float;

  -- DESCRIPTION :
  --   Returns the estimate to the distance to the nearest solution,
  --   based on the singular values in svh.

-- HESSE-PADE ON COEFFICIENT CIRCUITS :

  procedure Cached_Singular_Values
              ( cfs : in Standard_Coefficient_Circuits.Link_to_System;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes the singular values of the Hessians with caching for
  --   lower degree polynomials if not svh.first.
  --   The other parameters are the same as in the Hesse_Pade procedure.

  procedure Hesse_Pade
              ( cfs : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_LU_Predictor;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                vrblvl : in integer32 := 0 );
  procedure Hesse_Pade
              ( file : in file_type;
                cfs : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_LU_Predictor;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the singular values of the Hessians and estimates the
  --   distance to the closest path to determine the step size.

  -- ON ENTRY :
  --   file     to write the output to if verbose (optional);
  --   cfs      coefficient circuit system;
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   xr       real parts of the solution vector;
  --   xi       imaginary parts of the solution vector;
  --   vh       space allocated for dim matrices, dim = dimension,
  --            all matrices have 1..dim for range(1) and range(2);
  --   svls     svls(0) contains the singular values of s.jm, and
  --            svls(k) contains the singular values of vh(k),
  --            for k in vh'range.
  --   beta2    multiplication factor for the curvature step;
  --   verbose  flag for intermediate numerical output,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   res      solution error estimated by Pade approximants in prd,
  --            the range of res is cfs.crc'range;
  --   nrm      2-norm of res;
  --   step     computed curvature step.

  procedure Hesse_Pade
              ( cfs : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                vrblvl : in integer32 := 0 );
  procedure Hesse_Pade
              ( file : in file_type;
                cfs : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the singular values of the Hessians and estimates the
  --   distance to the closest path to determine the step size.

  -- ON ENTRY :
  --   file     to write the output to if verbose (optional);
  --   cfs      coefficient circuit system
  --   prd      predictor data for SVD Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   xr       real parts of the solution vector;
  --   xi       imaginary parts of the solution vector;
  --   vh       space allocated for dim matrices, dim = dimension,
  --            all matrices have 1..dim for range(1) and range(2);
  --   svls     svls(0) contains the singular values of s.jm, and
  --            svls(k) contains the singular values of vh(k),
  --            for k in vh'range.
  --   beta2    multiplication factor for the curvature step;
  --   verbose  flag for intermediate numerical output,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   res      solution error estimated by Pade approximants in prd,
  --            the range of res is cfs.crc'range;
  --   nrm      2-norm of res;
  --   step     computed curvature step.

-- HESSE-PADE ON COMPLEX CONVOLUTION CIRCUITS :

  procedure Second
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the Hessians for the convolution circuits in hom
  --   at the solution sol and stores the results in svh.

  procedure Hesse_Pade
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                vrblvl : in integer32 := 0 );
  procedure Hesse_Pade
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the singular values of the Hessians and estimates the
  --   distance to the closest path to determine the step size.

  -- ON ENTRY :
  --   file     to write the output to if verbose (optional);
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   sol      the leading coefficients of the solution series;
  --   beta2    multiplication factor for the curvature step;
  --   verbose  flag for intermediate numerical output,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   res      solution error estimated by Pade approximants in prd,
  --            the range of res is hom.crc'range;
  --   nrm      2-norm of res;
  --   step     computed curvature step.

  procedure Hesse_Pade
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                vrblvl : in integer32 := 0 );
  procedure Hesse_Pade
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the singular values of the Hessians and estimates the
  --   distance to the closest path to determine the step size.

  -- ON ENTRY :
  --   file     to write the output to if verbose (optional);
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for SVD Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   sol      the leading coefficients of the solution series;
  --   beta2    multiplication factor for the curvature step;
  --   verbose  flag for intermediate numerical output,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   res      solution error estimated by Pade approximants in prd,
  --            the range of res is hom.crc'range;
  --   nrm      2-norm of res;
  --   step     computed curvature step.

-- PREDICTOR FEEDBACK LOOPS ON COEFFICIENT SYSTEMS :

  procedure EvalCoeff
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                t : in double_float );

  -- DESCRIPIONT :
  --   Evaluates the power series coefficient at the circuits in hom
  --   at t and stores the evaluated coefficients in cfh.

  procedure AbsVal ( cfh : in Standard_Coefficient_Circuits.Link_to_System );

  -- DESCRIPTION :
  --   Replaces all coefficients in the circuits of cfh by their radii.

  procedure EvalCffRad
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                t : in double_float );

  -- DESCRIPIONT :
  --   Evaluates the power series coefficient at the circuits in hom
  --   at t and stores the evaluated coefficients in cfh.
  --   Updates the radii of the coefficients of the circuits in abh
  --   with the new coefficients of cfh.
  --   This procedure combines the above EvalCoeff and AbsVal procedures.

  procedure Predictor_Feedback
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                psv : in out Predictor_Vectors;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                step : in out double_float; minstep,alpha : in double_float;
                nrm,mixres : out double_float; nbfail : out integer32;
                vrblvl : in integer32 := 0 );
  procedure Predictor_Feedback
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                psv : in out Predictor_Vectors;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                step : in out double_float; minstep,alpha : in double_float;
                nrm,mixres : out double_float; nbfail : out integer32;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given a step size, runs a predictor feedback loop,
  --   cutting the step size each time in half,
  --   until the mixed residual is less than the tolerance,
  --   or when the step becomes smaller than minstep.

  -- REQUIRED : xr'range = xi'range = psv.sol'range.

  -- ON ENTRY :
  --   file     to write extra output to if verbose (optional);
  --   hom      homotopy convolution circuit system
  --   cfh      coefficient circuits for the homotopy;
  --   xr       work space allocated for the real parts of solution vectors;
  --   xi       work space allocated for the imag parts of solution vectors;
  --   psv      work space for solution vectors and residuals;
  --   numcff   coefficients of the numerator of the Pade approximants;
  --   dencff   coefficients of the denominator of the Pade approximants;
  --   step     current step size;
  --   minstep  the minimum step size;
  --   alpha    tolerance on the mixed predictor residual;
  --   verbose  flag for extra output if true,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   step     shorter step size if nbfail > 0;
  --   psv      psv.sol prediction by evaluation of Pade approximants,
  --            psv.radsol contains radii of the predicted solution,
  --            psv.res is the evaluation of the predicted solution, and
  --            psv.radres is the evaluation of psv.radsol;
  --   nrm      max norm of the components in res;
  --   mixres   mixed residual.

-- PREDICTOR FEEDBACK LOOPS ON CONVOLUTION SYSTEMS :

  procedure Predictor_Feedback
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Predictor_Vectors;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                step : in out double_float; minstep,alpha : in double_float;
                nrm,mixres : out double_float; nbfail : out integer32;
                vrblvl : in integer32 := 0 );
  procedure Predictor_Feedback
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Predictor_Vectors;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                step : in out double_float; minstep,alpha : in double_float;
                nrm,mixres : out double_float; nbfail : out integer32;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given a step size, runs a predictor feedback loop,
  --   cutting the step size each time in half,
  --   until the mixed residual is less than the tolerance,
  --   or when the step becomes smaller than minstep.

  -- ON ENTRY :
  --   file     to write extra output to if verbose (optional);
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coefficients, for mixed residuals;
  --   psv      work space for solution vectors and residuals;
  --   numcff   coefficients of the numerator of the Pade approximants;
  --   dencff   coefficients of the denominator of the Pade approximants;
  --   step     current step size;
  --   minstep  the minimum step size;
  --   alpha    tolerance on the mixed predictor residual;
  --   verbose  flag for extra output if true,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   step     shorter step size if nbfail > 0;
  --   psv      psv.sol prediction by evaluation of Pade approximants,
  --            psv.radsol contains radii of the predicted solution,
  --            psv.res is the evaluation of the predicted solution, and
  --            psv.radres is the evaluation of psv.radsol;
  --   nrm      max norm of the components in res;
  --   mixres   mixed residual.

-- MAIN PREDICTOR PROCEDURES ON COEFFICIENT CONVOLUTION CIRCUITS :

  procedure LU_Prediction
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                vrblvl : in integer32 := 0 );
  procedure LU_Prediction
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with LU,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution.

  -- ON ENTRY :
  --   file     to write extra output to if output and/or verbose (optional);
  --   hom      homotopy convolution circuit system;
  --   cfh      coefficient circuits for the homotopy; 
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   rx       work space for the real parts of the series coefficients;
  --   ix       work space for the imaginary parts of the series coefficients;
  --   xr       work space allocated for the real parts of solution vectors;
  --   xi       work space allocated for the imag parts of solution vectors;
  --   vh       space allocated for dim matrices, dim = dimension,
  --            all matrices have 1..dim for range(1) and range(2);
  --   svls     svls(0) contains the singular values of s.jm, and
  --            svls(k) contains the singular values of vh(k),
  --            for k in vh'range.
  --   psv      work space for solution vectors and residuals;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   beta2    multiplication factor for the curvature step;
  --   maxstep  the maximum step size;
  --   minstep  the minimum step size;
  --   endt     the end value for the homotopy continuation parameter t;
  --   acct     accumulated sum of all successful steps, equals the
  --            current value of the homotopy continuation parameter t;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate data output during computations,
  --            if a file is given on input;
  --   verbose  flag for intermediate numerical output,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   psv      psv.sol contains the predicted solution,
  --            psv.radsol has the radii of the psv.sol components,
  --            psv.res is the residual of psv.sol, and
  --            psv.radres contains the evaluation at psv.radsol;
  --   fail     indicates failure status;
  --   step     the step size;
  --   nbpole   updated number of times pole step was minimal;
  --   nbhess   updated number of times curve step was minimal;
  --   nbmaxm   updated number of times maximum step size was minimal.

  procedure SVD_Prediction
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                vrblvl : in integer32 := 0 );
  procedure SVD_Prediction
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution.

  -- ON ENTRY :
  --   file     to write extra output to if output and/or verbose (optional);
  --   hom      homotopy convolution circuit system;
  --   cfh      coefficient circuits for the homotopy; 
  --   prd      predictor data for SVD Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   rx       work space for the real parts of the series coefficients;
  --   ix       work space for the imaginary parts of the series coefficients;
  --   xr       work space allocated for the real parts of solution vectors;
  --   xi       work space allocated for the imag parts of solution vectors;
  --   vh       space allocated for dim matrices, dim = dimension,
  --            all matrices have 1..dim for range(1) and range(2);
  --   svls     svls(0) contains the singular values of s.jm, and
  --            svls(k) contains the singular values of vh(k),
  --            for k in vh'range.
  --   psv      work space for solution vectors and residuals;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   beta2    multiplication factor for the curvature step;
  --   maxstep  the maximum step size;
  --   minstep  the minimum step size;
  --   endt     the end value for the homotopy continuation parameter t;
  --   acct     accumulated sum of all successful steps, equals the
  --            current value of the homotopy continuation parameter t;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate data output during computations,
  --            if a file is given on input;
  --   verbose  flag for intermediate numerical output,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   psv      psv.sol contains the predicted solution,
  --            psv.radsol has the radii of the psv.sol components,
  --            psv.res is the residual of psv.sol, and
  --            psv.radres contains the evaluation at psv.radsol;
  --   fail     indicates failure status;
  --   step     the step size;
  --   nbpole   updated number of times pole step was minimal;
  --   nbhess   updated number of times curve step was minimal;
  --   nbmaxm   updated number of times maximum step size was minimal.

-- MAIN PREDICTOR PROCEDURES ON COMPLEX CONVOLUTIONS :

  procedure LU_Prediction
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                vrblvl : in integer32 := 0 );
  procedure LU_Prediction
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with LU,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution.

  -- ON ENTRY :
  --   file     to write extra output to if output and/or verbose (optional);
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coefficients, for mixed residuals;
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   psv      work space for solution vectors and residuals;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   beta2    multiplication factor for the curvature step;
  --   maxstep  the maximum step size;
  --   minstep  the minimum step size;
  --   endt     the end value for the homotopy continuation parameter t;
  --   acct     accumulated sum of all successful steps, equals the
  --            current value of the homotopy continuation parameter t;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate data output during computations,
  --            if a file is given on input;
  --   verbose  flag for intermediate numerical output,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   psv      psv.sol contains the predicted solution,
  --            psv.radsol has the radii of the psv.sol components,
  --            psv.res is the residual of psv.sol, and
  --            psv.radres contains the evaluation at psv.radsol;
  --   fail     indicates failure status;
  --   step     the step size;
  --   nbpole   updated number of times pole step was minimal;
  --   nbhess   updated number of times curve step was minimal;
  --   nbmaxm   updated number of times maximum step size was minimal.

  procedure SVD_Prediction
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                vrblvl : in integer32 := 0 );
  procedure SVD_Prediction
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution.

  -- ON ENTRY :
  --   file     to write extra output to if output and/or verbose (optional);
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coefficients, for mixed residuals;
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   psv      work space for solution vectors and residuals;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   beta2    multiplication factor for the curvature step;
  --   maxstep  the maximum step size;
  --   minstep  the minimum step size;
  --   endt     the end value for the homotopy continuation parameter t;
  --   acct     accumulated sum of all successful steps, equals the
  --            current value of the homotopy continuation parameter t;
  --   nbpole   number of times the pole step was minimal;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate data output during computations,
  --            if a file is given on input;
  --   verbose  flag for intermediate numerical output,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   psv      psv.sol contains the predicted solution,
  --            psv.radsol has the radii of the psv.sol components,
  --            psv.res is the residual of psv.sol, and
  --            psv.radres contains the evaluation at psv.radsol;
  --   fail     indicates failure status;
  --   step     the step size;
  --   nbpole   updated number of times pole step was minimal;
  --   nbhess   updated number of times curve step was minimal;
  --   nbmaxm   updated number of times maximum step size was minimal.

-- DESTRUCTORS :

  procedure Clear ( p : in out Link_to_LU_Predictor );
  procedure Clear ( p : in out Link_to_SVD_Predictor );

  -- DESCRIPTION :
  --   Deallocates the memory for the series, the work space,
  --   and the rational approximation.

  procedure Clear ( p : in out Predictor );
  procedure Clear ( p : in out Predictor_Array );

  -- DESCRIPTION :
  --   Deallocates the memory allocated to p.

  procedure Clear ( p : in out Link_to_Predictor_Vectors );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by p.

  procedure Clear ( h : in out SVD_Hessians );
  procedure Clear ( h : in out Link_to_SVD_Hessians );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by h.

  procedure Clear ( p : in out LU_Predictor_Array );
  procedure Clear ( p : in out SVD_Predictor_Array );
  procedure Clear ( p : in out Predictor_Vectors_Array );
  procedure Clear ( h : in out SVD_Hessians_Array );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the array.

end Standard_Predictor_Convolutions;
