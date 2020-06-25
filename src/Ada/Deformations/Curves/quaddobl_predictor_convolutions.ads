with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;
with QuadDobl_Speelpenning_Convolutions; use QuadDobl_Speelpenning_Convolutions;

package QuadDobl_Predictor_Convolutions is

-- DESCRIPTION :
--   Exports data structures and procedures to construct a rational prediction
--   for the next solution on a path defined by a polynomial homotopy,
--   represented by convolution circuits, in quad double precision.

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
    sol : QuadDobl_Complex_VecVecs.VecVec(1..dim);        -- degree deg series
    wrk : QuadDobl_Complex_Vectors.Link_to_Vector;        -- Newton work space
    numcff : QuadDobl_Complex_VecVecs.VecVec(1..dim);     -- numerator coeffs
    dencff : QuadDobl_Complex_VecVecs.VecVec(1..dim);     -- denominator coeffs
    newtpiv : Standard_Integer_Vectors.Vector(1..dim);    -- pivots for Newton
    padepiv : Standard_Integer_Vectors.Vector(1..dendeg); -- pivots for Pade
    mat : Matrix(1..dendeg,1..dendeg); -- matrix for rational approximation
    rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg);
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
    sol : QuadDobl_Complex_VecVecs.VecVec(1..dim);        -- degree deg series
    wrk : QuadDobl_Complex_Vectors.Link_to_Vector;        -- of range 1..neq
    dx : QuadDobl_Complex_VecVecs.VecVec(1..dim);         -- vector of series
    xd : QuadDobl_Complex_VecVecs.VecVec(0..deg);         -- series of vectors
    svl : QuadDobl_Complex_Vectors.Vector(1..dim1);       -- singular values
    U : QuadDobl_Complex_Matrices.Matrix(1..neq,1..neq);  -- U in SVD
    V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);  -- V in SVD
    ewrk : QuadDobl_Complex_Vectors.Link_to_Vector;       -- of range 1..dim
    numcff : QuadDobl_Complex_VecVecs.VecVec(1..dim);     -- numerator coeffs
    dencff : QuadDobl_Complex_VecVecs.VecVec(1..dim);     -- denominator coeffs
    padepiv : Standard_Integer_Vectors.Vector(1..dendeg); -- pivots for Pade
    mat : Matrix(1..dendeg,1..dendeg); -- matrix for the rational approximation
    rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg); -- right hand side
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
    sol : QuadDobl_Complex_Vectors.Vector(1..dim);    -- solution work space
    radsol : QuadDobl_Complex_Vectors.Vector(1..dim); -- radii of sol(k)
    res : QuadDobl_Complex_Vectors.Vector(1..neq);    -- residual work space
    radres : QuadDobl_Complex_Vectors.Vector(1..neq); -- evaluated at radsol
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
    H : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim); -- work Hessian
    U : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim); -- work U of SVD
    V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim); -- work V of SVD
    svl : QuadDobl_Complex_Vectors.Vector(1..dim1);      -- work values
    ewrk : QuadDobl_Complex_Vectors.Vector(1..dim);      -- work values
    vals : QuadDobl_Complex_Vectors.Vector(0..dim);      -- singular values
    work : QuadDobl_Complex_Vectors.Vector(1..dim);      -- work space for SVD
  end record;

  type Link_to_SVD_Hessians is access SVD_Hessians;

  type SVD_Hessians_Array is
    array ( integer32 range <> ) of Link_to_SVD_Hessians;

-- CONSTRUCTORS :

  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) return LU_Predictor;
  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 )
                  return Link_to_LU_Predictor;
  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) return SVD_Predictor;
  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
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
  --   dendeg   degree of the denominator of the rational approximation.

  -- ON RETURN :
  --   data allocated to run the Newton-Fabry-Pade predictor.

  function Create ( p : Link_to_LU_Predictor ) return Predictor;
  function Create ( p : Link_to_SVD_Predictor ) return Predictor;

  -- DESCRIPTION :
  --   Returns an instance of the predictor of the corresponding kind.

  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32;
                    kind : Predictor_Type ) return Predictor;

  -- DESCRIPTION :
  --   Given solution vector, dimensions, and the kind,
  --   returns the corresponding predictor.

  function Create ( nbr : integer32;
                    sol : QuadDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32;
                    kind : Predictor_Type ) return Predictor_Array;

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
                s : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sets the leading coefficients of the data in the predictor p
  --   with the values in the vector s. 
  --   Resets all higher order coefficients to zero.

  procedure Newton_Fabry_Report 
              ( file : in file_type;
                nbrit : in integer32; absdx : in quad_double;
                fail : in boolean;
                z : in QuadDobl_Complex_Numbers.Complex_Number;
                r,err,step : in quad_double;
                numcff,dencff : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean );

  -- DESCRIPTION :
  --   Writes the results of the Newton-Fabry-Pade predictor,
  --   in quad double precision.

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

  procedure Newton_Fabry
              ( hom : in Link_to_System; prd : in Link_to_LU_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out quad_double;
                fail : out boolean; z : out Complex_Number;
                rad,err : out quad_double );
  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Link_to_System; prd : in Link_to_LU_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out quad_double;
                fail : out boolean; z : out Complex_Number;
                rad,err : out quad_double; output : in boolean );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with LU,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in quad double precision.

  -- ON ENTRY :
  --   file     to write data to if output (optional);
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   output   flag to indicate data output during computations,
  --            if a file is given on input.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the last corrector step;
  --   fail     indicates failure status;
  --   z        closest singularity estimated by Fabry's theorem;
  --   rad      estimates radius of convergence of the series;
  --   err      error estimate on the location of z.

  procedure Newton_Fabry
              ( hom : in Link_to_System; prd : in Link_to_SVD_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out quad_double;
                fail : out boolean; z : out Complex_Number;
                rad,err : out quad_double );
  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Link_to_System; prd : in Link_to_SVD_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out quad_double;
                fail : out boolean; z : out Complex_Number;
                rad,err : out quad_double; output : in boolean );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in quad double precision.

  -- ON ENTRY :
  --   file     to write data to if output (optional);
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   output   flag to indicate data output during computations,
  --            if a file is given on input.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the last corrector step;
  --   rcond    estimate for the inverse of the condition number;
  --   fail     indicates failure status;
  --   z        closest singularity estimated by Fabry's theorem;
  --   rad      estimates radius of convergence of the series;
  --   err      error estimate on the location of z.

  procedure Second
              ( hom : in Link_to_System; svh : in Link_to_SVD_Hessians;
                sol : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the Hessians for the convolution circuits in hom
  --   at the solution sol and stores the results in svh.

  function Distance ( svh : in Link_to_SVD_Hessians ) return quad_double;

  -- DESCRIPTION :
  --   Returns the estimate to the distance to the nearest solution,
  --   based on the singular values in svh.

  procedure Hesse_Pade
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in QuadDobl_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in QuadDobl_Complex_Vectors.Vector;
                res : out QuadDobl_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out quad_double );
  procedure Hesse_Pade
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in QuadDobl_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in QuadDobl_Complex_Vectors.Vector;
                res : out QuadDobl_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out quad_double;
                verbose : in boolean := true );

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
  --            if a file is given on input.

  -- ON RETURN :
  --   res      solution error estimated by Pade approximants in prd,
  --            the range of res is hom.crc'range;
  --   nrm      2-norm of res;
  --   step     computed curvature step.

  procedure Hesse_Pade
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in QuadDobl_Complex_Vectors.Vector;
                res : out QuadDobl_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out quad_double );
  procedure Hesse_Pade
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in QuadDobl_Complex_Vectors.Vector;
                res : out QuadDobl_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out quad_double;
                verbose : in boolean := true );

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
  --            if a file is given on input.

  -- ON RETURN :
  --   res      solution error estimated by Pade approximants in prd,
  --            the range of res is hom.crc'range;
  --   nrm      2-norm of res;
  --   step     computed curvature step.

  procedure Predictor_Feedback
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out Predictor_Vectors;
                numcff,dencff : in QuadDobl_Complex_VecVecs.VecVec;
                step : in out quad_double; minstep,alpha : in double_float;
                nrm,mixres : out quad_double; nbfail : out integer32 );
  procedure Predictor_Feedback
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out Predictor_Vectors;
                numcff,dencff : in QuadDobl_Complex_VecVecs.VecVec;
                step : in out quad_double; minstep,alpha : in double_float;
                nrm,mixres : out quad_double; nbfail : out integer32;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Given a step size, runs a predictor feedback loop,
  --   cutting the step size each time in half,
  --   until the mixed residual is less than the tolerance,
  --   or until the step becomes smaller than the minimum step size.

  -- ON ENTRY :
  --   file     to write output to if verbose;
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coeffiecients, for mixed residuals;
  --   psv      work space for solution vectors and residuals;
  --   numcff   coefficients of the numerator of the Pade approximants;
  --   dencff   coefficients of the denominator of the Pade approximants;
  --   step     current step size;
  --   minstep  the minimum step size;
  --   alpha    tolerance on the mixed predictor residual;
  --   eva      work space of range 1..hom.dim for Pade evaluation;
  --   radsol   work space for the radii of eva;
  --   res      work space of hom.crc'range for evaluation of eva;
  --   absres   work space of abh.crc'range for evaluation of radsol;
  --   verbose  flag for extra output if true.

  -- ON RETURN :
  --   step     shorter step size if nbfail > 0;
  --   psv      psv.sol contains the predicted solution,
  --            psv.radsol has the radii of the psv.sol components,
  --            psv.res is the residual of psv.sol, and
  --            psv.radres contains the evaluation at psv.radsol;
  --   nrm      max norm of the components in res;
  --   mixres   mixed residual.

-- MAIN PREDICTOR PROCEDURES :

  procedure LU_Prediction
              ( hom,abh : in Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out quad_double;
                fail : out boolean; step : out quad_double;
                nbpole,nbhess,nbmaxm : in out natural32 );
  procedure LU_Prediction
              ( file : in file_type; hom,abh : in Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out quad_double;
                fail : out boolean; step : out quad_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution.

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

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
  --            if a file is given on input.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   psv      psv.sol contains the predicted solution,
  --            psv.radsol has the radii of the psv.sol components,
  --            psv.res is the residual of psv.sol, and
  --            psv.radres contains the evaluation at psv.radsol;
  --   fail     indicates failure status;
  --   step     the step size;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the curve step was minimal;
  --   nbmaxm   updated number of times the maximum step size was minimal.

  procedure SVD_Prediction
              ( hom,abh : in Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out quad_double;
                fail : out boolean; step : out quad_double;
                nbpole,nbhess,nbmaxm : in out natural32 );
  procedure SVD_Prediction
              ( file : in file_type; hom,abh : in Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out quad_double;
                fail : out boolean; step : out quad_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false; verbose : in boolean := false );

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
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate extra output during computations,
  --            if a file is given on input;
  --   verbose  flag for intermediate numerical output,
  --            if a file is given on input.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   psv      psv.sol contains the predicted solution,
  --            psv.radsol has the radii of the psv.sol components,
  --            psv.res is the residual of psv.sol, and
  --            psv.radres contains the evaluation at psv.radsol;
  --   fail     indicates failure status;
  --   step     the step size;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the curve step was minimal;
  --   nbmaxm   updated number of times the maximum step size was minimal.

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

  procedure Clear ( h : in out Link_to_SVD_Hessians );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by h.

  procedure Clear ( p : in out LU_Predictor_Array );
  procedure Clear ( p : in out SVD_Predictor_Array );
  procedure Clear ( p : in out Predictor_Vectors_Array );
  procedure Clear ( h : in out SVD_Hessians_Array );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the array.

end QuadDobl_Predictor_Convolutions;
