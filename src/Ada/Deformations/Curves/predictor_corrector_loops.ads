with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Homotopy_Continuation_Parameters;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;

package Predictor_Corrector_Loops is

-- DESCRIPTION :
--   A predictor-corrector loop does one predictor-corrector step,
--   with two feedback loops, one in the predictor, one in the corrector,
--   in case the predicted step size would be too large.

  procedure Predictor_Corrector_Loop
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_float;
                step,mixres : out double_float; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32; fail : out boolean );
  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_float;
                step,mixres : out double_float; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true );
  procedure Predictor_Corrector_Loop
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_double;
                step,mixres : out double_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32; fail : out boolean );
  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_double;
                step,mixres : out double_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true );
  procedure Predictor_Corrector_Loop
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out quad_double;
                step,mixres : out quad_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32; fail : out boolean );
  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out quad_double;
                step,mixres : out quad_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true );

  -- DESCRIPTION :
  --   Does one predictor-corrector step on one solution,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write the extra output to (optional);
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   homlead  leading coefficients for the circuits in hom;
  --   abhlead  leading coefficients for the circuits in abh;
  --   pars     values for the tolerances and parameters,
  --            if pars.corsteps = 0, then there is no corrector step;
  --   maxit    maximum number of steps in Newton's method on power series;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   prd      work space for the Newton-Fabry-Pade predictor;
  --   psv      work space vectors for the predictor,
  --            psv.sol contains a start solution;
  --   svh      work space for Hessian convolutions;
  --   endt     the end value for the homotopy continuation parameter t;
  --   acct     accumulated sum of all successful steps, equals the
  --            current value of the homotopy continuation parameter t;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the Hessian step was minimal;
  --   nbmaxm   number of times the maximum step was minimal;
  --   verbose  flag for extra output, if a file is given on input.

  -- ON RETURN :
  --   psv.sol  the corrected solution;
  --   dx       last update to the solution;
  --   ipvt     pivoting information for the LU Newton steps;
  --   acct     updated value for the homotopy continuation parameter t;
  --   step     the step size;
  --   mixres   mixed residual of corrector loop if pars.corsteps < 0;
  --   nbrit    number of iterations in the corrector loop;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the Hessian step was minimal;
  --   nbmaxm   updated number of times the maximum step was minimal;
  --   fail     true if the prescribed tolerance was not reached,
  --            false otherwise.

  procedure Track_One_Path
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out double_float;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float; fail : out boolean );
  procedure Track_One_Path
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out double_float;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean; verbose : in boolean := true );
  procedure Track_One_Path
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out double_double;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float; fail : out boolean );
  procedure Track_One_Path
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out double_double;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean; verbose : in boolean := true );
  procedure Track_One_Path
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out quad_double;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float; fail : out boolean );
  procedure Track_One_Path
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out quad_double;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean; verbose : in boolean := true );

  -- DESCRIPTION :
  --   Tracks one path in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write the extra output to, if verbose (optional);
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   homlead  leading coefficients for the circuits in hom;
  --   abhlead  leading coefficients for the circuits in abh;
  --   pars     values for the tolerances and parameters;
  --   maxit    maximum number of steps in Newton's method on power series;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   prd      work space for the Newton-Fabry-Pade predictor;
  --   psv      work space vectors for the predictor,
  --            psv.sol contains a start solution;
  --   svh      work space for Hessian convolutions;
  --   wrk      work space vector for power series coefficients
  --            during the shifting of the coefficients;
  --   acct     start value for the homotopy continuation parameter t;
  --   verbose  indicates if extra output is requested,
  --            if a file is given on input.

  -- ON RETURN :
  --   psv.sol  the corrected solution;
  --   dx       last update to the solution;
  --   ipvt     pivoting information for the LU Newton steps;
  --   acct     accumulated value of the homotopy continuation parameter t;
  --   mixres   mixed residual of corrector loop if pars.corsteps > 0;
  --   tnbrit   total number of corrector iterations;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the Hessian step was minimal;
  --   nbmaxm   updated number of times the maximum step was minimal;
  --   nbsteps  number of predictor-corrector steps done;
  --   minstpz  minimum step size;
  --   maxstpz  maximum step size;
  --   fail     true if the prescribed tolerance was not reached,
  --            false otherwise.

  procedure Track_All_Paths
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true );
  procedure Track_All_Paths
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true );
  procedure Track_All_Paths
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Tracks all paths starting at the solutions in sols,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write output information to;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

  procedure Write_Path_Statistics
              ( file : in file_type;
                nbrit,nbpole,nbhess,nbmaxm,nbsteps : in natural32;
                minstpz,maxstpz : in double_float );

  -- DESCRIPTION :
  --   Writes some counts on the path tracker to file.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbrit    the number of corrector iterations on the path;
  --   nbpole   number of times pole step was minimal;
  --   nbhess   number of times curvature step was minimal;
  --   nbmaxm   number of times maximum step was minimal;
  --   nbstep   total number of steps on the path;
  --   minstpz  minimum step size;
  --   maxstpz  maximum step size.

  procedure Write_Total_Path_Statistics
              ( file : in file_type;
                minnbrsteps,maxnbrsteps : in natural32;
                mincorsteps,maxcorsteps : in natural32;
                ratpole,rathess,ratmaxm : in double_float;
                minstpz,maxstpz : in double_float );

  -- DESCRIPTION :
  --   Writes some overall statistics to file.

  -- ON ENTRY :
  --   file          must be opened for output;
  --   minnbrsteps   smallest number of steps on a path;
  --   maxnbrsteps   largest number of steps on a path;
  --   mincorsteps   smallest number of corrector iterations on a path;
  --   maxcorsteps   largest number of corrector iterations on a path;
  --   ratpole       ratio of times pole step was minimal;
  --   rathess       ratio of times curvature step was minimal;
  --   ratmaxm       ratio of times curvature step was minimal;
  --   minstpz       minimum step size;
  --   maxstpz       maximum step size.

end Predictor_Corrector_Loops;
