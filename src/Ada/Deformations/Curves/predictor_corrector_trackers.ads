with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with Standard_Coefficient_Circuits;
with Standard_Coefficient_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Homotopy_Continuation_Parameters;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;

package Predictor_Corrector_Trackers is

-- DESCRIPTION :
--   A predictor-corrector tracker runs a predictor-corrector loop,
--   for one or on all paths defined by a polynomial homotopy,
--   in standard double, double double, or quad double precision.

-- ON COEFFICIENT CONVOLUTION CIRCUITS :

  procedure Track_One_Path
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                pwt : in Standard_Floating_Vectors.Link_to_Vector;
                acct,mixres : in out double_float;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean; vrblvl : in integer32 := 0 );
  procedure Track_One_Path
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                pwt : in Standard_Floating_Vectors.Link_to_Vector;
                acct,mixres : in out double_float;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks one path in double precision.

  -- ON ENTRY :
  --   file     to write the extra output to, if verbose (optional);
  --   hom      system of homotopy convolution circuits;
  --   rcf      real parts of the coefficients, for shifting circuits;
  --   icf      imaginary parts of the coefficients, for shifting circuits;
  --   cfh      coefficient circuits for the homotopy;
  --   abh      circuits with radii of cfh as coefficients,
  --            to compute the mixed residuals in the corrector loop;
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
  --   rx       work space for the real parts of the series coefficients;
  --   ix       work space for the imaginary parts of the series coefficients;
  --   xr       work space allocated for the real parts of solution vectors;
  --   xi       work space allocated for the imag parts of solution vectors;
  --   vh       space allocated for dim matrices, dim = dimension,
  --            all matrices have 1..dim for range(1) and range(2);
  --   svls     svls(0) contains the singular values of s.jm, and
  --            svls(k) contains the singular values of vh(k),
  --            for k in vh'range.
  --   pwt      work space the powers of the value used in the shift,
  --            pwt'range = 0..deg, where deg is the degree of the series;
  --   acct     start value for the homotopy continuation parameter t;
  --   verbose  indicates if extra output is requested,
  --            if a file is given on input;
  --   vrblvl   the verbose level.

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
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure Track_All_Paths
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                minpastp,maxpastp,ratpole,rathess,ratmaxm : out double_float;
                mincorsteps,maxcorsteps : out natural32;
                minnbrsteps,maxnbrsteps : out natural32;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks all solution paths defined by the homotopy in hom,
  --   starting at the solutions in sols, without output.

  -- ON ENTRY :
  --   file     optional output file, if verbose;
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
  --   vrblvl   the verbose level.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths;
  --   minpastp is the minimum step size on a path;
  --   maxpastp is the maximum step size on a path;
  --   ratpole  is the ratio of times the pole step size was minimal;
  --   rathess  is the ratio of times the curvature step size was minimal;
  --   ratmaxm  is the ratio of times the maximal step bound was minimal;
  --   mincorsteps is the minimum number of corrector steps on a path;
  --   maxcorsteps is the maximum number of corrector steps on a path;
  --   minnbrsteps is the minimum number of steps on a path;
  --   maxnbrsteps is the maximum number of steps on a path.

  procedure Track_All_Paths
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks all solution paths defined by the homotopy in hom,
  --   starting at the solutions in sols, writes the solution and
  --   statistics at the end of each path to file.

  -- ON ENTRY :
  --   file     optional output file, if verbose;
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
  --   verbose  indicates if extra output is requested;
  --   vrblvl   the verbose level.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

-- ON COMPLEX CONVOLUTION CIRCUITS :

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
                minstpz,maxstpz : out double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
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
                fail : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
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
                minstpz,maxstpz : out double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
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
                fail : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
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
                minstpz,maxstpz : out double_float; fail : out boolean;
                vrblvl : in integer32 := 0 );
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
                fail : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

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
  --            if a file is given on input;
  --   vrblvl   the verbose level.

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
                verbose : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Track_All_Paths
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );
  procedure Track_All_Paths
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );

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
  --   verbose  indicates if extra output is requested;
  --   vrblvl   the verbose level.
  
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

end Predictor_Corrector_Trackers;
