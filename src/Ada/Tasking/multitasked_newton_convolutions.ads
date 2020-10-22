with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with TripDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with PentDobl_Complex_VecVecs;
with OctoDobl_Complex_VecVecs;
with DecaDobl_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with TripDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with PentDobl_Speelpenning_Convolutions;
with OctoDobl_Speelpenning_Convolutions;
with DecaDobl_Speelpenning_Convolutions;

package Multitasked_Newton_Convolutions is

-- DESCRIPTION :
--   Runs Newton's method on power series
--   with the reverse mode of algorithmic differentation
--   and linearization to solve the matrix series equations,
--   in double, double double, triple double, quad double,
--   penta double, octo double, and deca double arithmetic,
--   with multitasking for shared memory parallel computers.

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                absdx : out triple_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                absdx : out penta_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                absdx : out octo_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                absdx : out deca_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Does one step with Newton's method using the LU factorization with
  --   multitasking in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        a system of convolution circuits with a parameter;
  --   x        vector of coefficients of power series;
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range 1..dim, dim = A(0)'last(1);
  --   output   if true, then the multitasked procedures write to screen,
  --            otherwise, the computations remain silent.

  -- ON RETURN :
  --   x        updated coefficients of the series solution;
  --   absdx    the absolute value of the maximal component update dx;
  --   info     info from the LU factorization;
  --   ipvt     pivoting information about the LU factorization.

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double; rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                absdx : out triple_double; rcond : out triple_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                absdx : out penta_double; rcond : out penta_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                absdx : out octo_double; rcond : out octo_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                absdx : out deca_double; rcond : out deca_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Does one step with Newton's method using the LU factorization with
  --   multitasking in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        a system of convolution circuits with a parameter;
  --   x        vector of coefficients of power series;
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range 1..dim, dim = A(0)'last(1);
  --   output   if true, then the multitasked procedures write to screen,
  --            otherwise, the computations remain silent.

  -- ON RETURN :
  --   x        updated coefficients of the series solution;
  --   absdx    the absolute value of the maximal component update dx;
  --   rcond    estimate for the inverse of the condition number;
  --   ipvt     pivoting information about the LU factorization.

-- SEVERAL NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_float; absdx : out double_float; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_float; absdx : out double_float; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_double; absdx : out double_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_double; absdx : out double_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in triple_double; absdx : out triple_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in triple_double; absdx : out triple_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in quad_double; absdx : out quad_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in quad_double; absdx : out quad_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in penta_double; absdx : out penta_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in penta_double; absdx : out penta_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in octo_double; absdx : out octo_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in octo_double; absdx : out octo_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in deca_double; absdx : out deca_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in deca_double; absdx : out deca_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Applies several Newton steps using the LU factorization with
  --   multitasking in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

  -- ON ENTRY :
  --   file     if provided, then info, absdx is written in every step;
  --   nbt      the number of tasks;
  --   s        a system of convolution circuits with a parameter;
  --   x        vector of coefficients of power series;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, using as stop criterium;
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range 1..dim, dim = A(0)'last(1);
  --   output   if true, then the multitasked procedures write to screen,
  --            otherwise, the computations remain silent.

  -- ON RETURN :
  --   x        updated coefficients of the series solution;
  --   nbrit    number of iterations done;
  --   absdx    the absolute value of the maximal component update dx;
  --   fail     true if absdx > tol after nbrit iterations;
  --   info     info from the LU factorization;
  --   ipvt     pivoting information about the LU factorization.

-- SEVERAL NEWTON STEPS WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_float; absdx : out double_float; 
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_float; absdx : out double_float; 
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_double; absdx : out double_double; 
                fail : out boolean; rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_double; absdx : out double_double; 
                fail : out boolean; rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in triple_double; absdx : out triple_double; 
                fail : out boolean; rcond : out triple_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in triple_double; absdx : out triple_double; 
                fail : out boolean; rcond : out triple_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in quad_double; absdx : out quad_double; 
                fail : out boolean; rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in quad_double; absdx : out quad_double; 
                fail : out boolean; rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in penta_double; absdx : out penta_double; 
                fail : out boolean; rcond : out penta_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in penta_double; absdx : out penta_double; 
                fail : out boolean; rcond : out penta_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in octo_double; absdx : out octo_double; 
                fail : out boolean; rcond : out octo_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in octo_double; absdx : out octo_double; 
                fail : out boolean; rcond : out octo_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in deca_double; absdx : out deca_double; 
                fail : out boolean; rcond : out deca_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );
  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in deca_double; absdx : out deca_double; 
                fail : out boolean; rcond : out deca_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Applies several Newton steps using the LU factorization with
  --   multitasking in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

  -- ON ENTRY :
  --   file     if provided, then info, absdx is written in every step;
  --   nbt      the number of tasks;
  --   s        a system of convolution circuits with a parameter;
  --   x        vector of coefficients of power series;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, using as stop criterium;
  --   wrk      work space as a vector of vectors of range 1..nbt,
  --            with every vector of range 1..dim, dim = A(0)'last(1);
  --   output   if true, then the multitasked procedures write to screen,
  --            otherwise, the computations remain silent.

  -- ON RETURN :
  --   x        updated coefficients of the series solution;
  --   nbrit    number of iterations done;
  --   absdx    the absolute value of the maximal component update dx;
  --   fail     true if absdx > tol after nbrit iterations;
  --   rcond    estimate for the inverse condition number;
  --   ipvt     pivoting information about the LU factorization.

end Multitasked_Newton_Convolutions;
