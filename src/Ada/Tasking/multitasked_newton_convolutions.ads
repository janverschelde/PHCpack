with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Multitasked_Newton_Convolutions is

-- DESCRIPTION :
--   Runs Newton's method on power series
--   with the reverse mode of algorithmic differentation
--   and linearization to solve the matrix series equations,
--   in double, double double, and quad double arithmetic,
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
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Does one step with Newton's method using the LU factorization with
  --   multitasking in double, double double, or quad double precision.

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
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Does one step with Newton's method using the LU factorization with
  --   multitasking in double, double double, or quad double precision.

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

  -- DESCRIPTION :
  --   Applies several Newton steps using the LU factorization with
  --   multitasking in double, double double, or quad double precision.

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

  -- DESCRIPTION :
  --   Applies several Newton steps using the LU factorization with
  --   multitasking in double, double double, or quad double precision.

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
