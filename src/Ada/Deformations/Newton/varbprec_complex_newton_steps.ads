with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with Standard_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_JacoMats;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_JacoMats;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_JacoMats;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Jaco_Matrices;
with Multprec_Complex_Laur_Systems;
with Multprec_Complex_Laur_JacoMats;

package Varbprec_Complex_Newton_Steps is

-- DESCRIPTION :
--   Applies one step of Newton's method to a polynomial system,
--   with condition number estimation.  
--   We have versions for four different levels of precision:
--   standard double, double double, quad double, and arbitrary 
--   multiprecision.  Variable precision methods take string
--   representations on input to evaluate the coefficients and
--   coordinates at sufficiently high precision to meet the
--   loss of accuracy computed by condition numbers.
--   We can run Newton steps on regular polynomial systems and on
--   systems of Laurent polynomials which may have negative exponents.

-- PART I : estimators and Newton steps at various level of precision

  procedure Estimate_Loss_in_Newton_Step
              ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
                z : in Standard_Complex_Vectors.Vector;
                jfz : out Standard_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Standard_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_float;
                jfzloss,fzloss : out integer32 );
  procedure Estimate_Loss_in_Newton_Step
              ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
                z : in DoblDobl_Complex_Vectors.Vector;
                jfz : out DoblDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out DoblDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_double;
                jfzloss,fzloss : out integer32 );
  procedure Estimate_Loss_in_Newton_Step
              ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
                z : in QuadDobl_Complex_Vectors.Vector;
                jfz : out QuadDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out QuadDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out quad_double;
                jfzloss,fzloss : out integer32 );
  procedure Estimate_Loss_in_Newton_Step
              ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
                jf : in Multprec_Complex_Jaco_Matrices.Jaco_Mat;
                z : in Multprec_Complex_Vectors.Vector;
                jfz : out Multprec_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Multprec_Complex_Vectors.Vector;
                jfzrco,fzrco : out Floating_Number;
                jfzloss,fzloss : out integer32 );
                
  -- DESCRIPTION :
  --   Estimates the condition numbers and determines the loss of decimal
  --   places when applying one Newton step on a system at some point.

  -- REQUIRED : the system has as many equations as unknowns,
  --   in particular: f'range = z'range = jf'range(1) = jf'range(2).

  -- ON ENTRY :
  --   f        a polynomial system in as many variables as equations;
  --   jf       the Jacobian matrix of f;
  --   z        current approximation for a solution of f(x) = 0.

  -- ON RETURN :
  --   jfz      output of lufco on the Jacobian matrix jf,
  --            evaluated at z, suitable for back substitution if nonsingular;
  --   piv      pivoting information computed by lufco;
  --   fz       evaluation of f at z;
  --   jfzrco   estimate for inverse condition number of jfz;
  --   fzrco    inverse condition number of evaluating f at z;
  --   jfzloss  10-logarithm of jfzrco, as indication for the loss of
  --            decimal places when computing the update in a Newton step;
  --   fzloss   10-logarithm of fzrco, as indication for the loss of
  --            decimal places when evaluating f at z, if absfz > 0.1,
  --            otherwise the 10-logarithm of the denominator of fzrco
  --            is used to predict the loss of decimal places.

  procedure Estimate_Loss_in_Newton_Step
              ( f : in Standard_Complex_Laur_Systems.Laur_Sys;
                jf : in Standard_Complex_Laur_JacoMats.Jaco_Mat;
                z : in Standard_Complex_Vectors.Vector;
                jfz : out Standard_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Standard_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_float;
                jfzloss,fzloss : out integer32 );
  procedure Estimate_Loss_in_Newton_Step
              ( f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                jf : in DoblDobl_Complex_Laur_JacoMats.Jaco_Mat;
                z : in DoblDobl_Complex_Vectors.Vector;
                jfz : out DoblDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out DoblDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_double;
                jfzloss,fzloss : out integer32 );
  procedure Estimate_Loss_in_Newton_Step
              ( f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Jaco_Mat;
                z : in QuadDobl_Complex_Vectors.Vector;
                jfz : out QuadDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out QuadDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out quad_double;
                jfzloss,fzloss : out integer32 );
  procedure Estimate_Loss_in_Newton_Step
              ( f : in Multprec_Complex_Laur_Systems.Laur_Sys;
                jf : in Multprec_Complex_Laur_JacoMats.Jaco_Mat;
                z : in Multprec_Complex_Vectors.Vector;
                jfz : out Multprec_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Multprec_Complex_Vectors.Vector;
                jfzrco,fzrco : out Floating_Number;
                jfzloss,fzloss : out integer32 );
                
  -- DESCRIPTION :
  --   Estimates the condition numbers and determines the loss of decimal
  --   places when applying one Newton step on a system at some point.

  -- REQUIRED : the system has as many equations as unknowns,
  --   in particular: f'range = z'range = jf'range(1) = jf'range(2).

  -- ON ENTRY :
  --   f        Laurent polynomial system in as many variables as equations;
  --   jf       the Jacobian matrix of f;
  --   z        current approximation for a solution of f(x) = 0.

  -- ON RETURN :
  --   jfz      output of lufco on the Jacobian matrix jf,
  --            evaluated at z, suitable for back substitution if nonsingular;
  --   piv      pivoting information computed by lufco;
  --   fz       evaluation of f at z;
  --   jfzrco   estimate for inverse condition number of jfz;
  --   fzrco    inverse condition number of evaluating f at z;
  --   jfzloss  10-logarithm of jfzrco, as indication for the loss of
  --            decimal places when computing the update in a Newton step;
  --   fzloss   10-logarithm of fzrco, as indication for the loss of
  --            decimal places when evaluating f at z, if absfz > 0.1,
  --            otherwise the 10-logarithm of the denominator of fzrco
  --            is used to predict the loss of decimal places.

  procedure do_Newton_Step
              ( z : in out Standard_Complex_Vectors.Vector;
                jfz : in Standard_Complex_Matrices.Matrix;
                piv : in Standard_Integer_Vectors.Vector;
                fz : in Standard_Complex_Vectors.Vector;
                err : out double_float );
  procedure do_Newton_Step
              ( z : in out DoblDobl_Complex_Vectors.Vector;
                jfz : in DoblDobl_Complex_Matrices.Matrix;
                piv : in Standard_Integer_Vectors.Vector;
                fz : in DoblDobl_Complex_Vectors.Vector;
                err : out double_double );
  procedure do_Newton_Step
              ( z : in out QuadDobl_Complex_Vectors.Vector;
                jfz : in QuadDobl_Complex_Matrices.Matrix;
                piv : in Standard_Integer_Vectors.Vector;
                fz : in QuadDobl_Complex_Vectors.Vector;
                err : out quad_double );
  procedure do_Newton_Step
              ( z : in out Multprec_Complex_Vectors.Vector;
                jfz : in Multprec_Complex_Matrices.Matrix;
                piv : in Standard_Integer_Vectors.Vector;
                fz : in Multprec_Complex_Vectors.Vector;
                err : out Floating_Number );

  -- DESCRIPTION :
  --   Updates the point z with one Newton step,
  --   using the output of Estimate_Loss_in_Newton_Step.

  -- ON ENTRY :
  --   z        current approximation for a solution;
  --   jfz      output of lufco, LU factorization of Jacobian matrix;
  --   piv      pivoting information of the LU factorization;
  --   fz       system evaluated at z.

  -- ON RETURN :
  --   z        updated approximation for a solution;
  --   err      magnitude of the correction term added to z.

  function Minimum ( a,b : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the minimum of a and b, used to compute the total loss,
  --   as the minimum of the loss in the accuracy of the linear solving
  --   and the loss in the accuracy of the polynomial evaluation.

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
                z : in out Standard_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_float;
                fz : out Standard_Complex_Vectors.Vector;
                fail : out boolean );
  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
                z : in out DoblDobl_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_double;
                fz : out DoblDobl_Complex_Vectors.Vector;
                fail : out boolean );
  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
                z : in out QuadDobl_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out quad_double;
                fz : out QuadDobl_Complex_Vectors.Vector;
                fail : out boolean );
  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
                jf : in Multprec_Complex_Jaco_Matrices.Jaco_Mat;
                z : in out Multprec_Complex_Vectors.Vector;
                prec,want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out Floating_Number;
                fz : out Multprec_Complex_Vectors.Vector;
                fail : out boolean );

  -- DESCRIPTION :
  --   Estimates the loss of accuracy and if the loss is small enough
  --   to achieve the desired wanted number of decimal places with
  --   the double, double double, quad double, or arbitrary precision, 
  --   then one Newton step is done.

  -- ON ENTRY :
  --   f        polynomial system in as many variables as equations;
  --   jf       Jacobian matrix of the system f;
  --   z        current approximation for a solution of f;
  --   prec     number of decimal places in the arbitrary multiprecision;
  --   want     wanted number of decimal places to be accurate.

  -- ON RETURN :
  --   z        if not fail, then updated approximation for a solution;
  --   loss     estimated loss of decimal places;
  --   jfzrco   inverse condition number of the Jacobian matrix at z;
  --   fzrco    inverse condition number of the evaluation problem of f at z;
  --   err      if not fail, then magnitude of the correction added to z;
  --   fz       system f evaluated at z;
  --   fail     if the precision is insufficient to meet the accuracy.

  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in Standard_Complex_Laur_Systems.Laur_Sys;
                jf : in Standard_Complex_Laur_JacoMats.Jaco_Mat;
                z : in out Standard_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_float;
                fz : out Standard_Complex_Vectors.Vector;
                fail : out boolean );
  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                jf : in DoblDobl_Complex_Laur_JacoMats.Jaco_Mat;
                z : in out DoblDobl_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_double;
                fz : out DoblDobl_Complex_Vectors.Vector;
                fail : out boolean );
  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Jaco_Mat;
                z : in out QuadDobl_Complex_Vectors.Vector;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out quad_double;
                fz : out QuadDobl_Complex_Vectors.Vector;
                fail : out boolean );
  procedure Newton_Step_to_Wanted_Accuracy
              ( f : in Multprec_Complex_Laur_Systems.Laur_Sys;
                jf : in Multprec_Complex_Laur_JacoMats.Jaco_Mat;
                z : in out Multprec_Complex_Vectors.Vector;
                prec,want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out Floating_Number;
                fz : out Multprec_Complex_Vectors.Vector;
                fail : out boolean );

  -- DESCRIPTION :
  --   Estimates the loss of accuracy and if the loss is small enough
  --   to achieve the desired wanted number of decimal places with
  --   the double, double double, quad double, or arbitrary precision, 
  --   then one Newton step is done.

  -- ON ENTRY :
  --   f        Laurent polynomial system in as many variables as equations;
  --   jf       Jacobian matrix of the system f;
  --   z        current approximation for a solution of f;
  --   prec     number of decimal places in the arbitrary multiprecision;
  --   want     wanted number of decimal places to be accurate.

  -- ON RETURN :
  --   z        if not fail, then updated approximation for a solution;
  --   loss     estimated loss of decimal places;
  --   jfzrco   inverse condition number of the Jacobian matrix at z;
  --   fzrco    inverse condition number of the evaluation problem of f at z;
  --   err      if not fail, then magnitude of the correction added to z;
  --   fz       system f evaluated at z;
  --   fail     if the precision is insufficient to meet the accuracy.

-- PART II : estimating loss of accuracy

  procedure Standard_Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out double_float; loss : out integer32 );
  procedure DoblDobl_Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out double_double; loss : out integer32 );
  procedure QuadDobl_Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out quad_double; loss : out integer32 );
  procedure Multprec_Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                prec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 );

  -- DESCRIPTION :
  --   Evaluates the system in f and vector in z in standard
  --   double, double double, or quad double  precision to determine the 
  --   condition numbers of the Jacobian matrix and evaluation problem.

  -- REQUIRED : f'range = number of equations and variables,
  --   so z will parse into a complex vector of that many variables.

  -- ON ENTRY :
  --   f        array of strings the contain a polynomial system;
  --   z        string of numbers, separated by newline symbols.

  -- ON RETURN :
  --   jfrco    inverse condition number of Jacobian matrix of f at z;
  --   fzrco    inverse condition number of evaluating f at z;
  --   loss     loss of number of decimal places, as negative number.

  procedure Standard_Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out double_float; loss : out integer32 );
  procedure DoblDobl_Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out double_double; loss : out integer32 );
  procedure QuadDobl_Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                jfrco,fzrco : out quad_double; loss : out integer32 );
  procedure Multprec_Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                prec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 );

  -- DESCRIPTION :
  --   Evaluates the system in f and vector in z in standard
  --   double, double double, or quad double  precision to determine the 
  --   condition numbers of the Jacobian matrix and evaluation problem.

  -- REQUIRED : f'range = number of equations and variables,
  --   so z will parse into a complex vector of that many variables.

  -- ON ENTRY :
  --   f        array of strings the contain a Laurent polynomial system;
  --   z        string of numbers, separated by newline symbols.

  -- ON RETURN :
  --   jfrco    inverse condition number of Jacobian matrix of f at z;
  --   fzrco    inverse condition number of evaluating f at z;
  --   loss     loss of number of decimal places, as negative number.

  procedure Estimate_Loss_for_Polynomial_System
              ( f : in Array_of_Strings; z : in string;
                maxprec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 );
  procedure Estimate_Loss_for_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in string;
                maxprec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 );

  -- DESCRIPTION :
  --   Evaluates the system in f and vector in z at various precision
  --   until sufficiently high to determine the condition numbers
  --   of the Jacobian matrix and evaluation problem.

  -- ON ENTRY :
  --   f        array of strings the contain a (Laurent) polynomial system;
  --   z        string of numbers, separated by newline symbols;
  --   maxprec  maximum number of decimal places that could be used.

  -- ON RETURN :
  --   jfrco    inverse condition number of Jacobian matrix of f at z;
  --   fzrco    inverse condition number of evaluating f at z;
  --   loss     loss of number of decimal places, as negative number.

  function Estimate_Loss_for_Polynomial_System
              ( f : Array_of_Strings; z : string; maxprec : natural32 )
              return integer32;
  function Estimate_Loss_for_Laurent_Polynomials
              ( f : Array_of_Strings; z : string; maxprec : natural32 )
              return integer32;

  -- DESCRIPTION :
  --   Returns estimated loss of decimal places as a negative number
  --   based on condition number estimates for one Newton step on f at z.
  --   This function encapsulates the same named procedure.

-- PART III : computing after parsing to the precision

  procedure Standard_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float );
  procedure DoblDobl_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float );
  procedure QuadDobl_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float );
  procedure Multprec_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                prcn : in natural32; err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Evaluates the system f and initial approximation z in standard
  --   double, double double, quad double, or arbitrary multiprecision.
  --   One Newton step is done on the evaluated system.

  -- ON ENTRY :
  --   f        string representation of a polynomial system;
  --   z        string representation of an initial approximation;
  --   prcn     number of decimal places in arbitrary multiprecision.

  -- ON RETURN :
  --   z        updated approximation for a solution of f;
  --   err      magnitude of the correction added to z;
  --   rco      estimate of the inverse condition number of Jacobian at z;
  --   res      magnitude of f evaluated at z.

  procedure Standard_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float );
  procedure DoblDobl_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float );
  procedure QuadDobl_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                err,rco,res : out double_float );
  procedure Multprec_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                prcn : in natural32; err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Evaluates the system f and initial approximation z in standard
  --   double, double double, quad double, or arbitrary multiprecision.
  --   One Newton step is done on the evaluated system.

  -- ON ENTRY :
  --   f        string representation of a Laurent polynomial system;
  --   z        string representation of an initial approximation;
  --   prcn     number of decimal places in arbitrary multiprecision.

  -- ON RETURN :
  --   z        updated approximation for a solution of f;
  --   err      magnitude of the correction added to z;
  --   rco      estimate of the inverse condition number of Jacobian at z;
  --   res      magnitude of f evaluated at z.

  procedure do_Newton_Step_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                loss,want : in integer32; err,rco,res : out double_float );
  procedure do_Newton_Step_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                loss,want : in integer32; err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Selects the precision to meet the wanted number of accurate
  --   decimal places taking the loss of decimal places into account
  --   when performing one Newton step on f at z.

  -- ON ENTRY :
  --   f        string representation of a (Laurent) polynomial system;
  --   z        string representation of an initial approximation;
  --   loss     estimated loss of decimal places as a negative number;
  --   want     wanted number of accurate decimal places.

  -- ON RETURN :
  --   z        updated approximation for a solution of f;
  --   err      magnitude of the correction added to z;
  --   rco      estimate of the inverse condition number of Jacobian at z;
  --   res      magnitude of f evaluated at z.

-- PART IV : sequence of Newton steps to the wanted accuracy

  procedure Newton_Steps_on_Polynomial_System
              ( f : in Array_of_Strings; z : in out Link_to_String;
                want : in integer32; maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float );
  procedure Newton_Steps_on_Laurent_Polynomials
              ( f : in Array_of_Strings; z : in out Link_to_String;
                want : in integer32; maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float );
  procedure Newton_Steps_on_Polynomial_System
              ( file : in file_type;
                f : in Array_of_Strings; z : in out Link_to_String;
                want : in integer32; maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float );
  procedure Newton_Steps_on_Laurent_Polynomials
              ( file : in file_type;
                f : in Array_of_Strings; z : in out Link_to_String;
                want : in integer32; maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Performs a sequence of Newton steps to approximate a root
  --   of a polynomial system up to a wanted number of accuracte
  --   decimal places taking the loss of decimal places into account
  --   when performing one Newton step on f at z.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted,
  --            then the procedure will be silent;
  --   f        string representation of a (Laurent) polynomial system;
  --   z        string representation of an initial approximation;
  --   want     wanted number of accurate decimal places;
  --   maxprc   maximum number of decimal places that can be used
  --            to estimate the loss of accuracy;
  --   maxitr   maximum number of Newton steps.

  -- ON RETURN :
  --   z        updated approximation for a solution of f;
  --   loss     estimated loss of decimal places as a negative number;
  --   err      magnitude of the correction added to z;
  --   rco      estimate of the inverse condition number of Jacobian at z;
  --   res      magnitude of f evaluated at z.

end Varbprec_Complex_Newton_Steps;
