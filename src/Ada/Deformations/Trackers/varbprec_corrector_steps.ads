with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;

package Varbprec_Corrector_Steps is

-- DESCRIPTION :
--   With condition number estimates at various precision levels,
--   we can predict the loss of accuracy and adjust the working precision to
--   correct the solution to the desired number of accurate decimal places.
--   The procedures in this package assume that Varbprec_Homotopy has
--   been initialized properly.

-- PART I : estimators and corrector steps at various level of precision

  procedure Estimate_Loss_in_Newton_Step
              ( z : in Standard_Complex_Vectors.Vector;
                t : in Standard_Complex_Numbers.Complex_Number;
                jfz : out Standard_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Standard_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_float;
                jfzloss,fzloss : out integer32 );
  procedure Estimate_Loss_in_Newton_Step
              ( z : in DoblDobl_Complex_Vectors.Vector;
                t : in DoblDobl_Complex_Numbers.Complex_Number;
                jfz : out DoblDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out DoblDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out double_double;
                jfzloss,fzloss : out integer32 );
  procedure Estimate_Loss_in_Newton_Step
              ( z : in QuadDobl_Complex_Vectors.Vector;
                t : in QuadDobl_Complex_Numbers.Complex_Number;
                jfz : out QuadDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out QuadDobl_Complex_Vectors.Vector;
                jfzrco,fzrco : out quad_double;
                jfzloss,fzloss : out integer32 );
  procedure Estimate_Loss_in_Newton_Step
              ( z : in Multprec_Complex_Vectors.Vector;
                t : in Multprec_Complex_Numbers.Complex_Number;
                deci : in natural32;
                jfz : out Multprec_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                fz : out Multprec_Complex_Vectors.Vector;
                jfzrco,fzrco : out Floating_Number; 
                jfzloss,fzloss : out integer32 );

  -- DESCRIPTION :
  --   Estimates the condition numbers and determines the loss of decimal
  --   places when applying one Newton step on a system at some point.

  -- REQUIRED : the package Varbprec_Homotopy has been initialized.

  -- ON ENTRY :
  --   z        current approximation for a solution of h(x,t) = 0,
  --            where h is defined by the package Varbprec_Homotopy;
  --   t        current value of the homotopy continuation parameter;
  --   deci     number of decimal places in working precision when
  --            arbitrary multiprecision complex arithmetic.

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

  procedure Newton_Step_to_Wanted_Accuracy
              ( z : in out Standard_Complex_Vectors.Vector;
                t : in Standard_Complex_Numbers.Complex_Number;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_float;
                fz : out Standard_Complex_Vectors.Vector;
                fail : out boolean );
  procedure Newton_Step_to_Wanted_Accuracy
              ( z : in out DoblDobl_Complex_Vectors.Vector;
                t : in DoblDobl_Complex_Numbers.Complex_Number;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out double_double;
                fz : out DoblDobl_Complex_Vectors.Vector;
                fail : out boolean );
  procedure Newton_Step_to_Wanted_Accuracy
              ( z : in out QuadDobl_Complex_Vectors.Vector;
                t : in QuadDobl_Complex_Numbers.Complex_Number;
                want : in integer32; loss : out integer32;
                jfzrco,fzrco,err : out quad_double;
                fz : out QuadDobl_Complex_Vectors.Vector;
                fail : out boolean );
  procedure Newton_Step_to_Wanted_Accuracy
              ( z : in out Multprec_Complex_Vectors.Vector;
                t : in Multprec_Complex_Numbers.Complex_Number;
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
  --   z        current approximation for a solution of h(x,t);
  --   t        current value of the homotopy continuation parameter;
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

  procedure Standard_Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in Standard_Complex_Numbers.Complex_Number;
                jfrco,fzrco : out double_float; loss : out integer32 );
  procedure DoblDobl_Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in DoblDobl_Complex_Numbers.Complex_Number;
                jfrco,fzrco : out double_double; loss : out integer32 );
  procedure QuadDobl_Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in QuadDobl_Complex_Numbers.Complex_Number;
                jfrco,fzrco : out quad_double; loss : out integer32 );
  procedure Multprec_Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in Multprec_Complex_Numbers.Complex_Number;
                prec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 );

  -- DESCRIPTION :
  --   Evaluates the homotopy at the vector in z in standard
  --   double, double double, or quad double  precision to determine the 
  --   condition numbers of the Jacobian matrix and evaluation problem.

  -- REQUIRED : the variable precision homotopy Varbprec_Homotopy 
  --   is initialized properly.

  -- ON ENTRY :
  --   z        string of numbers, separated by newline symbols;
  --   t        current value of the homotopy continuation parameter.

  -- ON RETURN :
  --   jfrco    inverse condition number of Jacobian matrix of f at z;
  --   fzrco    inverse condition number of evaluating f at z;
  --   loss     loss of number of decimal places, as negative number.

  procedure Estimate_Loss_for_Polynomial_Homotopy
              ( z : in string; t : in Standard_Complex_Numbers.Complex_Number;
                maxprec : in natural32;
                jfrco,fzrco : out Floating_Number; loss : out integer32 );

  -- DESCRIPTION :
  --   Evaluates the homotopy at vector in z at various precision
  --   until sufficiently high to determine the condition numbers
  --   of the Jacobian matrix and evaluation problem.

  -- ON ENTRY :
  --   z        string of numbers, separated by newline symbols;
  --   t        current value of the homotopy continuation parameter;
  --   maxprec  maximum number of decimal places that could be used.

  -- ON RETURN :
  --   jfrco    inverse condition number of Jacobian matrix of f at z;
  --   fzrco    inverse condition number of evaluating f at z;
  --   loss     loss of number of decimal places, as negative number.

  function Estimate_Loss_for_Polynomial_Homotopy
              ( z : string; t : Standard_Complex_Numbers.Complex_Number;
                maxprec : natural32 ) return integer32;

  -- DESCRIPTION :
  --   Returns estimated loss of decimal places as a negative number
  --   based on condition number estimates for one Newton step on f at z.
  --   This function encapsulates the same named procedure.

-- PART III : computing after parsing to the precision

  procedure Standard_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                err,rco,res : out double_float );
  procedure DoblDobl_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                err,rco,res : out double_float );
  procedure QuadDobl_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                err,rco,res : out double_float );
  procedure Multprec_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                prcn : in natural32; err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Evaluates the homotopy and the initial approximation z in standard
  --   double, double double, quad double, or arbitrary multiprecision.
  --   One Newton step is done on the evaluated system.

  -- ON ENTRY :
  --   z        string representation of an initial approximation;
  --   t        corresponding value of the homotopy continuation parameter;
  --   prcn     number of decimal places in arbitrary multiprecision.

  -- ON RETURN :
  --   z        updated approximation for a solution of the homotopy;
  --   err      magnitude of the correction added to z;
  --   rco      estimate of the inverse condition number of Jacobian at z;
  --   res      magnitude of homotopy evaluated at z and t.

  procedure do_Newton_Step_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                loss : in integer32; want : in natural32;
                err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Selects the precision to meet the wanted number of accurate
  --   decimal places taking the loss of decimal places into account
  --   when performing one Newton step on f at z.

  -- ON ENTRY :
  --   z        string representation of an initial approximation;
  --   t        correspondinf value of the homtopy continuation parameter;
  --   loss     estimated loss of decimal places as a negative number;
  --   want     wanted number of accurate decimal places.

  -- ON RETURN :
  --   z        updated approximation for a solution of a homotopy at t;
  --   err      magnitude of the correction added to z;
  --   rco      estimate of the inverse condition number of Jacobian at z;
  --   res      magnitude of the homotopy evaluated at z.

-- PART IV : sequence of Newton steps to the wanted accuracy

  procedure Newton_Steps_on_Polynomial_Homotopy
              ( z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                want,maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float );
  procedure Newton_Steps_on_Polynomial_Homotopy
              ( file : in file_type; z : in out Link_to_String;
                t : in Standard_Complex_Numbers.Complex_Number;
                want,maxprc,maxitr : in natural32;
                loss : out integer32; err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Performs a sequence of Newton steps to approximate a root
  --   of a polynomial homotopy up to a wanted number of accuracte
  --   decimal places taking the loss of decimal places into account
  --   when performing one Newton step on the homotopy at z.

  -- REQUIRED :
  --   The Varbprec_Homotopy is initialized properly.

  -- ON ENTRY :
  --   file     for intermediate output, if omitted,
  --            then the procedure will be silent;
  --   z        string representation of an initial approximation;
  --   t        corresponding value of the homotopy continuation parameter;
  --   want     wanted number of accurate decimal places;
  --   maxprc   maximum number of decimal places that can be used
  --            to estimate the loss of accuracy;
  --   maxitr   maximum number of Newton steps.

  -- ON RETURN :
  --   z        updated approximation for a solution of the homotopy at t;
  --   loss     estimated loss of decimal places as a negative number;
  --   err      magnitude of the correction added to z;
  --   rco      estimate of the inverse condition number of Jacobian at z;
  --   res      magnitude of the homotopy evaluated at z and at t.

end Varbprec_Corrector_Steps;
