with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Vectors;
with Multprec_Complex_Solutions;

package Verification_of_Solutions is

-- DESCRIPTION :
--   A solution is verified up to d decimal places of accuracy
--   if the forward error and backward error are both less than 10^(-d)
--   for a bounded conditioned number estimated with large enough precision.
--   The operations in this package interface to running sequences of
--   Newton steps in variable precision based on condittion number estimates.
--   The distinction between regular polynomials and Laurent polynomials
--   (those with negative exponents) are not made when parsing the string
--   representations of the systems.  The caller routines must decide
--   which strings contain Laurent systems and which ones do not.

  procedure Set_Default_Parameters
              ( wanted,maxitr,maxprc : out natural32;
                verbose : out boolean );

  -- DESCRIPTION :
  --   Sets default parameters to run variable precision Newton steps.

  -- ON RETURN :
  --   wanted   number of accurate decimal places wanted;
  --   maxitr   maximum number of Newton steps in the sequence;
  --   maxprc   maximum number of decimal places in the precision
  --            to estimate the condition numbers;
  --   verbose  for intermediate output and diagnostics during the steps.

  procedure Write_Parameters
              ( file : in file_type;
                wanted,maxitr,maxprc : in natural32;
                verbose : in boolean );

  -- DESCRIPTION :
  --   Writes the parameters of variable precision Newton method to file.

  -- ON ENTRY :
  --   wanted   number of accurate decimal places wanted;
  --   maxitr   maximum number of Newton steps in the sequence;
  --   maxprc   maximum number of decimal places in the precision
  --            to estimate the condition numbers;
  --   verbose  for intermediate output and diagnostics during the steps.

  procedure Menu_to_set_Parameters
              ( wanted,maxitr,maxprc : out natural32;
                verbose : out boolean );

  -- DESCRIPTION :
  --   Shows the default settings of the parameters to the user
  --   and then allows the user to change their values.

  -- ON RETURN :
  --   wanted   number of accurate decimal places wanted;
  --   maxitr   maximum number of Newton steps in the sequence;
  --   maxprc   maximum number of decimal places in the precision
  --            to estimate the condition numbers;
  --   verbose  for intermediate output and diagnostics during the steps.

  function to_strings ( sols : Multprec_Complex_Solutions.Solution_List )
                      return Array_of_Strings;

  -- DESCRIPTION :
  --   Returns the solution vectors in sols to strings.
  --   After the Verify, we can do the reverse operation,
  --   see the to_solutions function below.

  procedure Verify_Solutions_of_Polynomial_System
              ( p : in Array_of_Strings; z : in out Array_of_Strings;
                wanted,maxitr,maxprec : in natural32;
                err,rco,res : out Standard_Floating_Vectors.Vector );
  procedure Verify_Solutions_of_Laurent_Polynomials
              ( p : in Array_of_Strings; z : in out Array_of_Strings;
                wanted,maxitr,maxprec : in natural32;
                err,rco,res : out Standard_Floating_Vectors.Vector );
  procedure Verify_Solutions_of_Polynomial_System
              ( file : in file_type;
                p : in Array_of_Strings; z : in out Array_of_Strings;
                wanted,maxitr,maxprec : in natural32;
                err,rco,res : out Standard_Floating_Vectors.Vector );
  procedure Verify_Solutions_of_Laurent_Polynomials
              ( file : in file_type;
                p : in Array_of_Strings; z : in out Array_of_Strings;
                wanted,maxitr,maxprec : in natural32;
                err,rco,res : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps on the system in p,
  --   starting at the solution vectors in z, with respect to the
  --   parameters of the variable precision Newton steps.

  -- ON ENTRY :
  --   file     if provided then intermediate output is written to file,
  --            otherwise the procedure remains silent;
  --   p        string representations of a (Laurent) polynomial system;
  --   z        string representations of the coordinates of solutions;
  --   wanted   number of accurate decimal places wanted;
  --   maxitr   maximum number of Newton steps in the sequence;
  --   maxprc   maximum number of decimal places in the precision
  --            to estimate the condition numbers;

  -- ON RETURN :
  --   z        updated coordinates for the solutions;
  --   err      forward error of the k-th solution is in err(k);
  --   rco      estimate for the inverse of the condition number of
  --            the Jacobian matrix of the k-th solution is in rco(k);
  --   res      backward error of the k-th solution is in res(k).

  function to_Solutions
              ( z : Array_of_Strings;
                err,rco,res : Standard_Floating_Vectors.Vector )
              return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the list of solutions using as solution vectors
  --   the strings in z, and the corresponding (err,rco,res) tuple
  --   as the diagnostics for every solution.

  -- REQUIRED : z'range = err'range = rco'range = res'range.

  procedure Verify_Solutions_of_Polynomial_System
              ( p : in Array_of_Strings;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                wanted,maxitr,maxprc : in natural32 );
  procedure Verify_Solutions_of_Laurent_Polynomials
              ( p : in Array_of_Strings;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                wanted,maxitr,maxprc : in natural32 );
  procedure Verify_Solutions_of_Polynomial_System
              ( file : in file_type; p : in Array_of_Strings;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                wanted,maxitr,maxprc : in natural32 );
  procedure Verify_Solutions_of_Laurent_Polynomials
              ( file : in file_type; p : in Array_of_Strings;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                wanted,maxitr,maxprc : in natural32 );

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps on the system in p,
  --   starting at the solutions in sols, with respect to the
  --   parameters of the variable precision Newton steps.

  -- ON ENTRY :
  --   file     if provided then intermediate output is written to file,
  --            otherwise the procedure remains silent;
  --   p        string representations of a (Laurent) polynomial system;
  --   sols     initial approximations to the solutions of p;
  --   wanted   number of accurate decimal places wanted;
  --   maxitr   maximum number of Newton steps in the sequence;
  --   maxprc   maximum number of decimal places in the precision
  --            to estimate the condition numbers;

  -- ON RETURN :
  --   sols     updated list of solutions.

end Verification_of_Solutions;
