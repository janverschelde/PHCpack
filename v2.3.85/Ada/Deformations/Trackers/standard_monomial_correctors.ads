with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Continuation_Data;         use Standard_Continuation_Data;

package Standard_Monomial_Correctors is

-- DESCRIPTION :
--   This package offers implementation of a corrector method for use in
--   an increment-and-fix continuation method.
--   To evaluate the polynomials and their derivatives, a flattened sparse
--   monomial representation is expected to instantiate the routines.
--   The generic parameters of all correctors are
--     Norm : a vector norm,
--     V : evaluates the monomial vector of the homotopy;
--     F : given the monomial vector, evaluates the polynomials;
--     J : given the monomial vector, evaluates the derivatives.
--   Four different flavors of the correctors are currently implemented:
--   (normal,conditioned) : a normal corrector does not estimate the
--      inverse of the condition number, but a conditioned one does;
--   (silent,reporting) : a reporting corrector writes intermediate output
--       to file, whereas a silent one stays silent.
--   The other qualifiers mean the following :
--     affine : instead of projective space;
--     single : corrects only one solution;
--     severe : aborts as soon as divergence is observed.

  generic
     with function Norm ( x : Vector ) return double_float;
     with function V ( x : Vector; t : Complex_Number ) return Vector;
     with function F ( vxt : Vector ) return Vector;
     with function J ( vxt : Vector ) return Matrix;
  procedure Affine_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
     with function Norm ( x : Vector ) return double_float;
     with function V ( x : Vector; t : Complex_Number ) return Vector;
     with function F ( vxt : Vector ) return Vector;
     with function J ( vxt : Vector ) return Matrix;
  procedure Affine_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
     with function Norm ( x : Vector ) return double_float;
     with function V ( x : Vector; t : Complex_Number ) return Vector;
     with function F ( vxt : Vector ) return Vector;
     with function J ( vxt : Vector ) return Matrix;
  procedure Affine_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars );

  generic
     with function Norm ( x : Vector ) return double_float;
     with function V ( x : Vector; t : Complex_Number ) return Vector;
     with function F ( vxt : Vector ) return Vector;
     with function J ( vxt : Vector ) return Matrix;
  procedure Affine_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars );

  -- DESCRIPTION :
  --   Corrects the solution to the homotopy defined by F
  --   and Jacobian matrix at J.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   s        current solution;
  --   c        numerical settings for the parameters.

  -- ON RETURN :
  --   s        corrected solution.

end Standard_Monomial_Correctors;
