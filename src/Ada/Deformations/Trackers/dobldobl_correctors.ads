with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;           use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;          use DoblDobl_Complex_Matrices;
with Continuation_Parameters;            use Continuation_Parameters;
with DoblDobl_Continuation_Data;         use DoblDobl_Continuation_Data;

package DoblDobl_Correctors is

-- DESCRIPTION :
--   This package offers implementation of a corrector method for use in
--   an increment-and-fix continuation method.
--   The generic parameters of all correctors are
--     Norm : a vector norm,
--     Fun : given (x,t) returns the value of F at (x,t),
--     Jef : given (x,t) returns the value of the Jacobian matrix at (x,t),
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
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
     with function Norm ( x : Vector ) return double_double;
     with function Fun ( x : Vector; t : Complex_Number ) return Vector;
     with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
     with function Norm ( x : Vector ) return double_double;
     with function Fun ( x : Vector; t : Complex_Number ) return Vector;
     with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
     with function Norm ( x : Vector ) return double_double;
     with function Fun ( x : Vector; t : Complex_Number ) return Vector;
     with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars );

  generic
     with function Norm ( x : Vector ) return double_double;
     with function Fun ( x : Vector; t : Complex_Number ) return Vector;
     with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars );

  -- DESCRIPTION :
  --   Corrects the solution s to the homotopy defined by Fun
  --   and Jacobian matrix defined by Jef using parameters in c.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   s        current solution;
  --   c        numerical settings for the parameters.

  -- ON RETURN :
  --   s        corrected solution.

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  -- DESCRIPTION :
  --   The predicted solutions of the system Fun(x,t)=0 are corrected.
  --   With a multiple corrector, the correction starts at s(pivot).

  -- ON ENTRY :
  --   file       to write intermediate results on;
  --   s          are the predicted values for the solutions;
  --   pivot      is the index in the array where the correction
  --              process must start;
  --   dist_sols  two solutions x1 and x2 are different 
  --              if for some k in 1..n : | x1(k) - x2(k) | > dist_sols;
  --   c          the corrector parameters.
 
  -- ON RETURN :
  --   s          the computed solutions;
  --   pivot      if fail then pivot is the index in the array where
  --              a difficulty occured; otherwise it is the same pivot as
  --              on entry;
  --   fail       is false if all solutions are computed with the desired
  --              precision eps, within the maximum number of allowed
  --              iterations, and if all solutions are different.

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return double_double;
    with function Fun ( x : Vector; t : Complex_Number ) return Vector;
    with function Jef ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  -- DESCRIPTION :
  --   The predicted solutions of the system Fun(x,t)=0 are corrected.
  --   With a multiple corrector, the correction starts at s(pivot).

  -- ON ENTRY :
  --   file       to write intermediate results on;
  --   s          are the predicted values for the solutions;
  --   pivot      is the index in the array where the correction
  --              process must start;
  --   dist_sols  two solutions x1 and x2 are different 
  --              if for some k in 1..n : | x1(k) - x2(k) | > dist_sols;
  --   c          the corrector parameters.
 
  -- ON RETURN :
  --   s          the computed solutions;
  --   pivot      if fail then pivot is the index in the array where
  --              a difficulty occured; otherwise it is the same pivot as
  --              on entry;
  --   fail       is false if all solutions are computed with the desired
  --              precision eps, within the maximum number of allowed
  --              iterations, and if all solutions are different.

end DoblDobl_Correctors;
