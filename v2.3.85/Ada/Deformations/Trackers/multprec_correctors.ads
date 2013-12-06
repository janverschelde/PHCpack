with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;          use Multprec_Complex_Matrices;
with Multprec_Continuation_Data;         use Multprec_Continuation_Data;

package Multprec_Correctors is

-- DESCRIPTION :
--   This package contains implementations for the corrector in an
--   increment-and-fix continuation with multi-precision arithmetic.
--
--   The following options can be made :
--    (Affine,Projective)
--       An affine corrector works in affine space, while a projective
--       corrector is a projective-perpendicular corrector: it works in
--       projective space and corrects in a perpendicular way.
--    (Single,Multiple)
--       A single corrector only deals with one path at a time.
--       A multiple corrector corrects more than one path when it is called.
--    (Loose,Severe)
--       A loose corrector will stop when either one of the following
--       conditions is satisfied:
--         1. One of the desired accuracies has been met.
--         2. The maximum number of iterations is reached.
--         3. The Jacobian matrix is singular.
--       In addition to these stopping criteria, a severe corrector checks
--       the convergence during the iterations and stops when it notices
--       divergence is noticed.  A loose correctors allows divergence.
--    (Normal,Conditioned)
--       A normal corrector does not compute an estimate for the inverse of
--       the condition number of the Jacobian matrix.  This additional work
--       is done by a conditioned corrector.
--    (Silent,Reporting)
--       A silent corrector does not produce any output on file.
--       A reporting corrector allows to put intermediate results on file.
--
--   Based on these options, the following 32 different correctors 
--   are provided :
--
--     Affine_Single_Loose_Normal_Silent_Corrector
--     Affine_Single_Loose_Normal_Reporting_Corrector
--     Affine_Single_Loose_Conditioned_Silent_Corrector
--     Affine_Single_Loose_Conditioned_Reporting_Corrector
--     Affine_Single_Severe_Normal_Silent_Corrector
--     Affine_Single_Severe_Normal_Reporting_Corrector
--     Affine_Single_Severe_Conditioned_Silent_Corrector
--     Affine_Single_Severe_Conditioned_Reporting_Corrector
--     Affine_Multiple_Loose_Normal_Silent_Corrector
--     Affine_Multiple_Loose_Normal_Reporting_Corrector
--     Affine_Multiple_Loose_Conditioned_Silent_Corrector
--     Affine_Multiple_Loose_Conditioned_Reporting_Corrector
--     Affine_Multiple_Severe_Normal_Silent_Corrector
--     Affine_Multiple_Severe_Normal_Reporting_Corrector
--     Affine_Multiple_Severe_Conditioned_Silent_Corrector
--     Affine_Multiple_Severe_Conditioned_Reporting_Corrector
--     Projective_Single_Loose_Normal_Silent_Corrector
--     Projective_Single_Loose_Normal_Reporting_Corrector
--     Projective_Single_Loose_Conditioned_Silent_Corrector
--     Projective_Single_Loose_Conditioned_Reporting_Corrector
--     Projective_Single_Severe_Normal_Silent_Corrector
--     Projective_Single_Severe_Normal_Reporting_Corrector
--     Projective_Single_Severe_Conditioned_Silent_Corrector
--     Projective_Single_Severe_Conditioned_Reporting_Corrector
--     Projective_Multiple_Loose_Normal_Silent_Corrector
--     Projective_Multiple_Loose_Normal_Reporting_Corrector
--     Projective_Multiple_Loose_Conditioned_Silent_Corrector
--     Projective_Multiple_Loose_Conditioned_Reporting_Corrector
--     Projective_Multiple_Severe_Normal_Silent_Corrector
--     Projective_Multiple_Severe_Normal_Reporting_Corrector
--     Projective_Multiple_Severe_Conditioned_Silent_Corrector
--     Projective_Multiple_Severe_Conditioned_Reporting_Corrector
--
--   All these procedures have the following generic parameters:
--     a norm function, polynomial vector and Jacobian matrix function.
--   Note that the projective correctors require a homogeneous polynomial
--   vector function.

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Affine_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info; c : in Corr_Pars );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  generic
    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Projective_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  -- DESCRIPTION :
  --   The predicted solutions of the system H(x,t)=0 are corrected.
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

end Multprec_Correctors;
