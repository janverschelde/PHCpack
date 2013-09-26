with text_io;                            use text_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;          use Multprec_Complex_Matrices;
with Multprec_Complex_Solutions;         use Multprec_Complex_Solutions;

package Multprec_IncFix_Continuation is

-- DESCRIPTION :
--   The procedures below implement an increment-and-fix continuation
--   method with multi-precision arithmetic.  The generic parameters are 
--   a norm function, an evaluator and a differentiator of the homotopy.
--   There are two basic versions: a silent and a reporting one.
--   The silent continuation simply performs its calculations without output
--   of intermediate results.  The reporting continuation routine allows to
--   put various kinds of intermediate results on a file.
--   It is assumed that the continuation parameters are already determined
--   before calling these routines (see Continuation_Parameters).

  generic

    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Silent_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 target : in Complex_Number := Create(integer(1)) );

  generic

    with function Norm ( x : Vector ) return Floating_Number;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Reporting_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean;
                 target : in Complex_Number := Create(integer(1)) );

  -- DESCRIPTION :
  --   This routine implements the continuation strategy.

  -- ON ENTRY :
  --   file      to write intermediate results on (if Reporting_);
  --   sols      the start solutions;
  --   proj      for projective-perpendicular path following;
  --   target    value for the continuation parameter at the end.
 
  -- ON RETURN :
  --   sols      the computed solutions.

end Multprec_IncFix_Continuation;
