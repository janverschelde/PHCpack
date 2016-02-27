with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
 
package Standard_IncFix_Continuation is

-- DESCRIPTION :
--   This package provides a general implementation of an increment-and-fix
--   continuation method.  The generic parameters are a norm, an evaluator and
--   a differentiator of the homotopy.
--   There are two basic versions: a silent and a reporting one.
--   The silent continuation simply performs its calculations without output
--   of intermediate results.  The reporting continuation routine allows to
--   put various kinds of intermediate results on a file.
--   It is assumed that the continuation parameters are already determined
--   before calling these routines (see Continuation_Parameters).
--   For both the silent and the reporting version, the facility is added
--   to estimate the directions of the solution paths, useful in a polyhedral
--   end game.

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Silent_Small_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Silent_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
    with function Stop_Test ( s : Solution ) return boolean;

  procedure Silent_Continue_with_Stop
               ( sols : in out Solution_List; proj : in boolean;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Reporting_Small_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean; nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Reporting_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean; nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
    with function Stop_Test ( s : Solution ) return boolean;

  procedure Reporting_Continue_with_Stop
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean; nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) );

  -- DESCRIPTION :
  --   This routine implements the continuation strategy.
  --   The sequential tracking stops when the Stop_Test returns true.
  --   For large systems, use the "Large" option.

  -- ON ENTRY :
  --   file      to write intermediate results on (if Reporting_);
  --   sols      the start solutions;
  --   proj      for projective-perpendicular path following;
  --   nbq       number of equations to call the Gauss-Newton correctors;
  --   target    value for the continuation parameter at the end.
 
  -- ON RETURN :
  --   sols      the computed solutions.

-- WITH THE ESTIMATION OF THE PATH DIRECTIONS :

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Silent_Toric_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Standard_Floating_VecVecs.VecVec;
                 errv : in out Standard_Floating_Vectors.Vector;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Reporting_Toric_Continue
               ( file : in file_type;
                 sols : in out Solution_List; proj : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Standard_Floating_VecVecs.VecVec;
                 errv : in out Standard_Floating_Vectors.Vector;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(1.0) );

  -- DESCRIPTION :
  --   This routine implements the continuation strategy with the estimation
  --   of the directions of the solution paths at the end.

  -- ON ENTRY :
  --   file      to write intermediate results on (if Reporting_);
  --   sols      the start solutions;
  --   proj      for projective-perpendicular path following;
  --   w         w values for the winding number initialized at 1,
  --             w'range must be 1..Length_Of(sols);
  --   v         v must be initialized with zero vectors
  --             and v'range is 1..Length_Of(sols);
  --   errv      errors on the computed directions;
  --   nbq       number of equations to call the Gauss-Newton correctors;
  --   target    value for the continuation parameter at the end.

  -- ON RETURN :
  --   sols      the computed solutions;
  --   w         estimated values for the winding numbers of the paths.
  --   v         directions of the solution paths;
  --   errv      errors on the computed directions.

end Standard_IncFix_Continuation;
