with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Double_Double_Vectors;
with Double_Double_VecVecs;
with DoblDobl_Complex_Vectors;           use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;          use DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Solutions;         use DoblDobl_Complex_Solutions;

package DoblDobl_IncFix_Continuation is

-- DESCRIPTION :
--   The procedures below implement an increment-and-fix continuation
--   method with double double numbers.  The generic parameters are 
--   a norm function, an evaluator and a differentiator of the homotopy.
--   There are two basic versions: a silent and a reporting one.
--   The silent continuation simply performs its calculations without output
--   of intermediate results.  The reporting continuation routine allows to
--   put various kinds of intermediate results on a file.
--   It is assumed that the continuation parameters are already determined
--   before calling these routines (see Continuation_Parameters).

-- THE SILENT VERSIONS :

  generic

    with function Norm ( x : Vector ) return double_double;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Silent_Continue
               ( sols : in out Solution_List;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) );

  -- DESCRIPTION :
  --   Basic version of path tracking.

  -- ON ENTRY :
  --   sols      start solutions;
  --   nbq       number of equations to call the Gauss-Newton correctors;
  --   target    value of the continuation parameter at the end.

  -- ON RETURN :
  --   sols      the computed solutions.

  generic

    with function Norm ( x : Vector ) return double_double;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
    with function Stop_Test ( s : Solution ) return boolean;

  procedure Silent_Continue_with_Stop
               ( sols : in out Solution_List;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) );

  -- DESCRIPTION :
  --   Path tracking stops as soon as one of the end solutions
  --   meets the criterion provided by the Stop_Test.

  -- ON ENTRY :
  --   sols      start solutions;
  --   nbq       number of equations to call the Gauss-Newton correctors;
  --   target    value of the continuation parameter at the end.

  -- ON RETURN :
  --   sols      the computed solutions.

-- THE REPORTING VERSIONS :

  generic

    with function Norm ( x : Vector ) return double_double;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Reporting_Continue
               ( file : in file_type; sols : in out Solution_List;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) );

  -- DESCRIPTION :
  --   This routine implements the continuation strategy.

  -- ON ENTRY :
  --   file      to write intermediate results on (if Reporting_);
  --   sols      the start solutions;
  --   nbq       number of equations to call the Gauss-Newton correctors;
  --   target    value for the continuation parameter at the end.
 
  -- ON RETURN :
  --   sols      the computed solutions.

  generic

    with function Norm ( x : Vector ) return double_double;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
    with function Stop_Test ( s : Solution ) return boolean;

  procedure Reporting_Continue_with_Stop
               ( file : in file_type; sols : in out Solution_List;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) );

  -- DESCRIPTION :
  --   Path tracking stops as soon as one of the end solutions
  --   meets the criterion provided by the Stop_Test.

  -- ON ENTRY :
  --   file      to write intermediate results on (if Reporting_);
  --   sols      the start solutions;
  --   nbq       number of equations to call the Gauss-Newton correctors;
  --   target    value for the continuation parameter at the end.
 
  -- ON RETURN :
  --   sols      the computed solutions.

-- WITH THE ESTIMATION OF THE PATH DIRECTIONS :

  generic

    with function Norm ( x : Vector ) return double_double;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Silent_Toric_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Double_Double_VecVecs.VecVec;
                 errv : in out Double_Double_Vectors.Vector;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) );

  generic

    with function Norm ( x : Vector ) return double_double;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;

  procedure Reporting_Toric_Continue
               ( file : in file_type;
                 sols : in out Solution_List; proj : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Double_Double_VecVecs.VecVec;
                 errv : in out Double_Double_Vectors.Vector;
                 nbq : in integer32 := 0;
                 target : in Complex_Number := Create(integer(1)) );

  -- DESCRIPTION :
  --   This routine implements the continuation strategy with the estimation
  --   of the directions of the solution paths at the end.

  -- ON ENTRY :
  --   file      to write intermediate results on (if Reporting_);
  --   sols      the start solutions;
  --   proj      for projective-perpendicular path following;
  --   w         values for the winding numbers initialized at one,
  --             w'range must be 1..Length_Of(sols);
  --   v         v must be initialized with zero vectors
  --             and v'range is 1..Length_Of(sols);
  --   errv      errors on the computed directions;
  --   nbq       number of equations to call the Gauss-Newton correctors;
  --   target    value for the continuation parameter at the end.

  -- ON RETURN :
  --   sols      the computed solutions;
  --   w         estimated values for the winding numbers of the paths;
  --   v         directions of the solution paths;
  --   errv      errors on the computed directions.

end DoblDobl_IncFix_Continuation;
