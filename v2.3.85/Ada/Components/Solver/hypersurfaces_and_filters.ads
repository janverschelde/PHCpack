with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;

package Hypersurfaces_and_Filters is

-- DESCRIPTION :
--   This package provides utilities to compute representations
--   of witness sets for hypersurfaces and to filter points by 
--   means of equations of hypersurfaces.

  generic

    with function f ( x : Vector ) return Complex_Number;

  procedure RG_Hypersurface_Witness_Set
              ( file : in file_type; n,d : in natural32;
                b,v : in Vector; s : out Solution_List;
                res : out double_float );

  procedure RP_Hypersurface_Witness_Set
              ( file : in file_type; n,d : in natural32;
                p : in Eval_Poly; b,v : in Vector;
                s : out Solution_List; res : out double_float );

  -- DESCRIPTION :
  --   Reporting version of hypersurface witness set computation.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of variables;
  --   d        degree of the polynomial;
  --   p        polynomial function to evaluate;
  --   b        offset vector of a random line;
  --   v        direction of a random line.

  -- ON RETURN :
  --   s        witness set for p on line b + t*v;
  --   res      max norm of the residual vector.

  generic

    with function f ( x : Vector ) return Complex_Number;

  procedure SG_Hypersurface_Witness_Set
              ( n,d : in natural32; b,v : in Vector;
                s : out Solution_List; res : out double_float );

  procedure SP_Hypersurface_Witness_Set
              ( n,d : in natural32;
                p : in Eval_Poly; b,v : in Vector;
                s : out Solution_List; res : out double_float );

  -- DESCRIPTION :
  --   Silent version of hypersurface witness set computation.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        degree of the polynomial;
  --   p        polynomial function to evaluate;
  --   b        offset vector of a random line;
  --   v        direction of a random line.

  -- ON RETURN :
  --   s        witness set for p on line b + t*v;
  --   res      max norm of the residual vector.

  generic

    with function f ( k : integer32; x : Vector ) return Complex_Number;

  procedure RG_Filter ( file : in file_type; ne : in natural32;
                        sols : in out Solution_List;
                        b,v : in Vector; tol : in double_float );

  procedure RP_Filter ( file : in file_type; sols : in out Solution_List;
                        p : in Eval_Poly_Sys; b,v : in Vector;
                        tol : in double_float );

  -- DESCRIPTION :
  --   Reporting version of filter witness set with polynomial equations.

  -- ON ENTRY :
  --   file       for intermediate output;
  --   ne         number of equations in p, i.e.: p'range = 1..ne;
  --   sols       list of candidate witness points;
  --   p          polynomial functions to evaluate sols in;
  --   b          offset vector of a random line;
  --   v          direction vector of a random line;
  --   tol        tolerance to decide whether number is zero.

  -- ON RETURN :
  --   sols       only those solutions whose evaluation in all
  --              polynomial in p was larger than tol.

  generic

    with function f ( k : integer32; x : Vector ) return Complex_Number;

  procedure SG_Filter ( ne : in natural32; sols : in out Solution_List;
                        b,v : in Vector; tol : in double_float );

  procedure SP_Filter ( sols : in out Solution_List;
                        p : in Eval_Poly_Sys; b,v : in Vector;
                        tol : in double_float );

  -- DESCRIPTION :
  --   Silent version of filter witness set with polynomial equations.

  -- ON ENTRY :
  --   ne         number of equations in p, i.e.: p'range = 1..ne;
  --   sols       list of candidate witness points;
  --   p          polynomial functions to evaluate sols in;
  --   b          offset vector of a random line;
  --   v          direction vector of a random line;
  --   tol        tolerance to decide whether number is zero.

  -- ON RETURN :
  --   sols       only those solutions whose evaluation in all
  --              polynomial in p was larger than tol.

  generic

    with function f ( x : Vector ) return Complex_Number;

  procedure RG_Split_Filter
               ( file : in file_type; tol : in double_float;
                 p_sols : in out Solution_List; q_sols : out Solution_List;
                 plane : in Matrix );

  procedure RP_Split_Filter
               ( file : in file_type; p : in Eval_Poly; tol : in double_float;
                 p_sols : in out Solution_List; q_sols : out Solution_List;
                 plane : in Matrix );

  -- DESCRIPTION :
  --   Splits the given solution list into two lists: those which satisfy
  --   p and those which do not, with respect to the given tolerance.
  --   This reporting version writes extra output to file.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   p         polynomial function, to evaluate a new hypersurface;
  --   tol       tolerance to decide whether number is zero;
  --   p_sols    a witness set for previous equations;
  --   plane     plane which cuts out the solutions in p_sols.

  -- ON RETURN :
  --   p_sols    witness set which satisfies p;
  --   q_sols    solutions which do not satisfy p become start solutions.

  generic

    with function f ( x : Vector ) return Complex_Number;

  procedure SG_Split_Filter
               ( tol : in double_float;
                 p_sols : in out Solution_List; q_sols : out Solution_List;
                 plane : in Matrix );

  procedure SP_Split_Filter
               ( p : in Eval_Poly; tol : in double_float;
                 p_sols : in out Solution_List; q_sols : out Solution_List;
                 plane : in Matrix );

  -- DESCRIPTION :
  --   Splits the given solution list into two lists: those which satisfy
  --   p and those which do not, with respect to the given tolerance.
  --   This is a silent version without any output.

  -- ON ENTRY :
  --   p         polynomial function, to evaluate a new hypersurface;
  --   tol       tolerance to decide whether number is zero;
  --   p_sols    a witness set for previous equations;
  --   plane     plane which cuts out p_sols.

  -- ON RETURN :
  --   p_sols    witness set which satisfies p;
  --   q_sols    solutions which do not satisfy p become start solutions.

  generic
    with function Q ( x : Vector ) return Vector;
  procedure QSG_Filter ( s : in out Solution_List;
                         p : in Matrix; tol : in double_float );

  -- DESCRIPTION :
  --   Evaluates every solution of s in Q and removes those solutions
  --   whose evaluation vector in Q has a component less than tol.

  -- ON ENTRY :
  --   s         solutions in intrinsic coordinates;
  --   p         affine basis for a plane;
  --   tol       tolerance to decide whether a number is zero.

  generic

    with function Q ( x : Vector ) return Vector;

  procedure QRG_Filter ( file : in file_type; s : in out Solution_List;
                         p : in Matrix; tol : in double_float );

  -- DESCRIPTION :
  --   This is the reporting version of the Q_Filter, same parameters
  --   as above, only with a file for extra output.

end Hypersurfaces_and_Filters;
