with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;

package Hypersurfaces_and_Filters is

-- DESCRIPTION :
--   This package provides utilities to compute representations
--   of witness sets for hypersurfaces and to filter points by 
--   means of equations of hypersurfaces.

  generic

    with function f ( x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Numbers.Complex_Number;

  procedure RG_Hypersurface_Witness_Set
              ( file : in file_type; n,d : in natural32;
                b,v : in Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Solutions.Solution_List;
                res : out double_float );
  generic

    with function f ( x : DoblDobl_Complex_Vectors.Vector )
                    return DoblDobl_Complex_Numbers.Complex_Number;

  procedure RG_DoblDobl_Hypersurface_Witness_Set
              ( file : in file_type; n,d : in natural32;
                b,v : in DoblDobl_Complex_Vectors.Vector;
                s : out DoblDobl_Complex_Solutions.Solution_List;
                res : out double_double );
  generic

    with function f ( x : QuadDobl_Complex_Vectors.Vector )
                    return QuadDobl_Complex_Numbers.Complex_Number;

  procedure RG_QuadDobl_Hypersurface_Witness_Set
              ( file : in file_type; n,d : in natural32;
                b,v : in QuadDobl_Complex_Vectors.Vector;
                s : out QuadDobl_Complex_Solutions.Solution_List;
                res : out quad_double );

  -- DESCRIPTION :
  --   Reporting version of hypersurface witness set computation,
  --   where the function f gives the values of the polynomial,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of variables;
  --   d        degree of the polynomial;
  --   b        offset vector of a random line;
  --   v        direction of a random line.

  -- ON RETURN :
  --   s        witness set for p on line b + t*v;
  --   res      max norm of the residual vector.

  procedure RP_Hypersurface_Witness_Set
              ( file : in file_type; n,d : in natural32;
                p : in Standard_Complex_Poly_Functions.Eval_Poly;
                b,v : in Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Solutions.Solution_List;
                res : out double_float );
  procedure RP_Hypersurface_Witness_Set
              ( file : in file_type; n,d : in natural32;
                p : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                b,v : in DoblDobl_Complex_Vectors.Vector;
                s : out DoblDobl_Complex_Solutions.Solution_List;
                res : out double_double );
  procedure RP_Hypersurface_Witness_Set
              ( file : in file_type; n,d : in natural32;
                p : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                b,v : in QuadDobl_Complex_Vectors.Vector;
                s : out QuadDobl_Complex_Solutions.Solution_List;
                res : out quad_double );

  -- DESCRIPTION :
  --   Reporting version of hypersurface witness set computation,
  --   in standard double, double double, or quad double precision.

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

    with function f ( x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Numbers.Complex_Number;

  procedure SG_Hypersurface_Witness_Set
              ( n,d : in natural32;
                b,v : in Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Solutions.Solution_List;
                res : out double_float );

  procedure SP_Hypersurface_Witness_Set
              ( n,d : in natural32;
                p : in Standard_Complex_Poly_Functions.Eval_Poly;
                b,v : in Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Solutions.Solution_List;
                res : out double_float );

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

    with function f ( k : integer32;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Numbers.Complex_Number;

  procedure RG_Filter
              ( file : in file_type; ne : in natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                b,v : in Standard_Complex_Vectors.Vector;
                tol : in double_float );

  procedure RP_Filter
              ( file : in file_type;
                sols : in out Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                b,v : in Standard_Complex_Vectors.Vector;
                tol : in double_float );

  -- DESCRIPTION :
  --   Reporting version of filter witness set with polynomial equations.

  -- ON ENTRY :
  --   file     for intermediate output;
  --   ne       number of equations in p, i.e.: p'range = 1..ne;
  --   sols     list of candidate witness points;
  --   p        polynomial functions to evaluate sols in;
  --   b        offset vector of a random line;
  --   v        direction vector of a random line;
  --   tol      tolerance to decide whether number is zero.

  -- ON RETURN :
  --   sols     only those solutions whose evaluation in all
  --            polynomial in p was larger than tol.

  generic

    with function f ( k : integer32;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Numbers.Complex_Number;

  procedure SG_Filter
              ( ne : in natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                b,v : in Standard_Complex_Vectors.Vector;
                tol : in double_float );

  procedure SP_Filter
              ( sols : in out Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                b,v : in Standard_Complex_Vectors.Vector;
                tol : in double_float );

  -- DESCRIPTION :
  --   Silent version of filter witness set with polynomial equations.

  -- ON ENTRY :
  --   ne       number of equations in p, i.e.: p'range = 1..ne;
  --   sols     list of candidate witness points;
  --   p        polynomial functions to evaluate sols in;
  --   b        offset vector of a random line;
  --   v        direction vector of a random line;
  --   tol      tolerance to decide whether number is zero.

  -- ON RETURN :
  --   sols     only those solutions whose evaluation in all
  --            polynomials in p was larger than tol.

  generic

    with function f ( x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Numbers.Complex_Number;

  procedure RG_Split_Filter
              ( file : in file_type; tol : in double_float;
                p_sols : in out Standard_Complex_Solutions.Solution_List;
                q_sols : out Standard_Complex_Solutions.Solution_List;
                plane : in Standard_Complex_Matrices.Matrix );

  procedure RP_Split_Filter
              ( file : in file_type;
                p : in Standard_Complex_Poly_Functions.Eval_Poly;
                tol : in double_float;
                p_sols : in out Standard_Complex_Solutions.Solution_List;
                q_sols : out Standard_Complex_Solutions.Solution_List;
                plane : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Splits the given solution list into two lists: those which satisfy
  --   p and those which do not, with respect to the given tolerance.
  --   This reporting version writes extra output to file.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   p        polynomial function, to evaluate a new hypersurface;
  --   tol      tolerance to decide whether number is zero;
  --   p_sols   a witness set for previous equations;
  --   plane    plane which cuts out the solutions in p_sols.

  -- ON RETURN :
  --   p_sols   witness set which satisfies p;
  --   q_sols   solutions which do not satisfy p become start solutions.

  generic

    with function f ( x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Numbers.Complex_Number;

  procedure SG_Split_Filter
              ( tol : in double_float;
                p_sols : in out Standard_Complex_Solutions.Solution_List;
                q_sols : out Standard_Complex_Solutions.Solution_List;
                plane : in Standard_Complex_Matrices.Matrix );

  procedure SP_Split_Filter
              ( p : in Standard_Complex_Poly_Functions.Eval_Poly;
                tol : in double_float;
                p_sols : in out Standard_Complex_Solutions.Solution_List;
                q_sols : out Standard_Complex_Solutions.Solution_List;
                plane : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Splits the given solution list into two lists: those which satisfy
  --   p and those which do not, with respect to the given tolerance.
  --   This is a silent version without any output.

  -- ON ENTRY :
  --   p        polynomial function, to evaluate a new hypersurface;
  --   tol      tolerance to decide whether number is zero;
  --   p_sols   a witness set for previous equations;
  --   plane    plane which cuts out p_sols.

  -- ON RETURN :
  --   p_sols   witness set which satisfies p;
  --   q_sols   solutions which do not satisfy p become start solutions.

  generic

    with function Q ( x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector;

  procedure QSG_Filter
              ( s : in out Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Matrices.Matrix;
                tol : in double_float );

  -- DESCRIPTION :
  --   Evaluates every solution of s in Q and removes those solutions
  --   whose evaluation vector in Q has a component less than tol.

  -- ON ENTRY :
  --   s        solutions in intrinsic coordinates;
  --   p        affine basis for a plane;
  --   tol      tolerance to decide whether a number is zero.

  generic

    with function Q ( x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector;

  procedure QRG_Filter
              ( file : in file_type;
                s : in out Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Matrices.Matrix;
                tol : in double_float );

  -- DESCRIPTION :
  --   This is the reporting version of the Q_Filter, same parameters
  --   as above, only with a file for extra output.

end Hypersurfaces_and_Filters;
