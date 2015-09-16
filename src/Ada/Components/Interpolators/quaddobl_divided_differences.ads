with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;       use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;      use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Sample_Lists;              use QuadDobl_Sample_Lists;
with QuadDobl_Stacked_Sample_Grids;      use QuadDobl_Stacked_Sample_Grids;

package QuadDobl_Divided_Differences is

-- DESCRIPTION :
--   The Newton form of an interpolating polynomial in two variables
--   requires the interpolation points to lie on a grid in the plane.
--   Given a grid of samples on a planar curve, the operations of
--   this package create, evaluate and expand the Newton form of the
--   interpolating polynomial through those samples.
--   This package implements a bootstrapping technique to interpolate,
--   using quad double arithmetic.

  type Newton_Interpolator1 is private;    -- for curves (1D objects)
  type Newton_Form_Evaluator1 is private;
  type Newton_Taylor_Form1 is private;

  type Newton_Taylor_Form is private;      -- for general surfaces

-- CREATORS :

  function Create ( grid : Array_of_QuadDobl_Sample_Lists )
                  return Newton_Interpolator1;
  function Create ( grid : Array_of_QuadDobl_Sample_Lists;
                    c : Complex_Number ) return Newton_Interpolator1;

  -- DESCRIPTION :
  --   Creates the Newton form of the interpolating polynomial through
  --   the samples in the grid.  If c is not supplied as parameter,
  --   a random number for c is generated.  In the interpolation, c is
  --   used for normalization:  p(x(0),c) = 1 where x(0) is the first 
  --   coordinate of the rotated grid of samples (p is the interpolator).

  -- ASSUMED : we are interpolating a curve.

  -- REQUIRED :
  --   1) All samples in grid(i) lie in hyperplanes parallel to the
  --      hyperplanes used for the samples in grid(j);
  --   2) grid'range = 0..d, where d is the degree of the interpolator;
  --   3) Length_Of(grid(i)) = d, for all i in grid'range.

  function Create ( q : Newton_Interpolator1 ) return Newton_Form_Evaluator1;

  -- DESCRIPTION :
  --   Returns the encoding of a Horner scheme of the Newton form of
  --   the interpolating polynomial.

  function Create ( q : Newton_Interpolator1; y : Vector )
                  return Newton_Taylor_Form1;

  -- DESCRIPTION :
  --   Evaluates q at the (x,y) points and returns the Newton-Taylor
  --   form of the interpolating polynomial.

  function Create ( grid : Stacked_Sample_Grid; c : Complex_Number;
                    y : Vector ) return Newton_Taylor_Form;
  function Create ( file : file_type; grid : Stacked_Sample_Grid;
                    c : Complex_Number; y : Vector )
                  return Newton_Taylor_Form;

  -- DESCRIPTION :
  --   Returns the representation of the Newton interpolator for a surface
  --   of any degree and dimension.  Diagnostics are written to file.
  --   The constant "c" is the normalization constant used to create
  --   the intermediate Newton_Interpolator1, and y is a vector of range
  --   0..grid.d with y(0) = c of values for the second coordinate in
  --   the Newton_Taylor_Form1.  When the file is provided, intermediate
  --   results and output of test calculations are written on file.

  function Expand ( q : Newton_Form_Evaluator1 ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the symbolic expansion of the last row in the table 
  --   of divided differences.  The range of the system on return
  --   is 0..d, where d is the degree of the Newton Form Evaluator.
  --   Modulo roundoff, degree(res(i)) = i, for res on return.
  --   The polynomials on return are in two variables.

  function Expand ( q : Newton_Form_Evaluator1; dvd : Poly_Sys ) return Poly;
  function Expand ( q : Newton_Form_Evaluator1 ) return Poly;

  -- DESCRIPTION :
  --   Expands the table of divided differences symbolically into 
  --   an ordinary polynomial in n variables, n = Dimension(q)-1.
  --   Optionally, the parameter dvd = Expand(q).

-- SELECTORS :

  function Degree ( q : Newton_Interpolator1 ) return integer32;
  function Degree ( q : Newton_Form_Evaluator1 ) return integer32;

  -- DESCRIPTION :
  --   Returns the degree of the Newton form of the interpolator.

  function Dimension ( q : Newton_Interpolator1 ) return integer32;
  function Dimension ( q : Newton_Form_Evaluator1 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of variables after the embedding.

-- DIAGNOSTICS :

  function Errors ( q : Newton_Interpolator1;
                    grid : Array_of_QuadDobl_Sample_Lists ) return Matrix;

  -- DESCRIPTION :
  --   Returns the matrix of residuals after evaluation the interpolator
  --   at each sample point in the grid.  The range of the matrix on return
  --   is grid'range times 1..Length_Of(grid(i)).

  function Maximal_Error ( residuals : Matrix ) return double_float;
  function Maximal_Error
             ( q : Newton_Interpolator1;
               grid : Array_of_QuadDobl_Sample_Lists ) return double_float;

  -- DESCRIPTION :
  --   Returns the maximal error of the matrix of residuals after
  --   evaluation at all samples in the grid.

  function Maximal_Error
             ( q : Newton_Taylor_Form; grid : Stacked_Sample_Grid )
             return double_float;

  -- DESCRIPTION :
  --   Returns the maximal residual of q evaluated at the grid.

-- EVALUATORS :

  function Eval ( q : Newton_Interpolator1; x : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the table of divided differences for the given vector
  --   and returns the value of the interpolating polynomial.

  function Eval ( q : Newton_Form_Evaluator1; x : Vector )
                return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the interpolating polynomial at x,
  --   more efficiently than with the Eval on the Newton_Interpolator.

  function Eval ( q : Newton_Taylor_Form1; x : Vector ) return Complex_Number;
  function Eval ( q : Newton_Taylor_Form; x : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the interpolating polynomial at x.
  --   This is the most efficient way of evaluating the interpolator.

-- DESTRUCTORS :

  procedure Clear ( q : in out Newton_Interpolator1 );
  procedure Clear ( q : in out Newton_Form_Evaluator1 );
  procedure Clear ( q : in out Newton_Taylor_Form1 );
  procedure Clear ( q : in out Newton_Taylor_Form );

  -- DESCRIPTION :
  --   Releases all memory allocated to store the Newton form
  --   which ceases to exist after this Clear operation.

private

  type Newton_Interpolator1_Rep;
  type Newton_Interpolator1 is access Newton_Interpolator1_Rep;

  type Newton_Taylor_Form1_Rep;
  type Newton_Taylor_Form1 is access Newton_Taylor_Form1_Rep;

  type Newton_Form_Evaluator1_Rep;
  type Newton_Form_Evaluator1 is access Newton_Form_Evaluator1_Rep;

  type Newton_Taylor_Form_Rep;
  type Newton_Taylor_Form is access Newton_Taylor_Form_Rep;

end QuadDobl_Divided_Differences;
