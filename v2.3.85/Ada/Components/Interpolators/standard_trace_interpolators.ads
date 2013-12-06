with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Floating_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Standard_Stacked_Sample_Grids;      use Standard_Stacked_Sample_Grids;
with Standard_Nvariate_Interpolators;    use Standard_Nvariate_Interpolators;

package Standard_Trace_Interpolators is

-- DESCRIPTION :
--   This package offers a data abstraction and operations to create
--   and evaluate the trace form a polynomial in several variables,
--   interpolating through a grid of points.  The polynomials have
--   complex coefficients of standard machine precision.

  type Trace_Interpolator1 is private;   -- only for curves
  type Trace_Interpolator is private;    -- for general surfaces

-- CREATORS :

  function Create ( grid : Array_of_Standard_Sample_Lists )
                  return Trace_Interpolator1;
  function Create ( file : file_type; grid : Array_of_Standard_Sample_Lists )
                  return Trace_Interpolator1;

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through
  --   the samples in the grid.  With "file" as parameter, additional
  --   tests are performed and diagnostics are written to the file.

  -- ASSUMED : we are interpolating a curve.

  -- REQUIRED :
  --   1) grid'range = 0..d, where d is the degree of the curve;
  --   2) all samples in grid(i) lie in hyperplanes parallel to the
  --      hyperplanes used for the samples in grid(j);
  --   3) Length_Of(grid(i)) = d, for all i in grid'range.

  function Create ( grid : Array_of_Standard_Sample_Lists; i : integer32 )
                  return Vector;

  -- DESCRIPTION :
  --   Returns the coefficient vector for the i-th trace, for which only
  --   the samples lists in the range 0..i in the grid are needed.
  --   However: every sample list must contain d samples, where d is the
  --   degree of the curve.  The range of the vector on return is 0..i.
  --   The assumption is here also that we are interpolate a curve.

  function Create_on_Triangle
                  ( file : file_type; grid : Array_of_Standard_Sample_Lists )
                  return Trace_Interpolator1;
  function Create_on_Triangle
                  ( grid : Array_of_Standard_Sample_Lists )
                  return Trace_Interpolator1;

  -- DESCRIPTION :
  --   Creates the trace for of the interpolating polynomial through
  --   a triangular grid.  With "file" as parameter, additional tests
  --   are performed and diagnostics are written on file.

  -- ASSUMED : we are interpolating a curve.

  -- REQUIRED :
  --   1) grid'range = 0..d, where d is the degree of the curve
  --   2) all samples in grid(i) lie in hyperplanes parallel to the
  --      slicing planes in grid(j), j in 0..d;
  --   3) Length_Of(grid(i)) = d, for i = 0,1,
  --      and Length_Of(grid(i)) = d - i + 1, for i > 2.

  function Create ( grid : Stacked_Sample_Grid; d : integer32 ) 
                  return Trace_Interpolator;
  function Create ( file : file_type;
                    grid : Stacked_Sample_Grid; d : integer32 ) 
                  return Trace_Interpolator;

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through the
  --   samples in the grid, for any degree d and dimension.  With file,
  --   diagnostics and results of additional tests are written to file.

  function Create ( grid : Stacked_Sample_Grid; d,i : integer32 )
                  return Newton_Form;
  function Create ( file : file_type;
                    grid : Stacked_Sample_Grid; d,i : integer32 )
                  return Newton_Form;

  -- DESCRIPTION :
  --   Returns the i-th trace of the interpolator of degree d.
  --   With file as argument, extra tests are performed with output to file.

  function Expand ( t : Trace_Interpolator1 ) return Poly;
  function Expand ( t : Trace_Interpolator ) return Poly;

  -- DESCRIPTION :
  --   Expands the trace form into an ordinary polynomial in n variables.

-- SELECTORS :

  function Degree ( t : Trace_Interpolator1 ) return integer32;

  -- DESCRIPTION :
  --   Returns -1 if the trace interpolator is empty, 
  --   otherwise its degree is returned.

  function Dimension ( t : Trace_Interpolator1 ) return integer32;

  -- DESCRIPION :
  --   Returns -1 if the trace form of the interpolator is empty,
  --   otherwise the number of variables after the embedding is returned.

  function Trace ( t : Trace_Interpolator1; i : integer32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the coefficient vector of the trace of degree i.
  --   The vector of return has range 0..i and the j-th entry is
  --   the coordinate with the j-th power of the variable.

  -- REQUIRED : Degree(t) > 0 and 0 < i <= Degree(t).

  function Trace ( t : Trace_Interpolator; i : integer32 ) return Newton_Form;

  -- DESCRIPTION :
  --   Returns the Newton form for the trace of degree i.
  --   To evaluate the result, one must apply the projection first.

  -- REQUIRED : Degree(t) > 0 and 0 < i <= Degree(t).

-- EVALUATORS :

  function Eval ( t : Trace_Interpolator1; x : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the trace interpolator at the given vector.

  procedure Eval_Trace ( t : in Vector; d : in integer32;
                         sps : in Standard_Sample_List;
                         val,eva : out Complex_Number );

  -- DESCRIPTION :
  --   Evaluates a trace of degree d in the samples in the list.

  -- REQUIRED :
  --   The samples in sps all lie on the same slice.  If the trace
  --   is correct, then val - eva is of machine precision order.

  -- ON ENTRY :
  --   t        vector of range 0..d with coefficients of the trace;
  --   d        degree of the trace;
  --   sps      list of samples.

  -- ON RETURN :
  --   val      value at the trace;
  --   eva      evaluation at the sample list.

  procedure Eval_Trace ( t : in Newton_Form; d : in integer32;
                         sps : in Standard_Sample_List;
                         val,eva : out Complex_Number );

  -- DESCRIPTION :
  --   Higher dimensional analogue of the routine above.

  function Eval ( t : Trace_Interpolator; x : Vector ) return Complex_Number; 

  -- DESCRIPTION :
  --   Returns the value of the interpolating polynomial at the vector x.

  function Eval ( file : file_type;
                  t : Trace_Interpolator; x : Vector ) return Complex_Number; 

  -- DESCRIPTION :
  --   Writes intermediate output on file during evaluation.

-- DIAGNOSTICS :

  function Errors ( t : Trace_Interpolator1;
                    grid : Array_of_Standard_Sample_Lists )
                  return Standard_Floating_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the matrix of residuals after evaluation the interpolator
  --   at each sample point in the grid.  The range of the matrix on return
  --   is grid'range times 1..Length_Of(grid(i)).

  procedure Write_Errors ( file : in file_type; t : in Trace_Interpolator1;
                           grid : in Array_of_Standard_Sample_Lists );

  -- DESCRIPTION :
  --   Writes the results of the evaluation of the trace interpolator
  --   at the grid points to the file.

  function Maximal_Error ( residuals : Standard_Floating_Matrices.Matrix )
                         return double_float;
  function Maximal_Error ( t : Trace_Interpolator1;
                           grid : Array_of_Standard_Sample_Lists )
                         return double_float;

  -- DESCRIPTION :
  --   Returns the maximal error of the matrix of residuals after
  --   evaluation at all samples in the grid.

  function Maximal_Error ( t : Vector; d : integer32;
                           grid : Array_of_Standard_Sample_Lists ) 
                         return double_float;

  -- DESCRIPTION :
  --   Returns the maximal difference between the value at the trace and
  --   the actual computed power sum of degree d at the slices in the grid.

  procedure Write_Errors
              ( file : in file_type; t : in Trace_Interpolator;
                grid : in Stacked_Sample_Grid; maxerr : out double_float );

  -- DESCRIPTION :
  --   Writes all evaluations of t in the grid on file and returns in
  --   maxerr the maximal residual.

  function Maximal_Error
             ( t : Trace_Interpolator; grid : Stacked_Sample_Grid )
             return double_float;

  -- DESCRIPTION :
  --   Returns the maximal residual of t evaluated at the grid.

  function Maximal_Error
             ( file : file_type;
               t : Newton_Form; d,i : integer32; grid : Stacked_Sample_Grid )
             return double_float;

  -- DESCRIPTION :
  --   Returns the maximal difference between the values predicted by the
  --   trace t of degree i and the actual computed power sums of degree i
  --   at the slices in the full grid for degree d interpolation.
  --   The results of the comparisons are written to the file.

-- DESTRUCTORS :

  procedure Clear ( t : in out Trace_Interpolator1 );
  procedure Clear ( t : in out Trace_Interpolator );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the trace interpolator.

private

  type Trace_Interpolator1_Rep;
  type Trace_Interpolator1 is access Trace_Interpolator1_Rep;

  type Trace_Interpolator_Rep;
  type Trace_Interpolator is access Trace_Interpolator_Rep;

end Standard_Trace_Interpolators;
