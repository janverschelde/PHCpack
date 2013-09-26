with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Projection_Operators;               use Projection_Operators;
with Interpolation_Filters;              use Interpolation_Filters;
with Span_of_Component;                  use Span_of_Component;
with Standard_Divided_Differences;
with Multprec_Divided_Differences;
with Standard_Trace_Interpolators;
with Multprec_Trace_Interpolators;

package Irreducible_Component_Creators is

-- DESCRIPTION :
--   We create a representation for an irreducible component by
--   sampling and interpolating.  The incremental creators occur
--   in three groups :
--     1) plain linear projection operators;
--     2) plain linear projection with exploitation of the span;
--     3) central projection operators reduce degree of component.
--   Newton interpolation is used 
--     4) in bootstrapping fashion;
--     5) with the full trace form of the interpolator;
--     6) in the full trace form, exploiting the Newton identities.
--   Last but certainly not least is
--     7) the linear trace certification method.
--   Since linear polynomial are tolerant to roundoff, there is no
--   multi-precision version for (7).  All other interpolators are
--   implemented with standard and multi-precision floating arithmetic.

-- INTERPOLATION USING PLAIN LINEAR PROJECTION OPERATORS :

  procedure Standard_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Standard_Sample;
                   maxdeg : in natural32; tol : in double_float;
                   p : in Standard_Projector; f : out Standard_Filter );
  procedure Standard_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Standard_Sample; 
                   maxdeg : in natural32; tol : in double_float;
                   f : out Standard_Filter );
  procedure Multprec_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in out Standard_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   p : in Multprec_Projector; f : out Multprec_Filter );
  procedure Multprec_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Multprec_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   p : in Multprec_Projector; f : out Multprec_Filter );
  procedure Multprec_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in out Standard_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   f : out Multprec_Filter );
  procedure Multprec_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Multprec_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   f : out Multprec_Filter );

-- INTERPOLATION WITH LINEAR PROJECTIONS AND EXPLOITATION OF SPAN :

  procedure Standard_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Standard_Sample; 
                   maxdeg : in natural32; tol : in double_float;
                   p : in Standard_Projector;
                   f : out Standard_Filter; s : out Standard_Span );
  procedure Standard_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Standard_Sample; 
                   maxdeg : in natural32; tol : in double_float;
                   f : out Standard_Filter; s : out Standard_Span );
  procedure Multprec_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in out Standard_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   p : in Multprec_Projector;
                   f : out Multprec_Filter; s : out Multprec_Span );
  procedure Multprec_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in out Standard_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   f : out Multprec_Filter; s : out Multprec_Span );
  procedure Multprec_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Multprec_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   p : in Multprec_Projector;
                   f : out Multprec_Filter; s : out Multprec_Span );
  procedure Multprec_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Multprec_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   f : out Multprec_Filter; s : out Multprec_Span );

-- INTERPOLATION USING CENTRAL PROJECTIONS :

  procedure Standard_Central_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Standard_Sample; 
                   maxdeg : in natural32; tol : in double_float;
                   p : in Standard_Projector;
                   f : out Standard_Filter; s : out Standard_Span );
  procedure Standard_Central_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Standard_Sample; 
                   maxdeg : in natural32; tol : in double_float;
                   f : out Standard_Filter; s : out Standard_Span );
  procedure Multprec_Central_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in out Standard_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   p : in Multprec_Projector;
                   f : out Multprec_Filter; s : out Multprec_Span );
  procedure Multprec_Central_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in out Standard_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   f : out Multprec_Filter; s : out Multprec_Span );
  procedure Multprec_Central_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Multprec_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   p : in Multprec_Projector;
                   f : out Multprec_Filter; s : out Multprec_Span );
  procedure Multprec_Central_Interpolate
                 ( file : in file_type; full_output : in boolean;
                   spt : in Multprec_Sample;
                   maxdeg,size : in natural32; tol : in double_float;
                   f : out Multprec_Filter; s : out Multprec_Span );

  -- DESCRIPTION :
  --   Returns the filter and optionally also the span of the component.
  --   The "_Central_" interpolators use central projectors when the
  --   linear span of the component is strictly higher than the dimension
  --   of the component.

  -- REQUIRED : The sampling machine is initialized and properly tuned.

  -- ON ENTRY :
  --   file        for intermediate results and diagnostics;
  --   full_output indicates whether all intermediate sampling diagnostics
  --               are needed or whether only a summary will be written;
  --   spt         initial sample on the component, will be refined if
  --               multi-precision is invoked;
  --   maxdeg      upper bound on the degree, limits #samples;
  --   size        size of the numbers (for multi-precision arithmetic);
  --   tol         tolerance to for stopping the sampling;
  --   p           optional projector, if not provided, a random projector
  --               operator will be generated.

  -- ON RETURN :
  --   spt         when standard sample and Multprec, contains refinement;
  --   f           interpolation filter;
  --   s           span of the component.

-- NEWTON INTERPOLATION WITH DIVIDED DIFFERENCES :

  procedure Standard_Newton_Interpolate1
                 ( file : in file_type; sps : in Standard_Sample_List;
                   q : out Standard_Divided_Differences.Newton_Interpolator1;
                   eps,dist,gridres,testres : out double_float );
  procedure Standard_Newton_Interpolate
                 ( file : in file_type; sps : in Standard_Sample_List;
                   q : out Standard_Divided_Differences.Newton_Taylor_Form;
                   eps,dist,gridres,testres : out double_float );
  procedure Multprec_Newton_Interpolate1
                 ( file : in file_type; sps : in Standard_Sample_List;
                   size : in natural32;
                   q : out Multprec_Divided_Differences.Newton_Interpolator1;
                   eps,dist,gridres,testres : out double_float );
  procedure Multprec_Newton_Interpolate
                 ( file : in file_type; sps : in Standard_Sample_List;
                   size : in natural32;
                   q : out Multprec_Divided_Differences.Newton_Taylor_Form;
                   eps,dist,gridres,testres : out double_float );

  -- DESCRIPTION :
  --   Creates a Newton form of the interpolating filter,
  --   building a structured grid starting at the given list of samples.
  --   With the suffix "1", a curve is interpolated.

  -- ON ENTRY :
  --   file        for intermediate results and diagnostics;
  --   sps         list of starting sample points;
  --   size        size of the multi-precision numbers.

  -- ON RETURN :
  --   q           Newton form of the interpolator;
  --   eps         accuracy of the samples;
  --   dist        minimial distance between the samples;
  --   gridres     maximal residual of evaluation at all grid points;
  --   testres     maximal residual of evaluation at test points.

-- INTERPOLATION WITH TRACES :

  procedure Standard_Trace_Interpolate1
                 ( file : in file_type; sps : in Standard_Sample_List;
                   t : out Standard_Trace_Interpolators.Trace_Interpolator1;
                   eps,dist,gridres,testres : out double_float );
  procedure Standard_Power_Trace_Interpolate1
                 ( file : in file_type; sps : in Standard_Sample_List;
                   t : out Standard_Trace_Interpolators.Trace_Interpolator1;
                   eps,dist,gridres,testres : out double_float );
  procedure Standard_Trace_Interpolate
                 ( file : in file_type; sps : in Standard_Sample_List;
                   t : out Standard_Trace_Interpolators.Trace_Interpolator;
                   eps,dist,gridres,testres : out double_float );
  procedure Multprec_Trace_Interpolate1
                 ( file : in file_type; sps : in Standard_Sample_List;
                   size : in natural32;
                   t : out Multprec_Trace_Interpolators.Trace_Interpolator1;
                   eps,dist,gridres,testres : out double_float );
  procedure Multprec_Power_Trace_Interpolate1
                 ( file : in file_type; sps : in Standard_Sample_List;
                   size : in natural32;
                   t : out Multprec_Trace_Interpolators.Trace_Interpolator1;
                   eps,dist,gridres,testres : out double_float );
  procedure Multprec_Trace_Interpolate
                 ( file : in file_type; sps : in Standard_Sample_List;
                   size : in natural32;
                   t : out Multprec_Trace_Interpolators.Trace_Interpolator;
                   eps,dist,gridres,testres : out double_float );

  -- DESCRIPTION :
  --   Creates a trace form of the interpolating filter on a structured
  --   grid of samples.  With the suffix "1", a curve is interpolated.
  --   With "_Power_", the Newton identities are exploited.

  -- ON ENTRY :
  --   file        for intermediate results and diagnostics;
  --   sps         list of starting sample points;
  --   size        size of multi-precision numbers.

  -- ON RETURN :
  --   t           trace form of the interpolator;
  --   eps         accuracy of the samples;
  --   dist        minimial distance between the samples;
  --   gridres     maximal residual of evaluation at all grid points;
  --   testres     maximal residual of evaluation at test points.

  procedure Standard_Linear_Trace_Interpolate
                 ( file : in file_type; sps : in Standard_Sample_List;
                   t : out Standard_Complex_Vectors.Vector;
                   eps,dist,gridres,testres : out double_float );

  -- DESCRIPTION :
  --   Creates the linear trace through the samples.  This linear trace
  --   predicts the sum of the second component of the solution vectors
  --   on the same slice.

  -- ON ENTRY :
  --   file        for intermediate results and diagnostics;
  --   sps         list of starting sample points.

  -- ON RETURN :
  --   t           linear trace, predicts sum of the roots;
  --   eps         accuracy of the samples;
  --   dist        minimial distance between the samples;
  --   gridres     difference between value of trace and sum at grid;
  --   testres     difference between value of trace and sum at test points.

end Irreducible_Component_Creators;
