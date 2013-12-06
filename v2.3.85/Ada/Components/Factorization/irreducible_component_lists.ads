with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Matrices;
with Standard_Complex_Vectors;
with Multprec_Complex_Vectors;
with generic_lists;
with Irreducible_Components;             use Irreducible_Components;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Sample_Point_Grids;                 use Sample_Point_Grids;

package Irreducible_Component_Lists is

-- DESCRIPTION :
--   This package provides data types and interpolators to represent
--   the irreducible components of a pure dimensional solution set.
--   We represent a component by a projector and an interpolator.
--   There are two major creator types :
--     1) incremental interpolation, for increasing degrees;
--     2) monodromy followed by some kind of Newton interpolation.

-- DATA STRUCTURES :

  package Lists_of_Standard_Irreducible_Components is
    new Generic_Lists(Standard_Irreducible_Component);
  type Standard_Irreducible_Component_List is
    new Lists_of_Standard_Irreducible_Components.List;

  Standard_Null_List : constant Standard_Irreducible_Component_List
    := Standard_Irreducible_Component_List
         (Lists_of_Standard_Irreducible_Components.Null_List);

  package Lists_of_Multprec_Irreducible_Components is
    new Generic_Lists(Multprec_Irreducible_Component);
  type Multprec_Irreducible_Component_List is
    new Lists_of_Multprec_Irreducible_Components.List;

  Multprec_Null_List : constant Multprec_Irreducible_Component_List
    := Multprec_Irreducible_Component_List
         (Lists_of_Multprec_Irreducible_Components.Null_List);

-- CREATORS WITH INCREMENTAL INTERPOLATION :

  procedure Standard_Massive_Interpolate
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 stoptol,membtol : in double_float;
                 deco,deco_last : out Standard_Irreducible_Component_List );
  procedure Standard_Incremental_Interpolate
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 stoptol,membtol : in double_float;
                 deco,deco_last : out Standard_Irreducible_Component_List );
  procedure Standard_Incremental_Interpolate_with_Span
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 stoptol,membtol : in double_float;
                 deco,deco_last : out Standard_Irreducible_Component_List );
  procedure Standard_Incremental_Central_Interpolate
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 stoptol,membtol : in double_float;
                 deco,deco_last : out Standard_Irreducible_Component_List );

  procedure Multprec_Massive_Interpolate
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 size : in natural32; stoptol,membtol : in double_float;
                 deco,deco_last : out Multprec_Irreducible_Component_List );
  procedure Multprec_Incremental_Interpolate
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 size : in natural32; stoptol,membtol : in double_float;
                 deco,deco_last : out Multprec_Irreducible_Component_List );
  procedure Multprec_Incremental_Interpolate_with_Span
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 size : in natural32; stoptol,membtol : in double_float;
                 deco,deco_last : out Multprec_Irreducible_Component_List );
  procedure Multprec_Incremental_Central_Interpolate
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 size : in natural32; stoptol,membtol : in double_float;
                 deco,deco_last : out Multprec_Irreducible_Component_List );

  -- DESCRIPTION :
  --   Creates a list of irreducible components from the list of samples.
  --   The "Massive" version creates as many interpolators as the degree,
  --   whereas "Incremental" versions check first whether points satisfy
  --   the already created filters before interpolating.
  --   The "with_Span" routines also determine the linear span of components.
  --   The "_Central_" routines use central projections.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   full_output indicates whether all sample diagnostics are needed,
  --             or whether only a summary will be written;
  --   sps       list of generic points on the solution set;
  --   size      size of the multi-precision numbers;
  --   stoptol   tolerance to stop interpolation;
  --   membtol   tolerance to decide membership.

  -- ON RETURN :
  --   deco      list of irreducible components;
  --   deco_last points to the last element of the list deco.

-- CREATORS TO PREDICT BREAKUP WITH MONODROMY LOOPS :

  procedure Monodromy_Breakup
               ( file : in file_type; sps : in Standard_Sample_List;
                 threshold : in natural32; tol : in double_float;
                 deco,deco_last : out Standard_Irreducible_Component_List;
                 nit : out natural32 );

  procedure New_Monodromy_Breakup
               ( file : in file_type; sps : in Standard_Sample_List;
                 threshold : in natural32; tol : in double_float;
                 deco,deco_last : out Standard_Irreducible_Component_List;
                 nit : out natural32 );

  procedure Monodromy_Breakup
               ( file : in file_type; sps : in Standard_Sample_List;
                 threshold : in natural32; tol : in double_float;
                 deco,deco_last : out Standard_Irreducible_Component_List;
                 nit : out natural32;
                 grid,grid_last : in out Standard_Sample_Grid );

  -- DESCRIPTION :
  --   Applies the monodromy group actions to break up a solution set.
  --   The second routine accumulates the lists of sample points used.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   sps       list of generic points;
  --   threshold is the stabilizing treshold to stop the iterations;

  -- ON RETURN :
  --   deco      list of irreducible components, contains labels to
  --             generic points and the points on the same component;
  --   deco_last points to the last element of the list deco;
  --   nit       total number of iterations;
  --   grid      list of all sample points used in the breakup;
  --   grid_last points to the last list in grid.

  procedure Distribute_Points
               ( deco : in out Standard_Irreducible_Component_List;
                 sps : in Standard_Sample_List );
  procedure Distribute_Points
               ( deco : in out Multprec_Irreducible_Component_List;
                 sps : in Standard_Sample_List );

  -- DESCRIPTION :
  --   Distributes the generic points in the sample list sps along the
  --   components in the list deco.  This operation only makes sense
  --   after the labels have been read from file into the list.

-- CREATORS WITH BOOTSTRAPPING NEWTON :

  procedure Standard_Newton_Interpolate
               ( file : in file_type;
                 deco : in out Standard_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix );
  procedure Multprec_Newton_Interpolate
               ( file : in file_type; size : in natural32;
                 deco : in out Multprec_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix );

  -- DESCRIPTION :
  --   Creates filters by Newton interpolation.

  -- REQUIRED : The sampling machine is initialized and properly tuned.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   size      size of the multi-precision numbers;
  --   deco      list of components, with labels and generic points for
  --             every component in the list.

  -- ON RETURN :
  --   deco      updated with filters;
  --   numres    i-th row contains numerical results of i-th component:
  --              (i,1) : degree of component;
  --              (i,2) : accuracy of the samples;
  --              (i,3) : minimal distance between samples;
  --              (i,4) : maximal residual at all grid points;
  --              (i,5) : maximal residual at test points.

-- CREATORS WITH TRACE FORMS :

  procedure Standard_Trace_Interpolate
               ( file : in file_type;
                 deco : in out Standard_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix );
  procedure Multprec_Trace_Interpolate
               ( file : in file_type; size : in natural32;
                 deco : in out Multprec_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix );

  -- DESCRIPTION :
  --   Creates filters using the trace form of the interpolator.
  --   The parameters have the same meaning as with Newton interpolation.

  procedure Standard_Power_Trace_Interpolate
               ( file : in file_type;
                 deco : in out Standard_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix );
  procedure Multprec_Power_Trace_Interpolate
               ( file : in file_type; size : in natural32;
                 deco : in out Multprec_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix );

  -- DESCRIPTION :
  --   Creates filters using Newton identities in the trace form.
  --   Parameters have same meaning as above.

  procedure Standard_Linear_Trace_Interpolate
               ( file : in file_type;
                 deco : in out Standard_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix );

  -- DESCRIPTION :
  --   Creates only the linear trace for each component.
  --   The parameters have the same role as above.

-- SELECTORS :

  function On_Component ( L : Standard_Irreducible_Component_List;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return natural32;
  function On_Component ( file : file_type;
                          L : Standard_Irreducible_Component_List;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return natural32;
  function On_Component ( L : Multprec_Irreducible_Component_List;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return natural32;
  function On_Component ( file : file_type;
                          L : Multprec_Irreducible_Component_List;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return natural32;
  function On_Component ( L : Multprec_Irreducible_Component_List;
                          x : Multprec_Complex_Vectors.Vector;
                          tol : double_float ) return natural32;
  function On_Component ( file : file_type;
                          L : Multprec_Irreducible_Component_List;
                          x : Multprec_Complex_Vectors.Vector;
                          tol : double_float ) return natural32;

  -- DESCRIPTION :
  --   If the point x belongs to one of the component in the list,
  --   then the index of the container component in the list is
  --   returned, otherwise the number of return is zero.
  --   The optional file is for intermediate output and diagnostics.

  function Homotopy_Filter ( L : Standard_Irreducible_Component_List;
                             x : Standard_Complex_Vectors.Vector;
                             tol : double_float ) return natural32;
  function Homotopy_Filter ( file : file_type;
                             L : Standard_Irreducible_Component_List;
                             x : Standard_Complex_Vectors.Vector;
                             tol : double_float ) return natural32;

  -- DESCRIPTION :
  --   These membership tests rely on homotopy continuation,
  --   and ought to be used in conjunction with the monodromy breakup.
  -- REQUIRED : Sampling Machine is initialized and properly tuned.

end Irreducible_Component_Lists;
