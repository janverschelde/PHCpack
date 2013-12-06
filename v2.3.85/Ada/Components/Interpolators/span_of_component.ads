with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;

package Span_of_Component is

-- DESCRIPTION :
--   This packages provides an abstraction for the space spanned by
--   a solution component of a polynomial system.

  type Standard_Span is private;
  type Multprec_Span is private;

  Null_Standard_Span : constant Standard_Span;
  Null_Multprec_Span : constant Multprec_Span;

-- CREATORS :

  function Create ( s : Standard_Sample_List; tol : double_float )
                  return Standard_Span;
  function Create ( s : Multprec_Sample_List; size : natural32;
                    tol : double_float ) return Multprec_Span;

  -- DESCRIPTION :
  --   Creates the span from a list s of sample points.
  --   The parameter size is the size of the multi-precision numbers.
  --   The parameter tol is used to decide whether a number is zero.
  --   The span on return may be empty if either the points span the
  --   entire ambient space or if there are insufficient points given
  --   to determine the dimension of the space.

  function Create ( n,d : natural32; frv : Standard_Integer_Vectors.Vector;
                    equ : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Standard_Span;
  function Create ( n,d : natural32; frv : Standard_Integer_Vectors.Vector;
                    equ : Multprec_Complex_Poly_Systems.Poly_Sys )
                  return Multprec_Span;

  -- DESCRIPTION :
  --   Returns the span defined by the parameters.

  -- ON ENTRY :
  --   n          dimension of the ambient space;
  --   d          dimension of the space;
  --   frv        free variables, i.e.: d naturals in 1..n;
  --   equ        n-d equations in n variables defining the space.

-- SELECTORS :

  function Empty ( sp : Standard_Span ) return boolean;
  function Empty ( sp : Multprec_Span ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the space is undetermined.

  function Ambient_Dimension ( sp : Standard_Span ) return natural32;
  function Ambient_Dimension ( sp : Multprec_Span ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the ambient space.
  --   Returns zero if Empty(sp).

  function Dimension ( sp : Standard_Span ) return natural32;
  function Dimension ( sp : Multprec_Span ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the space spanned by the component.
  --   Returns zero if Empty(sp).

  function Free_Variables ( sp : Standard_Span )
                          return Standard_Integer_Vectors.Vector;
  function Free_Variables ( sp : Multprec_Span )
                          return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of d naturals in the range 1..n,
  --   where n = Ambient_Dimension(sp) and d = Dimension(sp).
  -- REQUIRED : not Empty(sp).

  function Equations ( sp : Standard_Span )
                     return Standard_Complex_Poly_Systems.Poly_Sys;
  function Equations ( sp : Multprec_Span )
                     return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns n-d equations that cut out the space,
  --   where n = Ambient_Dimension(sp) and d = Dimension(sp).
  --   The polynomials have n variables.
  -- REQUIRED : not Empty(sp).

  function In_Span ( sp : Standard_Span; tol : double_float;
                     x : Standard_Complex_Vectors.Vector ) return boolean;
  function In_Span ( file : file_type;
                     sp : Standard_Span; tol : double_float;
                     x : Standard_Complex_Vectors.Vector ) return boolean;
  function In_Span ( sp : Multprec_Span; tol : double_float;
                     x : Standard_Complex_Vectors.Vector ) return boolean;
  function In_Span ( file : file_type;
                     sp : Multprec_Span; tol : double_float;
                     x : Standard_Complex_Vectors.Vector ) return boolean;
  function In_Span ( sp : Multprec_Span; tol : double_float;
                     x : Multprec_Complex_Vectors.Vector ) return boolean;
  function In_Span ( file : file_type;
                     sp : Multprec_Span; tol : double_float;
                     x : Multprec_Complex_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Decides whether the point x belongs to the span of the component.
  --   Returns true if the residual of the point x evaluated in the equations
  --   is within the given tolerance.  Otherwise false is returned.
  --   If the span is empty, then true is returned.
  --   The residual is written on file if the file is input parameter.

-- SUBSPACE RESTRICTIONS :

  function Restrict_Hyperplane
              ( sp : Standard_Span; L : natural32; tol : double_float;
                hyp : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector;
  function Restrict_Hyperplane
              ( sp : Standard_Span; L : natural32; tol : double_float;
                hyp : Standard_Complex_VecVecs.VecVec )
              return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Restricts the hyperplane(s) hyp to the span.

  function Restrict ( sp : Standard_Span; L : in natural32; tol : double_float;
                      p : Standard_Complex_Poly_Systems.Poly_Sys )
                    return Standard_Complex_Poly_Systems.Poly_Sys;
  function Restrict ( sp : Multprec_Span; L : in natural32; tol : double_float;
                      p : Multprec_Complex_Poly_Systems.Poly_Sys )
                    return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the restriction of the polynomial system p to the span.

  -- ON ENTRY :
  --   sp        span of component;
  --   l         level, number of added sections on p;
  --   tol       tolerance to decide whether number is zero;
  --   p         polynomial system to restrict to sp.

  function Restrict ( sp : Standard_Span; L : in natural32;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector;
  function Restrict ( sp : Multprec_Span; L : in natural32;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector;
  function Restrict ( sp : Multprec_Span; L : in natural32;
                      x : Multprec_Complex_Vectors.Vector )
                    return Multprec_Complex_Vectors.Vector;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      sol : Standard_Complex_Solutions.Solution )
                    return Standard_Complex_Solutions.Solution;
  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sol : Standard_Complex_Solutions.Solution )
                    return Standard_Complex_Solutions.Solution;
  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sol : Multprec_Complex_Solutions.Solution )
                    return Multprec_Complex_Solutions.Solution;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      spt : Standard_Sample ) return Standard_Sample;
  function Restrict ( sp : Multprec_Span; L : in natural32;
                      spt : Standard_Sample ) return Standard_Sample;
  function Restrict ( sp : Multprec_Span; L : in natural32;
                      spt : Multprec_Sample ) return Multprec_Sample;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      sps : Standard_Sample_List ) return Standard_Sample_List;
  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sps : Multprec_Sample_List ) return Multprec_Sample_List;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      sol : Standard_Complex_Solutions.Solution_List )
                    return Standard_Complex_Solutions.Solution_List;
  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sol : Standard_Complex_Solutions.Solution_List )
                    return Standard_Complex_Solutions.Solution_List;
  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sol : Multprec_Complex_Solutions.Solution_List )
                    return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the restriction of the solutions to the free variables,
  --   as defined by the span.  The last l variables are copied.

-- DESTRUCTORS :

  procedure Clear ( sp : in out Standard_Span );
  procedure Clear ( sp : in out Multprec_Span );

  -- DESCRIPTION :
  --   Deallocation of the occupied memory space.

private

  type Standard_Span_Rep;
  type Standard_Span is access Standard_Span_Rep;

  Null_Standard_Span : constant Standard_Span := null;

  type Multprec_Span_Rep;
  type Multprec_Span is access Multprec_Span_Rep;

  Null_Multprec_Span : constant Multprec_Span := null;

end Span_of_Component;
