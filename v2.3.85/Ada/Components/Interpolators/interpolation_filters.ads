with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Polynomials;
with Multprec_Complex_Polynomials;
with Sample_Points;                       use Sample_Points;
with Sample_Point_Lists;                  use Sample_Point_Lists;
with Interpolation_Point_Lists;           use Interpolation_Point_Lists;
with Projection_Operators;                use Projection_Operators;

package Interpolation_Filters is

-- DESCRIPTION :
--   This packages manages the data to identify a solution component
--   of a polynomial system.

  type Standard_Filter is private;
  type Multprec_Filter is private;

-- CREATORS :

  function Create ( p : Standard_Projector ) return Standard_Filter;
  function Create ( p : Multprec_Projector ) return Multprec_Filter;

  -- DESCRIPTION :
  --   Returns an encapsulated projector as degree 0 filter.

  procedure Sample_Update
               ( file : in file_type; f : in out Standard_Filter;
                 s : in Standard_Sample_List; d : in natural32 );
  procedure Sample_Update 
               ( file : in file_type; f : in out Multprec_Filter;
                 s : in Multprec_Sample_List; d : in natural32 );

  -- DESCRIPTION :
  --   Updates to a degree d filter with new sample points.
  --   Writes estimate for inverse condition number on file.
  -- REQUIRED :
  --   The filter may not be empty and s must have enough additional points.

  procedure Central_Update
               ( f : in out Standard_Filter; basept : in Standard_Sample;
                 basehyp : in Standard_Complex_Vectors.Vector );
  procedure Central_Update
               ( f : in out Multprec_Filter; basept : in Multprec_Sample;
                 basehyp : in Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   With this update, one central point is added to the projector.
  --   This leads to fewer samples needed to construct the interpolator.
  -- REQUIRED :
  --   The sample point basept does not lie on the hyperplane basehyp.

-- SELECTORS :

  function Dimension ( f : Standard_Filter ) return natural32;
  function Dimension ( f : Multprec_Filter ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the component, 0 if empty.
  --   Note that Dimension(f) = Number_of_Unknowns(Interpolator(f)) - 1.

  function Degree ( f : Standard_Filter ) return natural32;
  function Degree ( f : Multprec_Filter ) return natural32;

  -- DESCRIPTION :
  --   Returns the degree of the interpolating polynomial, 0 if empty.

  function Centrality ( f : Standard_Filter ) return natural32;
  function Centrality ( f : Multprec_Filter ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of central points in the projector.
  --   This number should be added to the degree of the interpolator
  --   to get to the actual degree of the component.

  function Sample_Nodes ( f : Standard_Filter ) 
                        return Standard_Sample_Node_List;
  function Sample_Nodes ( f : Multprec_Filter ) 
                        return Multprec_Sample_Node_List;

  -- DESCRIPTION :
  --   Returns the interpolation points used in the filter.

  function Projector ( f : Standard_Filter ) return Standard_Projector;
  function Projector ( f : Multprec_Filter ) return Multprec_Projector;

  -- DESCRIPTION :
  --   Returns the projection operator used in the filter.

  function Interpolator ( f : Standard_Filter )
                        return Standard_Complex_Polynomials.Poly;
  function Interpolator ( f : Multprec_Filter )
                        return Multprec_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the interpolator representation of the filter.
  --   Note that this does not include the projector.

  function Evaluate ( f : Standard_Filter;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Numbers.Complex_Number;
  function Evaluate ( f : Multprec_Filter;
                      x : Standard_Complex_Vectors.Vector )
                    return Multprec_Complex_Numbers.Complex_Number;
  function Evaluate ( f : Multprec_Filter;
                      x : Multprec_Complex_Vectors.Vector )
                    return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the evaluated form of the interpolator at x.

  function On_Component ( f : Standard_Filter; tol : double_float;
                          x : Standard_Complex_Vectors.Vector ) 
                        return boolean;
  function On_Component ( file : file_type;
                          f : Standard_Filter; tol : double_float;
                          x : Standard_Complex_Vectors.Vector ) 
                        return boolean;
  function On_Component ( f : Multprec_Filter; tol : double_float;
                          x : Standard_Complex_Vectors.Vector ) 
                        return boolean;
  function On_Component ( file : file_type;
                          f : Multprec_Filter; tol : double_float;
                          x : Standard_Complex_Vectors.Vector ) 
                        return boolean;
  function On_Component ( f : Multprec_Filter; tol : double_float;
                          x : Multprec_Complex_Vectors.Vector ) 
                        return boolean;
  function On_Component ( file : file_type;
                          f : Multprec_Filter; tol : double_float;
                          x : Multprec_Complex_Vectors.Vector ) 
                        return boolean;

  -- DESCRIPTION :
  --   Returns true if the point belongs to the component, when
  --   the evaluation in the filter is within the given tolerance.
  --   If the filter is empty, then true is returned as well.
  --   The residual is written to the file when the file is provided.

-- DESTRUCTORS :

  procedure Shallow_Clear ( f : in out Standard_Filter );
  procedure Shallow_Clear ( f : in out Multprec_Filter );

  -- DESCRIPTION :
  --   Destroys only the interpolator, not the interpolation points.

  procedure Deep_Clear ( f : in out Standard_Filter );
  procedure Deep_Clear ( f : in out Multprec_Filter );

  -- DESCRIPTION :
  --   Destroys interpolation points and filtering polynomial.

private

  type Standard_Filter_Rep;
  type Standard_Filter is access Standard_Filter_Rep;
  type Multprec_Filter_Rep;
  type Multprec_Filter is access Multprec_Filter_Rep;

end Interpolation_Filters;
