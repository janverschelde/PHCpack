with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Multprec_Complex_Vectors;
with Sample_Points;                    use Sample_Points;
with Sample_Point_Lists;               use Sample_Point_Lists;
with Span_of_Component;                use Span_of_Component;
with Interpolation_Filters;            use Interpolation_Filters;

package Irreducible_Components is

-- DESCRIPTION :
--   This package offers a data structure to represent irreducible
--   solution components of polynomial systems.
--   There are three parts in this representation :
--     1) an interpolation filter = projector + polynomial;
--     2) the span of the component = linear equations;
--     3) a set of labels to generic points on the component.

  type Standard_Irreducible_Component is private;
  type Multprec_Irreducible_Component is private;

-- CREATORS :

  function Create ( f : Standard_Filter )
                  return Standard_Irreducible_Component;
  function Create ( f : Standard_Filter; s : Standard_Span )
                  return Standard_Irreducible_Component;
  function Create ( f : Multprec_Filter )
                  return Multprec_Irreducible_Component;
  function Create ( f : Multprec_Filter; s : Multprec_Span )
                  return Multprec_Irreducible_Component;

  -- DESCRIPTION :
  --   Creates an irreducible component from the filter and
  --   optionally from the space spanned by the component.

  procedure Initialize_Labels ( c : in out Standard_Irreducible_Component;
                                d : in integer32 );
  procedure Initialize_Labels ( c : in out Standard_Irreducible_Component;
                                lab : in Standard_Natural_Vectors.Vector );
  procedure Initialize_Labels ( c : in out Multprec_Irreducible_Component;
                                d : in integer32 );
  procedure Initialize_Labels ( c : in out Multprec_Irreducible_Component;
                                lab : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Initializes the vector of labels of generic points,
  --   where d ought to correspond with the degree of the component.
  --   The component is initialized with the labels if given as parameter.

  -- REQUIRED : Labels are sorted in increasing order.

  procedure Add_Label ( c : in out Standard_Irreducible_Component;
                        i : in natural32; ind : out integer32);
  procedure Add_Label ( c : in out Multprec_Irreducible_Component;
                        i : in natural32; ind : out integer32 );

  -- DESCRIPTION :
  --   When the i-th generic point is found to lie on c,
  --   then its label should be added with Add_Label(c,i).
  --   The number on return is the index of i as a label to c.

  procedure Add_Point ( c : in out Standard_Irreducible_Component;
                        spt : in Standard_Sample );
  procedure Add_Point ( c : in out Multprec_Irreducible_Component;
                        spt : in Standard_Sample );

  -- DESCRIPTION :
  --   This is the companion operation to "Add_Label": the generic
  --   point found to lie on c is added to the component.

  procedure Add_Points ( c : in out Standard_Irreducible_Component;
                         sps,sps_last : in Standard_Sample_List );
  procedure Add_Points ( c : in out Multprec_Irreducible_Component;
                         sps,sps_last : in Standard_Sample_List );

  -- DESCRIPTION :
  --   Adds all Degree(c) samples all at once to the component, where sps
  --   points to the first and sps_last to the last element in the list.
  --   Note that there is sharing of data.

  procedure Select_Labeled_Points
               ( c : in out Standard_Irreducible_Component;
                 sps : in Standard_Sample_List );
  procedure Select_Labeled_Points
               ( c : in out Multprec_Irreducible_Component;
                 sps : in Standard_Sample_List );

  -- DESCRIPTION :
  --   Selects those points from the list sps that have been initialized
  --   as labels in the component.

  procedure Add_Filter ( c : in out Standard_Irreducible_Component;
                         f : in Standard_Filter );

  -- DESCRIPTION :
  --   Adds the filter to the internal representation of c.

-- SELECTORS :

  function Dimension ( c : Standard_Irreducible_Component ) return natural32;
  function Dimension ( c : Multprec_Irreducible_Component ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the irreducible component.

  function Degree ( c : Standard_Irreducible_Component ) return natural32;
  function Degree ( c : Multprec_Irreducible_Component ) return natural32;

  -- DESCRIPTION :
  --   Returns the degree of the irreducible component.
  --   Note that this is the degree of the interpolating filter plus
  --   the number of central points used in the projection.

  function Filter ( c : Standard_Irreducible_Component )
                  return Standard_Filter;
  function Filter ( c : Multprec_Irreducible_Component )
                  return Multprec_Filter;

  -- DESCRIPTION :
  --   Returns the filter representing the component.

  function Span ( c : Standard_Irreducible_Component ) return Standard_Span;
  function Span ( c : Multprec_Irreducible_Component ) return Multprec_Span;

  -- DESCRIPTION :
  --   Returns the space spanned by the component.

  function Labels ( c : Standard_Irreducible_Component )
                  return Standard_Natural_Vectors.Link_to_Vector;
  function Labels ( c : Multprec_Irreducible_Component )
                  return Standard_Natural_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns labels to the generic points that have been found to
  --   lie on the component.

  function Points ( c : Standard_Irreducible_Component )
                  return Standard_Sample_List;
  function Points ( c : Multprec_Irreducible_Component )
                  return Standard_Sample_List;

  -- DESCRIPTION :
  --   Returns the pointer to the list of generic points on the component.

-- MEMBERSHIP TESTS :

  function On_Component ( c : Standard_Irreducible_Component;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return boolean;
  function On_Component ( file : file_type;
                          c : Standard_Irreducible_Component;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return boolean;
  function On_Component ( c : Multprec_Irreducible_Component;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return boolean;
  function On_Component ( file : file_type;
                          c : Multprec_Irreducible_Component;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return boolean;
  function On_Component ( c : Multprec_Irreducible_Component;
                          x : Multprec_Complex_Vectors.Vector;
                          tol : double_float ) return boolean;
  function On_Component ( file : file_type;
                          c : Multprec_Irreducible_Component;
                          x : Multprec_Complex_Vectors.Vector;
                          tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the point x belongs to the component.
  --   The optional file is for intermediate output and diagnostics.

  function Homotopy_Membership_Test
              ( c : Standard_Irreducible_Component;
                x : Standard_Complex_Vectors.Vector; tol : double_float )
              return boolean;
  function Homotopy_Membership_Test
              ( file : file_type; c : Standard_Irreducible_Component;
                x : Standard_Complex_Vectors.Vector; tol : double_float )
              return boolean;

  -- DESCRIPTION :
  --   Performs a homotopy membership test for the point x to belong to c.
  -- REQUIRED :
  --   not Is_Null(Points(c)) and Sampling Machine tuned and initialized.

-- DESTRUCTORS :

  procedure Clear ( c : in out Standard_Irreducible_Component );
  procedure Clear ( c : in out Multprec_Irreducible_Component );

  -- DESCRIPTION :
  --   Memory is freed and the objects are destroyed.

private

  type Standard_Irreducible_Component is record
    f : Standard_Filter;
    s : Standard_Span;
    g : Standard_Natural_Vectors.Link_to_Vector;
    p,plast : Standard_Sample_List;
  end record;

  type Multprec_Irreducible_Component is record
    f : Multprec_Filter;
    s : Multprec_Span;
    g : Standard_Natural_Vectors.Link_to_Vector;
    p,plast : Standard_Sample_List;
  end record;

end Irreducible_Components;
