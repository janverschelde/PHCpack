with generic_lists;
with Standard_Complex_VecVecs;
with Multprec_Complex_VecVecs;
with Sample_Point_Lists;                    use Sample_Point_Lists;
with Interpolation_Points;                  use Interpolation_Points;

package Interpolation_Point_Lists is

-- DESCRIPTION :
--   Lists of interpolation points are projected sample point lists.
--   This package provides an encapsulation to manage those points.

  package Lists_of_Standard_Sample_Nodes is
    new generic_lists(Standard_Sample_Node);
  type Standard_Sample_Node_List is new Lists_of_Standard_Sample_Nodes.List;

  package Lists_of_Multprec_Sample_Nodes is
    new generic_lists(Multprec_Sample_Node);
  type Multprec_Sample_Node_List is new Lists_of_Multprec_Sample_Nodes.List;

-- CREATORS :

  procedure Create ( spl : in Standard_Sample_List;
                     prospl : in Standard_Complex_VecVecs.VecVec;
                     first,last : in out Standard_Sample_Node_List );
  procedure Create ( spl : in Multprec_Sample_List;
                     prospl : in Multprec_Complex_VecVecs.VecVec;
                     first,last : in out Multprec_Sample_Node_List );

  -- DESCRIPTION :
  --   Returns a list of interpolation points from a list of samples
  --   and their projections.
  -- REQUIRED : prospl'range = 1..Length_Of(spl)

-- SELECTORS :

-- DESTRUCTORS :

  procedure Shallow_Clear ( l : in out Standard_Sample_Node_List );
  procedure Shallow_Clear ( l : in out Multprec_Sample_Node_List );

  -- DESCRIPTION :
  --   Clears only the projected samples, not the samples.

  procedure Deep_Clear ( l : in out Standard_Sample_Node_List );
  procedure Deep_Clear ( l : in out Multprec_Sample_Node_List );

  -- DESCRIPTION :
  --   Deallocation of all occupied memory, also the sample point.

end Interpolation_Point_Lists;
