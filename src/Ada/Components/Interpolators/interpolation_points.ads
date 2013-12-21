with Standard_Complex_Vectors;
with Multprec_Complex_Vectors;
with Sample_Points;                         use Sample_Points;

package Interpolation_Points is

-- DESCRIPTION :
--   Interpolation points are projected sample points.
--   This package provides an encapsulation to manage those points.

  type Standard_Sample_Node is private;
  type Multprec_Sample_Node is private;

-- CREATORS :

  function Create ( spt : Standard_Sample;
                    prospt : Standard_Complex_Vectors.Vector )
                  return Standard_Sample_Node;
  function Create ( spt : Multprec_Sample;
                    prospt : Multprec_Complex_Vectors.Vector )
                  return Multprec_Sample_Node;

  -- DESCRIPTION :
  --    An interpolation point is created from a sample and 
  --    the projected solution vector.

  procedure Update ( snd : in out Standard_Sample_Node;
                     prospt : in Standard_Complex_Vectors.Vector );
  procedure Update ( snd : in out Multprec_Sample_Node;
                     prospt : in Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Updates the node with a new projected vector.

-- SELECTORS :

  function Empty ( sn : Standard_Sample_Node ) return boolean;
  function Empty ( sn : Multprec_Sample_Node ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the node is empty, false otherwise.

  function Sample_Node ( sn : Standard_Sample_Node )
                       return Standard_Complex_Vectors.Vector;
  function Sample_Node ( sn : Multprec_Sample_Node )
                       return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the projected vector in the sample node.
  -- REQUIRED : not Empty(sn).

  function Sample_Point ( sn : Standard_Sample_Node ) return Standard_Sample;
  function Sample_Point ( sn : Multprec_Sample_Node ) return Multprec_Sample;

  -- DESCRIPTION :
  --   Returns the sample from which the node was projected from.
  -- REQUIRED : not Empty(sn).

  function Sample_Vector ( sn : Standard_Sample_Node ) 
                         return Standard_Complex_Vectors.Vector;
  function Sample_Vector ( sn : Multprec_Sample_Node ) 
                         return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector representation of the sample.

-- DESTRUCTORS :

  procedure Shallow_Clear ( sn : in out Standard_Sample_Node );
  procedure Shallow_Clear ( sn : in out Multprec_Sample_Node );

  -- DESCRIPTION :
  --   Clears only the projected samples, not the sample.

  procedure Deep_Clear ( sn : in out Standard_Sample_Node );
  procedure Deep_Clear ( sn : in out Multprec_Sample_Node );

  -- DESCRIPTION :
  --   Deallocation of all occupied memory, also the sample point.

private

  type Standard_Sample_Node_Rep;
  type Standard_Sample_Node is access Standard_Sample_Node_Rep;
  type Multprec_Sample_Node_Rep;
  type Multprec_Sample_Node is access Multprec_Sample_Node_Rep;

end Interpolation_Points;
