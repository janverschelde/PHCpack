with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Sample_Points;                       use Sample_Points;
with Sample_Point_Lists;                  use Sample_Point_Lists;
with Interpolation_Points;                use Interpolation_Points;
with Interpolation_Point_Lists;           use Interpolation_Point_Lists;

package Projection_Operators is

-- DESCRIPTION :
--   This package manages linear and central projector operators to
--   project k dimensional solution components to a general k-space.
--   Projectors turn sample points into interpolation points.

  type Standard_Projector is private;
  type Multprec_Projector is private;

-- CREATORS :

  function Create ( hyps : Standard_Complex_VecVecs.VecVec )
                  return Standard_Projector;
  function Create ( hyps : Multprec_Complex_VecVecs.VecVec )
                  return Multprec_Projector;
  function Create ( k,n : natural32 ) return Standard_Projector;
  function Create ( k,n,size : natural32 ) return Multprec_Projector;

  -- DESCRIPTION :
  --   Creates a linear projector that maps points in n-space to k-space
  --   when range of hyps is 1..k and range of hyps(i) is 0..n.
  --   If only the dimensions are given, then a random projector is returned.
  --   The parameter size indicates the size of the multi-precision numbers.

  procedure Update ( p : in out Standard_Projector;
                     basept : in Standard_Sample;
                     basehyp : in Standard_Complex_Vectors.Vector);
  procedure Update ( p : in out Multprec_Projector;
                     basept : in Multprec_Sample;
                     basehyp : in Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Updates the projector with a new base point.  After the first
  --   update, the projector becomes a central projector operator.
  --   The hyperplane basehyp should not contain the base point.
  -- REQUIRED : not Empty(p).

-- SELECTORS :

  function Empty ( p : Standard_Projector ) return boolean;
  function Empty ( p : Multprec_Projector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the projector has been cleared or not been created.

  function Origin_Dimension ( p : Standard_Projector ) return natural32;
  function Origin_Dimension ( p : Multprec_Projector ) return natural32;

  -- DESCRIPTION :
  --   The "origin dimension" of a projector is the dimension of the
  --   points accepted as input to the projector.

  function Target_Dimension ( p : Standard_Projector ) return natural32;
  function Target_Dimension ( p : Multprec_Projector ) return natural32;

  -- DESCRIPTION :
  --   The "target dimension" of a projector is the dimension of the
  --   points on the output of the projector.

  function Hyperplanes ( p : Standard_Projector )
                       return Standard_Complex_VecVecs.VecVec;
  function Hyperplanes ( p : Multprec_Projector )
                       return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the hyperplanes that define the linear operator.

  function Is_Central ( p : Standard_Projector ) return boolean;
  function Is_Central ( p : Multprec_Projector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if p has been updated with base points, false otherwise.

  function Centrality ( p : Standard_Projector ) return natural32;
  function Centrality ( p : Multprec_Projector ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of central points used in the projector.

  function Base_Points ( p : Standard_Projector ) return Standard_Sample_List;
  function Base_Points ( p : Multprec_Projector ) return Multprec_Sample_List;

  -- DESCRIPTION :
  --   Returns the base points that define the central projection operator.

-- PROJECTORS ON BASIC DATA :

  function Project ( p : Standard_Projector;
                     x : Standard_Complex_Vectors.Vector )
                   return Standard_Complex_Vectors.Vector;
  function Project ( p : Standard_Projector;
                     x : Standard_Complex_VecVecs.VecVec )
                   return Standard_Complex_VecVecs.VecVec;
  function Project ( p : Multprec_Projector;
                     x : Multprec_Complex_Vectors.Vector )
                   return Multprec_Complex_Vectors.Vector;
  function Project ( p : Multprec_Projector;
                     x : Standard_Complex_Vectors.Vector )
                   return Multprec_Complex_Vectors.Vector;
  function Project ( p : Multprec_Projector;
                     x : Standard_Complex_VecVecs.VecVec )
                   return Multprec_Complex_VecVecs.VecVec;
  function Project ( p : Multprec_Projector;
                     x : Multprec_Complex_VecVecs.VecVec )
                   return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the linear or central projection of the points.
  --   If the projection operator is empty, then Project returns x.

-- PROJECTORS ON ENCAPSULATED DATA :

  function Project ( p : Standard_Projector; s : Standard_Sample )
                   return Standard_Sample_Node;
  function Project ( p : Multprec_Projector; s : Multprec_Sample )
                   return Multprec_Sample_Node;

  -- DESCRIPTION :
  --   An interpolation point or sample node is a projected sample.

  procedure Project ( p : in Standard_Projector; s : in Standard_Sample_List;
                      sn,sn_last : in out Standard_Sample_Node_List );
  procedure Project ( p : in Multprec_Projector; s : in Multprec_Sample_List;
                      sn,sn_last : in out Multprec_Sample_Node_List );

  -- DESCRIPTION :
  --   The projectors on encapsulated data are essential creators 
  --   for interpolation points, turning sample points into sample nodes.
  --   The projected samples are added to the list sn, with sn_last a
  --   pointer to the last element in the list sn.

-- DESTRUCTORS :

  procedure Shallow_Clear ( p : in out Standard_Projector );
  procedure Shallow_Clear ( p : in out Multprec_Projector );
  procedure Deep_Clear ( p : in out Standard_Projector );
  procedure Deep_Clear ( p : in out Multprec_Projector );

  -- DESCRIPTION :
  --   Releases the occupied memory space.  A shallow clear only
  --   destroys the links, whereas a deep clear also clears the
  --   hyperplanes and base points.

private

  type Standard_Projector_Rep;
  type Standard_Projector is access Standard_Projector_Rep;

  type Multprec_Projector_Rep;
  type Multprec_Projector is access Multprec_Projector_Rep;

end Projection_Operators;
