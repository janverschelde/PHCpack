with Abstract_Ring;
with Generic_Vectors;
with Generic_VecVecs;
with Generic_Lists;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);

package Generic_Lists_of_Vectors is

-- DESCRIPTION :
--   An implementation of lists of links to vectors.

  use Vectors,VecVecs;

-- DATA STRUCTURE : a list of pointers to vectors

  package Vector_Lists is new Generic_Lists(Link_to_Vector);
  type List is new Vector_Lists.List;

-- CREATORS :

  function Deep_Create    ( v : VecVec ) return List;
  function Shallow_Create ( v : VecVec ) return List;

  function Deep_Create    ( L : List ) return VecVec;
  function Shallow_Create ( L : List ) return VecVec;

  -- DESCRIPTION :
  --   L := Create(v) equals v := Create(l).
  --   There is no sharing of pointers with a deep Create.
  --   With a shallow Create, both structure share the pointers.

-- COMPARISON and COPYING :

  procedure Copy ( L1 : in List; L2 : in out List );

  -- DESCRIPTION :
  --   After Copy(L1,L2), Equal(L1,L2) holds.
  --   Of course, this is a deep copy, a shallow copy is given by L2 := L1.

  function Equal ( L1,L2 : List ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both lists have the same vectors, false otherwise.

-- SELECTORS :

  function Is_In ( L : List; v : Vector ) return boolean;
  function Is_In ( L : List; v : Link_to_Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the vector belongs to l, false otherwise.

  function Sub_List ( L1,L2 : List ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all elements in L1 occur in L2, false otherwise.

-- CONSTRUCTORS :

  procedure Append ( first,last : in out List; v : in Vector );

  -- DESCRIPTION :
  --   The vector will be appended to the list first,
  --   last is a pointer to the last element of the list first.

  procedure Append_Diff ( first,last : in out List; v : in Vector );
  procedure Append_Diff ( first,last : in out List; v : in Link_to_Vector );

  -- DESCRIPTION :
  --   Only when v does not already belong to first, v will be added.

  procedure Deep_Concat    ( first,last : in out List; L : in List );
  procedure Shallow_Concat ( first,last : in out List; L : in List );

  -- DESCRIPTION :
  --   The list l will be concatenated to the list first,
  --   last is a pointer to the last element of the list first.
  --   With a deep concatenation, no pointers are shared.

  procedure Deep_Concat_Diff    ( first,last : in out List; L : in List );
  procedure Shallow_Concat_Diff ( first,last : in out List; L : in List );

  -- DESCRIPTION :
  --   Only those vectors of l will be concatenated that are not already
  --   in the list first.
  --   With a deep concatenation, no pointers are shared.

  procedure Remove ( L : in out List; x : in Vector );
  procedure Remove ( L : in out List; x : in Link_to_Vector );

  -- DESCRIPTION :
  --   Removes the point x from the list l.

  procedure Swap_to_Front ( L : in out List; x : in Vector );
  procedure Swap_to_Front ( L : in out List; x : in Link_to_Vector );

  -- DESCRIPTION :
  --   The point x belongs to the list l,
  --   Its content will be place in front of l and the first element
  --   of l will be swapped to the place of x.

-- DESTRUCTORS :

  procedure Deep_Clear    ( L : in out List );
  procedure Shallow_Clear ( L : in out List );

  -- DESCRIPTION :
  --   Frees all allocated memory space.
  --   A deep clear deallocates also the points to the integer vectors,
  --   while a shallow clear only removes the list structure.

end Generic_Lists_of_Vectors;
