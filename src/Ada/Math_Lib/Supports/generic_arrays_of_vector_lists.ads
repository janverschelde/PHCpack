with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors;
with Generic_VecVecs;
with Generic_Lists_of_Vectors;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);
  with package Lists is new Generic_Lists_of_Vectors(Ring,Vectors,VecVecs);

package Generic_Arrays_of_Vector_Lists is

-- DESCRIPTION :
--   This package defines arrays of lists of links to vectors.

  use VecVecs,Lists;

-- DATA STRUCTURES :

  type Array_of_Lists is array ( integer32 range <> ) of List;
  type Link_to_Array_of_Lists is access Array_of_Lists;

-- CREATORS :

  function Deep_Create    ( v : Array_of_VecVecs ) return Array_of_Lists;
  function Shallow_Create ( v : Array_of_VecVecs ) return Array_of_Lists;

  function Deep_Create    ( L : Array_of_Lists ) return Array_of_VecVecs;
  function Shallow_Create ( L : Array_of_Lists ) return Array_of_VecVecs;

  -- DESCRIPTION :
  --   L := Create(v) equals v := Create(L).
  --   There is no sharing of pointers with a deep Create.
  --   With a shallow Create, both structure share the pointers.


-- COMPARISON and COPYING :

  function Equal ( L1,L2 : Array_of_Lists ) return boolean;
  function Equal ( L1,L2 : Link_to_Array_of_Lists ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both arrays have the same lists.

  procedure Copy ( L1 : in Array_of_Lists; L2 : in out Array_of_Lists );

  -- DESCRIPTION :
  --   After Copy(L1,L2), Equal(L1,L2) holds.
  --   Of course, this is a deep copy, in constrast to L2 := L1.

-- SELECTOR :

  function Length_Of ( L : Array_of_Lists ) return natural32;

  -- DESCRIPTION :
  --   Returns the sum of all lengths of the lists in L.

-- DESTRUCTORS :

  procedure Deep_Clear    ( L : in out Array_of_Lists );
  procedure Shallow_Clear ( L : in out Array_of_Lists );
  procedure Deep_Clear    ( L : in out Link_to_Array_of_Lists );
  procedure Shallow_Clear ( L : in out Link_to_Array_of_Lists );

  -- DESCRIPTION :
  --   Frees allocated memory space.
  --   With a deep clear, also the content of the lists are cleared,
  --   while with a shallow clear, only the lists structures will be
  --   destroyed, the points in the lists will remain.

end Generic_Arrays_of_Vector_Lists;
