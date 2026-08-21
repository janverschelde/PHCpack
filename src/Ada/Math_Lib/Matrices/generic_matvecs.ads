with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);

package Generic_MatVecs is 

-- DESCRIPTION :
--   A MatVec data structure is a matrix of vectors.

  use Vectors;

  type MatVec is
    array ( integer32 range <>, integer32 range <> ) of Link_to_Vector;

  type Link_to_MatVec is access MatVec;

  type MatVec_Array is array ( integer32 range <> ) of Link_to_MatVec;

  procedure Copy ( A : in MatVec; B : in out MatVec );

  -- DESCRIPTION :
  --   Makes a deep copy of all vectors in A to the Vectors in B,
  --   after performing a deep clear on B.

  -- REQUIRED : A'range(1) = B'range(1) and A'range(2) = B'range(2).

  procedure Clear ( A : in out MatVec );
  procedure Shallow_Clear ( A : in out Link_to_MatVec );
  procedure Deep_Clear ( A : in out Link_to_MatVec );

  -- DESCRIPTION :
  --   A shallow clear deallocates only the pointers.
  --   A deep clear deallocates both pointers and the content.
  --   By default a clear is always deep.

  procedure Clear ( A : in out MatVec_Array );

  -- DESCRIPTION :
  --   Deallocates the space occupied by A.

end Generic_MatVecs;
