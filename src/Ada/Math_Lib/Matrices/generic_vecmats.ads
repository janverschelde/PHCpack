with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors,Generic_Matrices;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Matrices is new Generic_Matrices(Ring,Vectors);

package Generic_VecMats is 

  use Matrices;

  type VecMat is array ( integer32 range <> ) of Link_to_Matrix;
  type Link_to_VecMat is access VecMat;

  type VecMat_Array is array ( integer32 range <> ) of Link_to_VecMat;

  procedure Copy ( v : in VecMat; w : in out VecMat );

  -- DESCRIPTION :
  --   Makes a deep copy of all matrices in v to the matrices in w,
  --   after performing a deep clear on w.

  -- REQUIRED : v'range = w'range.

  procedure Clear ( v : in out VecMat );
  procedure Shallow_Clear ( v : in out Link_to_VecMat );
  procedure Deep_Clear ( v : in out Link_to_VecMat );

  -- DESCRIPTION :
  --   A shallow clear deallocates only the pointers.
  --   A deep clear deallocates both pointers and the content.
  --   By default a clear is always deep.

  procedure Clear ( v : in out VecMat_Array );

  -- DESCRIPTION :
  --   Deallocates the space occupied by v.

end Generic_VecMats;
