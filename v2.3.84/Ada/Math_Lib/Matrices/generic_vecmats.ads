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

  procedure Clear ( v : in out VecMat );
  procedure Shallow_Clear ( v : in out Link_to_VecMat );
  procedure Deep_Clear ( v : in out Link_to_VecMat );

  -- DESCRIPTION :
  --   A shallow clear on deallocates the pointers, whereas a deep clear
  --   also releases the content.  By default a clear is always deep.

end Generic_VecMats;
