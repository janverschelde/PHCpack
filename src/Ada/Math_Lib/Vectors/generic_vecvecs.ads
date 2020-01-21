with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Abstract_Ring;
with Generic_Vectors;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);

package Generic_VecVecs is 

  use Vectors;

  type VecVec is array ( integer32 range <> ) of Link_to_Vector;
  type Link_to_VecVec is access VecVec;
  type Array_of_VecVecs is array ( integer32 range <> ) of Link_to_VecVec;
  type Link_to_Array_of_VecVecs is access Array_of_VecVecs;

  function Create_Copy ( v : VecVec ) return VecVec;
  function Create_Copy ( v : Array_of_VecVecs ) return Array_of_VecVecs;

  -- DESCRIPTION :
  --   The object on return is a copy of the given v.

  procedure Copy ( v : in VecVec; w : in out VecVec );

  -- DESCRIPTION :
  --   Makes a deep copy of v into w, after clearing w first.

  -- REQUIRED : v'range = w'range.

  procedure Clear ( v : in out VecVec );
  procedure Clear ( v : in out Array_of_VecVecs );
  procedure Shallow_Clear ( v : in out Link_to_VecVec );
  procedure Shallow_Clear ( v : in out Link_to_Array_of_VecVecs );
  procedure Deep_Clear ( v : in out Link_to_VecVec );
  procedure Deep_Clear ( v : in out Link_to_Array_of_VecVecs );

  -- DESCRIPTION :
  --   A shallow clear deallocates the pointers, whereas a deep clear
  --   also releases the content.  By default a clear is always deep.

end Generic_VecVecs;
