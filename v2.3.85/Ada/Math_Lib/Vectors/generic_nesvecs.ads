with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring,Generic_Vectors;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);

package Generic_NesVecs is 

-- DESCRIPTION :
--   A nested vector represents n-dimensional matrices.

-- DATA STRUCTURE :

  use Vectors;

  type NesVec;
  type Link_to_NesVec is access NesVec;
  type Array_of_NesVecs is array ( integer32 range <> ) of Link_to_NesVec;

  type NesVec ( n : natural32; a,b : integer32 ) is record
    case n is
      when 1 => v : Vector(a..b);
      when others => w : Array_of_NesVecs(a..b);
    end case;
  end record;

-- DEEP COPY :

  procedure Copy ( v : in NesVec; w : in out NesVec );
  procedure Copy ( v : in Link_to_NesVec; w : in out Link_to_NesVec );
  procedure Copy ( v : in Array_of_NesVecs; w : in out Array_of_NesVecs );

  -- DESCRIPTION :
  --   Makes a deep copy of v to w.

  -- REQUIRED : v and w have the same dimensions.

-- DESTRUCTORS :
 
  procedure Clear ( v : in out NesVec );
  procedure Clear ( v : in out Link_to_NesVec );
  procedure Clear ( v : in out Array_of_NesVecs );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the nested vector.

end Generic_NesVecs;
