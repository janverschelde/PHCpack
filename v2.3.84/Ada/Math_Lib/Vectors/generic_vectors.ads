with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Abstract_Ring;

generic

  with package Ring is new Abstract_Ring(<>);

package Generic_Vectors is

-- DESCRIPTION :
--   An abstraction for vectors with coefficients over any ring.
--   The same functionality is provided for pointers to vectors.

  use Ring;

  type Vector is array ( integer32 range <> ) of number;
  type Link_to_Vector is access Vector;

-- COMPARISON AND COPYING :

  function Equal ( v1,v2 : Vector ) return boolean;
  function Equal ( v1,v2 : Link_to_Vector ) return boolean;

  procedure Copy ( v1 : in Vector; v2 : in out Vector );
  procedure Copy ( v1 : in Link_to_Vector; v2 : in out Link_to_Vector );

-- ARITHMETIC OPERATIONS AS FUNCTIONS :

  function "+" ( v1,v2 : Vector ) return Vector;            -- return v1+v2;
  function "+" ( v1,v2 : Link_to_Vector ) return Link_to_Vector;
  function "+" ( v : Vector ) return Vector;                -- return +v;
  function "+" ( v : Link_to_Vector ) return Link_to_Vector;
  function "-" ( v : Vector ) return Vector;                -- return -v;
  function "-" ( v : Link_to_Vector ) return Link_to_Vector;
  function "-" ( v1,v2 : Vector ) return Vector;            -- return v1-v2;
  function "-" ( v1,v2 : Link_to_Vector ) return Link_to_Vector;
  function "*" ( v : Vector; a : number ) return Vector;   -- return v*a;
  function "*" ( v : Link_to_Vector; a : number ) return Link_to_Vector; 
  function "*" ( a : number; v : Vector ) return Vector;   -- return a*v;
  function "*" ( a : number; v : Link_to_Vector ) return Link_to_Vector;

  function "*" ( v1,v2 : Vector ) return number;           
  function "*" ( v1,v2 : Link_to_Vector ) return number; 
     -- returns the inner product of the vectors v1 and v2;
     -- requires v1 /= null /= v2, in case of link_to_vector.

  function "*" ( v1,v2 : Vector ) return Vector;
  function "*" ( v1,v2 : Link_to_Vector ) return Link_to_Vector;
     -- returns the vector v, with v(k) = v1(k)*v2(k);
     -- requires v1 /= null /= v2, in case of link_to_vector.

  function Sum ( v : Vector ) return number;
  function Sum ( v : Link_to_Vector ) return number;
     -- returns the sum of all components of v;

-- ARITHMETIC OPERATIONS AS PROCEDURES :

  procedure Add ( v1 : in out Vector; v2 : in Vector );        -- v1 := v1+v2
  procedure Add ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector ); 
  procedure Min ( v  : in out Vector);                            -- v  := -v
  procedure Min ( v  : in Link_to_Vector );
  procedure Sub ( v1 : in out Vector; v2 : in Vector );        -- v1 := v1-v2
  procedure Sub ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector );
  procedure Mul ( v  : in out Vector; a : in number );           -- v  := v*a
  procedure Mul ( v  : in Link_to_Vector; a : in number );
  procedure Mul ( v1 : in out Vector; v2 : in Vector );        -- v1 := v1*v2
  procedure Mul ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector );

-- DESTRUCTORS :

  procedure Clear ( v : in out Vector );
  procedure Clear ( v : in out Link_to_Vector );

end Generic_Vectors;
