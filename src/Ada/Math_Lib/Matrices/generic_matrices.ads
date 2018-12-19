with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring,Generic_Vectors;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);

package Generic_Matrices is

-- DESCRIPTION :
--   An abstraction for matrices with coefficients over any ring.

  use Ring;  use Vectors;

  type Matrix is array ( integer32 range <>, integer32 range <> ) of number;
  type Link_to_Matrix is access Matrix;

-- COMPARISON AND COPYING :

  function Equal ( a,b : Matrix ) return boolean;

  -- DESCRIPTION :
  --   Comparing by using the comparison operations of the Ring.

  procedure Copy ( a : in Matrix; b : in out Matrix );

  -- DESCRIPTION :
  --   Makes a deep copy of the matrix a to the matrix b.

  -- REQUIRED : a'range(1) = b'range(1) and a'range(2) = b'range(2).

-- MATRIX-MATRIX OPERATIONS :

  function Transpose ( a : Matrix ) return Matrix;         -- return a'

  function "+" ( a,b : Matrix ) return Matrix;             -- return a+b
  function "+" ( a : Matrix ) return Matrix;               -- copies a
  function "-" ( a,b : Matrix ) return Matrix;             -- return a-b
  function "-" ( a : Matrix ) return Matrix;               -- return -a
  function "*" ( a,b : Matrix ) return Matrix;             -- return a*b

  procedure Add ( a : in out Matrix; b : in Matrix );      -- a := a+b
  procedure Sub ( a : in out Matrix; b : in Matrix );      -- a := a-b

  procedure Mul1 ( a : in out Matrix; b : in Matrix );     -- a := a*b    
  procedure Mul2 ( a : in Matrix; b : in out Matrix );     -- b := a*b

-- MATRIX-VECTOR OPERATIONS :

  function "*" ( a : Matrix; v : Vector ) return Vector;   -- return a*v
  function "*" ( v : Vector; a : Matrix ) return Vector;   -- return v*a

  procedure Mul ( a : in Matrix; v : in out Vector );      -- v := a*v
  procedure Mul ( v : in out Vector; a : in Matrix );      -- v := a*v

-- SCALING A MATRIX :

  function "*" ( a : Matrix; x : number ) return Matrix;   -- return x*a
  function "*" ( x : number; a : Matrix ) return Matrix;   -- return a*x

  procedure Mul ( a : in out Matrix; x : in number );      -- a := a*x

-- DESTRUCTORS :

  procedure Clear ( a : in out Matrix );
  procedure Clear ( a : in out Link_to_Matrix );

  -- DESCRIPTION :
  --   Deallocation of memory.

end Generic_Matrices;
