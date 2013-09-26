with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer64_Vectors;         use Standard_Integer64_Vectors;
with Standard_Integer64_VecVecs;         use Standard_Integer64_VecVecs;
with Standard_Integer64_Matrices;        use Standard_Integer64_Matrices;
with Standard_Complex_Vectors;

package Standard_Integer64_Transformations is

-- DESCRIPTION :
--   This package implements monomial transformations.

  type Transfo is private;

-- CONSTRUCTORS :

  function Id ( n : natural32 ) return Transfo;

  -- DESCRIPTION :
  --   Returns the identical transformation.

  function Create ( v : Vector; i : integer32 ) return Transfo;

  -- DESCRIPTION :
  --   returns a transformation T such that after w = Apply(T,v),
  --   there exists a j /= i: w(i) = 1 and w(j) = 0;
  --   if v(i) = 0, then the identity transformation will be returned;
  --   T is an isomorphism.

  function Create ( v : VecVec ) return Transfo;

  -- DESCRIPTION :
  --   returns a transformation T, where v(i) is the image of
  --   the i-th basis vector, for i in da'range.

  function Create ( m : Matrix ) return Transfo;

  -- DESCRIPTION :
  --   Returns the transformation defined by the matrix m.
  --   Apply(Create(m),v) = m*v.

  function Rotate ( v : Vector; i : integer32 ) return Transfo;
  function Rotate ( v : Link_to_Vector; i : integer32 ) return Transfo;

  -- DESCRIPTION :
  --   Returns a transformation T which will reduce v into e_i;
  --   if v(i) /= 0; T is an isomorphism.

  function Build_Transfo ( v : Vector; i : integer32 ) return Transfo;
  function Build_Transfo ( v : Link_to_Vector; i : integer32 ) return Transfo;

  -- DESCRIPTION :
  --   Returns a transformation T so that for all vectors x
  --   < x , v > = pv, after application of T on all x, the
  --   following holds: x(i) = pv / d, where d is the gcd of
  --   the components of v; T is an isomorphism.

-- SELECTORS :

  function Dimension ( t : Transfo ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the transformation.

  function Sign ( t : Transfo ) return integer32;

  -- DESCRIPTION :
  --   Returns the determinant of the matrix representation of t,
  --   which is either +1 or -1.

-- OPERATIONS :

  function Transpose ( t : Transfo ) return Transfo;

  -- DESCRIPTION :
  --   Returns the transposed transformation.

  function Invert ( t : Transfo ) return Transfo;

  -- DESCRIPTION :
  --   Computes the inverse transformation of t;
  --   after t1 := Invert(t), t1*t = t*t1 = Id.

  -- REQUIRED : t is an isomorphism.

  function "*" ( t1,t2 : Transfo ) return Transfo;

  -- DESCRIPTION :
  --   Returns t1 after t2

  procedure Mult1 ( t1 : in out Transfo; t2 : in Transfo );

  -- DESCRIPTION :
  --   t1 := t1 * t2;  but with an efficient memory management

  procedure Mult2 ( t1 : in Transfo; t2 : in out Transfo );

  -- DESCRIPTION :
  --   t2 := t1 * t2;  but with an efficient memory management

  function "*"( t : Transfo; v : Vector ) return Vector;  -- return t*v;
  function "*"( t : Transfo; v : Link_to_Vector ) return Link_to_Vector;
  procedure Apply ( t : in Transfo; v : in out Vector );  -- v := t*v;
  procedure Apply ( t : in Transfo; v : in Link_to_Vector );

  function "*" ( t : Transfo; v : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector;
  procedure Apply ( t : in Transfo;
                    v : in out Standard_Complex_Vectors.Vector ); -- v := t*v;

-- DESTRUCTOR :

  procedure Clear ( t : in out Transfo );

  -- DESCRIPTION :
  --   Frees all allocated memory space

private

  type Transfo_tp;
  type Transfo is access Transfo_tp;
  type Transfo_tp is new matrix;

end Standard_Integer64_Transformations;
