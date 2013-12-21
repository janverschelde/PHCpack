with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Multihomogeneous_Solutions is

-- DESCRIPTION :
--   A multi-homogeneous k-dimensional component is represented by an array
--   of solution lists and a corresponding vector of special directions.
--   This package offers tools to manage such a representation.

-- DATA STRUCTURES :

  type Multihomogeneous_Solution ( n,nb : natural ) is record
    k : natural;                         -- co-dimension of the solution set
    b : Vector(1..n);                    -- offset vector for affine planes
    w : Array_of_VecVecs(1..nb);         -- directions spanning planes
    s : Array_of_Solution_Lists(1..nb);  -- linear combination coefficients
  end record;

  type Link_to_Multihomogeneous_Solution is access Multihomogeneous_Solution;

-- CREATORS :

  function Create ( n : natural; b : Vector; w : Array_of_VecVecs;
                    s : Array_of_Solution_Lists )
                  return Multihomogeneous_Solution;
  function Create ( n : natural; b : Vector; w : Array_of_VecVecs;
                    s : Array_of_Solution_Lists )
                  return Link_to_Multihomogeneous_Solution;

  -- DESCRIPTION :
  --   Removes empty entries in w and s, orthogonalizes the directions in w,
  --   and rewrites the solution lists in terms of the new orthogonal bases.

-- DESTRUCTORS :

  procedure Clear ( mhs : in out Multihomogeneous_Solution );
  procedure Clear ( mhs : in out Link_to_Multihomogeneous_Solution );

  -- DESCRIPTION :
  --   Releases the allocated memory.

end Multihomogeneous_Solutions;