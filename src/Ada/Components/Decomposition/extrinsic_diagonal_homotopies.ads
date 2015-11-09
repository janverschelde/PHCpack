with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Permutations;                      use Permutations;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Extrinsic_Diagonal_Homotopies is

-- DESCRIPTION :
--   This package contains tools to deal with diagonal homotopies,
--   which are used to intersect pure dimensional varieties.
--   The operations in this package primarily support the extrinsic
--   version of the diagonal homotopies.

  function Cascade_Dimension ( n1,n2,a,b : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the cascade to intersect two sets of
  --   dimensions a and b, in spaces of respective dimensions n1 and n2.

  function Cascade_Dimension
             ( p1e,p2e : Poly_Sys; a,b : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the cascade to compute all components
  --   of the intersection of two components of dimensions a and b
  --   defined by two embedded polynomial systems p1e and p2e.
  --   For a minimal dimension, we must have that a >= b.

  procedure Cascade1 ( p1e,p2e : in Poly_Sys; a,b : in natural32;
                       start,target : out Poly_Sys );

  -- DESCRIPTION :
  --   Given two polynomial systems, properly embedded to sample from
  --   positive dimensional solution components, this function returns
  --   the cascade to get all components of their intersection.

  -- ON ENTRY :
  --   p1e       1st system with a-dimensional solution component;
  --   p2e       2nd system with b-dimensional solution component;
  --   a         dimension of component of 1st system;
  --   b         dimension of component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start     start system to start the cascade;
  --   target    embedded cascade to compute all components of dimension < b
  --             of the intersection of two components defined by p1e and p2e.

  procedure Cascade2 ( p1e,p2e : in Poly_Sys; a,b : in natural32;
                       start,target : out Poly_Sys );

  -- DESCRIPTION :
  --   Given two polynomial systems, properly embedded to sample from
  --   positive dimensional solution components, this function returns
  --   the cascade to get all components of their intersection.

  -- ON ENTRY :
  --   p1e       1st system with a-dimensional solution component;
  --   p2e       2nd system with b-dimensional solution component;
  --   a         dimension of solution component of 1st system;
  --   b         dimension of solution component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start     start system to start the cascade;
  --   target    embedded cascade to compute all components of dimension < b
  --             of the intersection of two components defined by p1e and p2e.

  procedure Extrinsic_Cascade_Homotopy
               ( p1e,p2e : in Poly_Sys; a,b : in natural32;
                 start,target : out Poly_Sys );

  -- DESCRIPTION :
  --   Returns homotopy to start the cascade, using Cascade1 and Cascade2,
  --   to find all components of the intersection of two witness sets, of
  --   dimensions a and b defined by p1e and p2e respectively.

  -- ON ENTRY :
  --   p1e       1st system with a-dimensional solution component;
  --   p2e       2nd system with b-dimensional solution component;
  --   a         dimension of solution component of 1st system;
  --   b         dimension of solution component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start     start system to start the cascade;
  --   target    embedded cascade to compute all components of dimension < b
  --             of the intersection of two components defined by p1e and p2e.

  function Extrinsic_Product
               ( a,b : natural32; s1,s2 : Solution ) return Solution;

  -- DESCRIPTION :
  --   Returns the product of the two solutions, embedded properly to
  --   intersect sets of dimensions a and b.

  function Extrinsic_Product
               ( a,b,k : natural32; sols1,sols2 : Solution_List )
               return Solution_List;

  -- DESCRIPTION :
  --   Returns the product of the two solution lists, embedded properly
  --   for use in the homotopy to start the cascade in the diagonal homotopy.

  -- ON ENTRY :
  --   a         dimension of the first solution set;
  --   b         dimension of the second solution set;
  --   k         ambient dimension;
  --   sols1     witness points in the 1st witness set, without the embedding;
  --   sols2     witness points in the 2nd witness set, without the embedding.

  procedure Extrinsic_Cascade_Homotopy
               ( p1e,p2e : in Poly_Sys; a,b : in natural32;
                 sols1,sols2 : in Solution_List;
                 start,target : out Poly_Sys; esols : out Solution_List );

  -- DESCRIPTION :
  --   Uses Cascade1 and Cascade2 to create the homotopy to start the cascade
  --   to find all components of the intersection of two witness sets, of
  --   dimensions a and b defined by (p1e,sols1) and (p2e,sols2) respectively.

  -- ON ENTRY :
  --   p1e       1st system with a-dimensional solution component;
  --   p2e       2nd system with b-dimensional solution component;
  --   sols1     witness points in 1st witness set, without the embedding;
  --   sols2     witness points in 2nd witness set, without the embedding;
  --   a         dimension of solution component of 1st system;
  --   b         dimension of solution component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start     start system to start the cascade;
  --   target    embedded cascade to compute all components of dimension < b
  --             of the intersection of two components defined by p1e and p2e;
  --   esols     embedded solutions of the start system in the cascade.

end Extrinsic_Diagonal_Homotopies;
