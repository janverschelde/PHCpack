with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Solutions;

package Extrinsic_Diagonal_Homotopies is

-- DESCRIPTION :
--   This package contains tools to deal with diagonal homotopies,
--   which are used to intersect pure dimensional varieties.
--   The operations in this package primarily support the extrinsic
--   version of the diagonal homotopies, at several precision levels:
--   for standard double, double double, or quad double precision.

  function Cascade_Dimension ( n1,n2,a,b : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the cascade to intersect two sets of
  --   dimensions a and b, in spaces of respective dimensions n1 and n2.

  function Cascade_Dimension
             ( p1e,p2e : Standard_Complex_Poly_Systems.Poly_Sys;
               a,b : natural32 ) return natural32;
  function Cascade_Dimension
             ( p1e,p2e : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               a,b : natural32 ) return natural32;
  function Cascade_Dimension
             ( p1e,p2e : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               a,b : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the cascade to compute all components
  --   of the intersection of two components of dimensions a and b
  --   defined by two embedded polynomial systems p1e and p2e.
  --   For a minimal dimension, we must have that a >= b.

  function Cascade_Dimension
             ( p1e,p2e : Standard_Complex_Laur_Systems.Laur_Sys;
               a,b : natural32 ) return natural32;
  function Cascade_Dimension
             ( p1e,p2e : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               a,b : natural32 ) return natural32;
  function Cascade_Dimension
             ( p1e,p2e : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               a,b : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the cascade to compute all components
  --   of the intersection of two components of dimensions a and b
  --   defined by two embedded Laurent polynomial systems p1e and p2e.
  --   For a minimal dimension, we must have that a >= b.

  procedure Cascade1
              ( p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Cascade1
              ( p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Cascade1
              ( p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Given two polynomial systems, properly embedded to sample from
  --   positive dimensional solution components, this function returns
  --   the cascade to get all components of their intersection,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   p1e      1st system with a-dimensional solution component;
  --   p2e      2nd system with b-dimensional solution component;
  --   a        dimension of component of 1st system;
  --   b        dimension of component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start    start system to start the cascade;
  --   target   embedded cascade to compute all components of dimension < b
  --            of the intersection of two components defined by p1e and p2e.

  procedure Cascade1
              ( p1e,p2e : in Standard_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Cascade1
              ( p1e,p2e : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Cascade1
              ( p1e,p2e : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Given two Laurent systems, properly embedded to sample from
  --   positive dimensional solution components, this function returns
  --   the cascade to get all components of their intersection,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   p1e      1st system with a-dimensional solution component;
  --   p2e      2nd system with b-dimensional solution component;
  --   a        dimension of component of 1st system;
  --   b        dimension of component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start    start system to start the cascade;
  --   target   embedded cascade to compute all components of dimension < b
  --            of the intersection of two components defined by p1e and p2e.

  procedure Cascade2
              ( p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Cascade2
              ( p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Cascade2
              ( p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Given two polynomial systems, properly embedded to sample from
  --   positive dimensional solution components, this function returns
  --   the cascade to get all components of their intersection,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   p1e      1st system with a-dimensional solution component;
  --   p2e      2nd system with b-dimensional solution component;
  --   a        dimension of solution component of 1st system;
  --   b        dimension of solution component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start    start system to start the cascade;
  --   target   embedded cascade to compute all components of dimension < b
  --            of the intersection of two components defined by p1e and p2e.

  procedure Cascade2
              ( p1e,p2e : in Standard_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Cascade2
              ( p1e,p2e : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Cascade2
              ( p1e,p2e : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Given two Laurent polynomial systems, properly embedded to sample
  --   from positive dimensional solution components, this function returns
  --   the cascade to get all components of their intersection,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   p1e      1st system with a-dimensional solution component;
  --   p2e      2nd system with b-dimensional solution component;
  --   a        dimension of solution component of 1st system;
  --   b        dimension of solution component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start    start system to start the cascade;
  --   target   embedded cascade to compute all components of dimension < b
  --            of the intersection of two components defined by p1e and p2e.

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Returns homotopy to start the cascade, using Cascade1 and Cascade2,
  --   to find all components of the intersection of two witness sets, of
  --   dimensions a and b defined by p1e and p2e respectively,
  --   in standard double, double double, or quad double precision,
  --   for witness sets defined by ordinary polynomial systems.

  -- ON ENTRY :
  --   p1e      1st system with a-dimensional solution component;
  --   p2e      2nd system with b-dimensional solution component;
  --   a        dimension of solution component of 1st system;
  --   b        dimension of solution component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start    start system to start the cascade;
  --   target   embedded cascade to compute all components of dimension < b
  --            of the intersection of two components defined by p1e and p2e.

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Standard_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Returns homotopy to start the cascade, using Cascade1 and Cascade2,
  --   to find all components of the intersection of two witness sets, of
  --   dimensions a and b defined by p1e and p2e respectively,
  --   in standard double, double double, or quad double precision,
  --   for witness sets defined by Laurent polynomial systems.

  -- ON ENTRY :
  --   p1e      1st system with a-dimensional solution component;
  --   p2e      2nd system with b-dimensional solution component;
  --   a        dimension of solution component of 1st system;
  --   b        dimension of solution component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start    start system to start the cascade;
  --   target   embedded cascade to compute all components of dimension < b
  --            of the intersection of two components defined by p1e and p2e.

  function Extrinsic_Product
              ( a,b : natural32;
                s1,s2 : Standard_Complex_Solutions.Solution )
              return Standard_Complex_Solutions.Solution;
  function Extrinsic_Product
              ( a,b : natural32;
                s1,s2 : DoblDobl_Complex_Solutions.Solution )
              return DoblDobl_Complex_Solutions.Solution;
  function Extrinsic_Product
              ( a,b : natural32;
                s1,s2 : QuadDobl_Complex_Solutions.Solution )
              return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns the product of the two solutions,
  --   embedded properly to intersect sets of dimensions a and b,
  --   in standard double, double double, or quad double precision.

  function Extrinsic_Product
              ( a,b,k : natural32;
                sols1,sols2 : Standard_Complex_Solutions.Solution_List )
              return Standard_Complex_Solutions.Solution_List;
  function Extrinsic_Product
              ( a,b,k : natural32;
                sols1,sols2 : DoblDobl_Complex_Solutions.Solution_List )
              return DoblDobl_Complex_Solutions.Solution_List;
  function Extrinsic_Product
              ( a,b,k : natural32;
                sols1,sols2 : QuadDobl_Complex_Solutions.Solution_List )
              return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the product of the two solution lists, embedded properly
  --   for use in the homotopy to start the cascade in the diagonal homotopy,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   a        dimension of the first solution set;
  --   b        dimension of the second solution set;
  --   k        ambient dimension;
  --   sols1    witness points in the 1st witness set, without the embedding;
  --   sols2    witness points in the 2nd witness set, without the embedding.

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                sols1,sols2 : in Standard_Complex_Solutions.Solution_List;
                start,target : out Standard_Complex_Poly_Systems.Poly_Sys;
                esols : out Standard_Complex_Solutions.Solution_List );
  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                sols1,sols2 : in DoblDobl_Complex_Solutions.Solution_List;
                start,target : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                esols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                sols1,sols2 : in QuadDobl_Complex_Solutions.Solution_List;
                start,target : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                esols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Uses Cascade1 and Cascade2 to create the homotopy to start the cascade
  --   to find all components of the intersection of two witness sets, of
  --   dimensions a and b defined by (p1e,sols1) and (p2e,sols2) respectively,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   p1e      1st system with a-dimensional solution component;
  --   p2e      2nd system with b-dimensional solution component;
  --   sols1    witness points in 1st witness set, without the embedding;
  --   sols2    witness points in 2nd witness set, without the embedding;
  --   a        dimension of solution component of 1st system;
  --   b        dimension of solution component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start    start system to start the cascade;
  --   target   embedded cascade to compute all components of dimension < b
  --            of the intersection of two components defined by p1e and p2e;
  --   esols    embedded solutions of the start system in the cascade.

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Standard_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                sols1,sols2 : in Standard_Complex_Solutions.Solution_List;
                start,target : out Standard_Complex_Laur_Systems.Laur_Sys;
                esols : out Standard_Complex_Solutions.Solution_List );
  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                sols1,sols2 : in DoblDobl_Complex_Solutions.Solution_List;
                start,target : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                esols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                sols1,sols2 : in QuadDobl_Complex_Solutions.Solution_List;
                start,target : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                esols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Uses Cascade1 and Cascade2 to create the homotopy to start the cascade
  --   to find all components of the intersection of two witness sets, of
  --   dimensions a and b defined by (p1e,sols1) and (p2e,sols2) respectively,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   p1e      1st system with a-dimensional solution component;
  --   p2e      2nd system with b-dimensional solution component;
  --   sols1    witness points in 1st witness set, without the embedding;
  --   sols2    witness points in 2nd witness set, without the embedding;
  --   a        dimension of solution component of 1st system;
  --   b        dimension of solution component of 2nd system.

  -- REQUIRED :
  --   a >= b and a+b < k, k = ambient dimension;
  --   start'range = target'range = 1..Cascade_Dimension(p1e,p2e,a,b).

  -- ON RETURN :
  --   start    start system to start the cascade;
  --   target   embedded cascade to compute all components of dimension < b
  --            of the intersection of two components defined by p1e and p2e;
  --   esols    embedded solutions of the start system in the cascade.

end Extrinsic_Diagonal_Homotopies;
