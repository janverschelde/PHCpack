with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Extrinsic_Diagonal_Continuation is

-- DESCRIPTION :
--   A diagonal homotopy allows to compute witness sets for all components
--   of the intersection of two positive dimensional solution sets.
--   The extrinsic version applies the equations for the witness sets.

-- IMPORTANT REQUIREMENT :
--   The polynomial systems and system functions assume the right
--   number of equations, i.e.: complete intersections.
--   Randomize whenever necessary before applying diagonal homotopies.

  function Minimal_Intersection_Dimension
             ( n,a,b : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the minimal dimension of the intersection of 
  --   an a-dimensional with a b-dimensional component in n-space.

  procedure Start_Extrinsic_Cascade
              ( file : in file_type; report : in boolean;
                p,q : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Does the path following to start the extrinsic cascade.

  -- ON ENTRY :
  --   file     for output and intermediate diagnostics;
  --   report   if reporting version of path trackers is needed;
  --   p        target polynomial system;
  --   q        system to start the cascade;
  --   sols     solutions of the start system q.

  -- ON RETURN :
  --   sols     solutions at the end of converging paths.

  procedure Down_Extrinsic_Cascade
              ( file : in file_type; name : in string; report : in boolean;
                q : in Poly_Sys; sols : in out Solution_List;
                p1,p2 : in Poly_Sys; sli : in VecVec;
                tol : in double_float; n,k : in integer32 );

  -- DESCRIPTION :
  --   Does one path following stage down in the extrinsic cascade.

  -- ON ENTRY :
  --   file     for diagnostics and intermediate results;
  --   name     to write the (super) witness sets to;
  --   report   flag to ask for intermediate output of path trackers;
  --   q        embedded system, serves as start system;
  --   sols     start solutions, all with nonzero slack variable;
  --   p1       original system of the 1st witness set;
  --   p2       original system of the 2nd witness set;
  --   sli      set of k hyperplanes for the top embedding;
  --   tol      tolerance used by the filter;
  --   n        dimension of the ambient space;
  --   k        top dimension.

  -- ON RETURN :
  --   sols     solutions at the end of the paths.

  procedure Extrinsic_Diagonal_Homotopy
              ( file : in file_type; name : in string; report : in boolean;
                p1e,p2e : in Poly_Sys; a,b : in natural32;
                sols1,sols2 : in Solution_List );

  -- DESCRIPTION :
  --   Runs the diagonal homotopy algorithm in extrinsic coordinates
  --   to intersect a-dimensional witness set (p1e,sols1)
  --   with the b-dimensional witness set (p2e,sols2).

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   name     name of the output file, will be used as prefix to write
  --            the (super) witness sets;
  --   report   flag to ask for intermediate output during path following;
  --   p1e      1st embedded polynomial system;
  --   p2e      2nd embedded polynomial system;
  --   a        dimension of the 1st witness set;
  --   b        dimension of the 2nd witness set;
  --   sols1    witness points on the 1st component, without embedding;
  --   sols2    witness points on the 2nd component, without embedding.

  -- REQUIRED : a >= b.

end Extrinsic_Diagonal_Continuation;
