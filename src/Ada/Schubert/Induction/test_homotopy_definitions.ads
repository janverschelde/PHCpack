with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Natural_Matrices;          use Standard_Natural_Matrices;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Checker_Posets;                     use Checker_Posets;
with Intersection_Posets;                use Intersection_Posets;

package Test_Homotopy_Definitions is

-- DESCRIPTION :
--   Interactively tests the setup and the running of the
--   Littlewood-Richardson homotopies.

  procedure Write_Short ( x : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the matrix in short format, without the imaginary part,
  --   and as an integer.

  procedure Generalizing_Moving_Flags ( n : in integer32 );

  -- DESCRIPTION :
  --   Shows all boards, flags and transformations in a generalizing
  --   moving flag homotopy towards a general triangular matrix.

  procedure Test_Coordinates ( n,k : in integer32 );
 
  -- DESCRIPTION :
  --   Reads in a checker board and sets up the coordinates
  --   associated with the positions of white and black checkers.

  procedure Specializing_Homotopies ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Walks through the poset starting at the root, printing information
  --   about the deformations involving a k-plane in n-space.

  procedure Descend_Poset
              ( ps : in Poset; n,k,leaf : in integer32;
                nd : in Link_to_Node );

  -- DESCRIPTION :
  --   Interactive descent in the poset from a leaf to the root.

  -- ON ENTRY :
  --   ps       poset created to deform k-planes in n-space;
  --   n        dimension of the ambient space;
  --   k        dimension of the subspace;
  --   leaf     label of a leaf to start at;
  --   nd       node at the leaf, should not be null.

  procedure Start_Descent_in_Poset ( n,k : in integer32; ps : in Poset );

  -- DESCRIPTION :
  --   Prompts the user to give a leaf in the poset to start the descent.

  procedure Generalizing_Homotopies ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Creates a poset to deform k-planes in n-space.
  --   After prompting the user for the number of a leaf, when moving to 
  --   the root, information is printed about the homotopies.

  procedure Descend_Intersection_Poset
              ( n,k : in integer32; ips : in Intersection_Poset;
                leaf : in Link_to_Poset_Node );

  -- DESCRIPTION :
  --   Interactive descent in the intersection poset to the root,
  --   starting at the given leaf.

  procedure Intersecting_Homotopies ( n,k,m : in integer32 );

  -- DESCRIPTION :
  --   Prompts for two intersection conditions at the start
  --   of an intersection poset for m-1 levels.

  procedure Localization_along_Paths ( n,k : integer32; ps : in Poset );

  -- DESCRIPTION :
  --   Writes for every node along the paths in the poset
  --   the localization pattern.

  procedure Run_along_Paths_in_Poset ( n,k : integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a configuration of white checkers,
  --   creates a poset and writes all paths through the poset.
  --   For each node the localization pattern is written.

  procedure real_small_put ( A : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the real part of the matrix A with only one decimal place
  --   of precision, enough to see the structure of the matrix.

  procedure real_small_put ( A : in DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the real part of the matrix A with only one decimal place
  --   of precision, enough to see the structure of the matrix.

  procedure real_small_put ( A : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the real part of the matrix A with only one decimal place
  --   of precision, enough to see the structure of the matrix.

  procedure Standard_Test_Flag_Transformation
              ( f1,f2,g1,g2,A,T1,T2 : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Tests whether A*f1 = g1*T1 and A*f2 = g2*T2,
  --   using standard double precision arithmetic.

  procedure DoblDobl_Test_Flag_Transformation
              ( f1,f2,g1,g2,A,T1,T2 : in DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Tests whether A*f1 = g1*T1 and A*f2 = g2*T2,
  --   using double double precision arithmetic.

  procedure QuadDobl_Test_Flag_Transformation
              ( f1,f2,g1,g2,A,T1,T2 : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Tests whether A*f1 = g1*T1 and A*f2 = g2*T2,
  --   using quad double precision arithmetic.

  procedure Standard_Test_Random_Flags2Flags ( n : in integer32 );

  -- DESCRIPTION :
  --   Generates two pairs of random flags (F1, F2) and (G1, G2) in n-space
  --   and computes the matrix A, upper triangular matrices T1 and T2
  --   so that A*F1 = G1*T1 and A*F2 = G2*T2.
  --   The test is done in standard double precision.

  procedure DoblDobl_Test_Random_Flags2Flags ( n : in integer32 );

  -- DESCRIPTION :
  --   Generates two pairs of random flags (F1, F2) and (G1, G2) in n-space
  --   and computes the matrix A, upper triangular matrices T1 and T2
  --   so that A*F1 = G1*T1 and A*F2 = G2*T2.
  --   The test is done in double double precision.

  procedure QuadDobl_Test_Random_Flags2Flags ( n : in integer32 );

  -- DESCRIPTION :
  --   Generates two pairs of random flags (F1, F2) and (G1, G2) in n-space
  --   and computes the matrix A, upper triangular matrices T1 and T2
  --   so that A*F1 = G1*T1 and A*F2 = G2*T2.
  --   The test is done in quad double precision.

  procedure Standard_Test_Specific_Flags2Flags ( n : in integer32 );

  -- DESCRIPTION :
  --   For two pairs of flags (M, I) and (I, F) in n-space
  --   where M is the moved flags, I the identity matrix,
  --   and F some random n-dimensional matrix,
  --   computes the matrix A, upper triangular matrices T1 and T2
  --   so that A*F1 = G1*T1 and A*F2 = G2*T2.
  --   The test is executed in standard double precision.

  procedure DoblDobl_Test_Specific_Flags2Flags ( n : in integer32 );

  -- DESCRIPTION :
  --   For two pairs of flags (M, I) and (I, F) in n-space
  --   where M is the moved flags, I the identity matrix,
  --   and F some random n-dimensional matrix,
  --   computes the matrix A, upper triangular matrices T1 and T2
  --   so that A*F1 = G1*T1 and A*F2 = G2*T2.
  --   The test is executed in double double precision.

  procedure QuadDobl_Test_Specific_Flags2Flags ( n : in integer32 );

  -- DESCRIPTION :
  --   For two pairs of flags (M, I) and (I, F) in n-space
  --   where M is the moved flags, I the identity matrix,
  --   and F some random n-dimensional matrix,
  --   computes the matrix A, upper triangular matrices T1 and T2
  --   so that A*F1 = G1*T1 and A*F2 = G2*T2.
  --   The test is executed in double double precision.

  procedure Test_Flag_Transformations ( n : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for the level of precision and then generates
  --   random flags in n-space and test the transformations.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.
 
end Test_Homotopy_Definitions;
