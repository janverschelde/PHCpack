with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_natural_Numbers;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Pieri_Root_Counts;                  use Pieri_Root_Counts;

package Pieri_Deformations is

-- DESCRIPTION :
--   Deformations that start at pairs of leaves satisfying Pieri's conditions.

  procedure Deform_Pair
              ( file : in file_type; pnd : in Paired_Nodes; id : in natural32;
                f1,f2 : in Standard_Complex_Matrices.Matrix;
                l1,l2 : in VecMat; ln : in Matrix; report,outlog : in boolean;
                sol : in out Matrix );

  -- DESCRIPTION :
  --   Does one step in the Pieri deformation at the current pair and
  --   moves down to the node below in the chains until the solution sol
  --   meets all planes in l1,l2 and ln nontrivially.

  -- ON ENTRY :
  --   file     for writing intermediate results;
  --   pnd      pair of nodes in a chain that starts at a pair of leaves
  --            for which Pieri's condition is satisfied;
  --   id       identity number of the pair in the list of paired nodes;
  --   f1       random upper triangular matrix with 1's on its anti-diagonal;
  --   f2       random lower triangular matrix with 1's on its diagonal;
  --   l1       first sequence of input planes;
  --   l2       second sequence of input planes;
  --   ln       last input plane;
  --   report   indicates whether intermediate output during path tracking;
  --   sol      solution at the pair of nodes above the current pair.

  -- ON ENTRY :
  --   sol      updated solution plane.

  procedure Deform_Pairs
              ( file : in file_type; n,d : in natural32;
                lp : in List_of_Paired_Nodes; l1,l2 : in VecMat;
                ln : in Matrix; report,outlog : in boolean; sols : out VecMat );

  -- DESCRIPTION :
  --   Performs the deformation of pairs for every pair in the list.

  -- ON ENTRY :
  --   file     to write intermediate results;
  --   n        dimension of the space the planes all live in;
  --   d        dimension of the output planes;
  --   lp       list of paired nodes, satisfying Pieri's condition;
  --   l1       first sequence of input planes,
  --            l1(0) is spanned by first standard basis vectors;
  --   l2       second sequence of input planes,
  --            l2(0) is spanned by last standard basis vectors;
  --   ln       last input plane;
  --   report   indicates whether intermediate output during path tracking;
  --   outlog   if switched on, writes moving cycles and polynomial systems.

  -- ON RETURN :
  --   sols     sequence of solution planes.

end Pieri_Deformations;
