with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Sample_Point_Lists;                use Sample_Point_Lists; 
with Sample_Point_Grids;                use Sample_Point_Grids;
with Monodromy_Group_Actions;           use Monodromy_Group_Actions;

package Monodromy_Actions_Breakup is

-- DESCRIPTION :
--   This routine offers some basic routines to break up an equidimensional
--   solution set into irreducible components.

  function Compare_Labels ( file : in file_type; tol : in double_float;
                            sps1,sps2 : in Standard_Sample_List )
                          return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   For two lists of the same samples, the labels of the solutions
  --   are compared with each other.  The i-th component of the vector
  --   on return indicates the position of the i-th sample of sps1 in
  --   the list sps2.  When the i-th sample of sps1 does not occur,
  --   then an error message is printed, and the bug is "fixed" by
  --   assigning i in the i-th entry on return.

  procedure Data_Management
                ( file : in file_type;
                  rel : in Standard_Natural_Vectors.Vector;
                  ic : in out Irreducible_Components;
                  n1,n2,cnt : in out natural32 );

  -- DESCRIPTION :
  --   Manages the result of one monodromy iteration.

  -- ON ENTRY :
  --   file       for intermediate output and diagnostics;
  --   rel        rel(i) is related to point i;
  --   ic         current grouping of generic points;
  --   n1         previous number of representative sets in ic;
  --   n2         current number of sets in ic;
  --   cnt        number of consecutive stabilizations.

  -- ON RETURN :
  --   ic         new grouping of the generic points;
  --   n1         equals the value of the incoming n2;
  --   n2         cardinality of the newly computed ic;
  --   cnt        updated counter.

  procedure Generic_Points
                    ( sps : in Standard_Sample_List;
                      labels : in Standard_Natural_Vectors.Vector;
                      gp,gp_last : in out Standard_Sample_List );

  -- DESCRIPTION :
  --   Returns the list of generic points, selected from the list sps,
  --   that lie on the component corresponding to the set of labels.

  procedure Breakup ( file : in file_type; sps : in Standard_Sample_List;
                      threshold : in natural32; tol : in double_float;
                      ic : out Irreducible_Components; nbit : out natural32 );

  -- DESCRIPTION :
  --   Generates new samples and groups points connected by paths.

  -- REQUIRED : Sampling Machine is initialized and tuned properly.

  -- ON ENTRY :
  --   file           for diagonostics and intermediate output;
  --   sps            list of samples, must be regular solutions;
  --   threshold      stabilizing threshold on iterations;
  --   tol            tolerance to decide equality of solution vectors.

  -- ON RETURN :
  --   ic             contains sets of labels of generic points that
  --                  are certain to belong to the same component;
  --   nit            total number of iterations used.

  procedure Breakup ( file : in file_type; sps : in Standard_Sample_List;
                      threshold : in natural32; tol : in double_float;
                      ic : out Irreducible_Components; nbit : out natural32;
                      grid,grid_last : in out Standard_Sample_Grid );

  -- DESCRIPTION :
  --   This procedure also returns all lists of samples that are generated.
  --   The samples here are generated on parallel hyperplanes.

end Monodromy_Actions_Breakup;
