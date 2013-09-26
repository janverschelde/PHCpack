with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Multitasking_Volume_Computation is

-- DESCRIPTION :
--   Offers routines to compute the mixed volume with multiple tasks.
--   There are two types of routines:
--   (1) silent or reporting the progress of the computations;
--   (2) static or dynamic load balancing:
--   with static, cell k is process by task 1 + (k mod #tasks),
--   with dynamic, tasks see the collection of cells as a queue.

  procedure Extract_Support_Vectors
              ( mic : in Mixed_Cell;
                A : out Standard_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given the points in the mixed cell, fills the rows of A
  --   with the support vectors.

  -- REQUIRED :
  --   A is a square matrix of dimension equal to the dimension
  --   of the ambient space, that is: before lifting.

  -- ON INPUT :
  --   mic      a mixed cell.
 
  -- ON RETURN :
  --   A        the rows of A contain the vectors computed as the
  --            differences between the points in the mixed cell
  --            with the first point in each support of the cell,
  --            |det(A)| equals the mixed volume of the cell.

  procedure Reporting_Static_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32 );
  procedure Reporting_Dynamic_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32 );
  procedure Reporting_Static_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                vol : out Standard_Integer_Vectors.Vector;
                mixvol : out natural32 );
  procedure Reporting_Dynamic_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                vol : out Standard_Integer_Vectors.Vector;
                mixvol : out natural32 );

  -- DESCRIPTION :
  --   Given a mixed cell configuration, computes the mixed volume.
  --   Intermediate output is written to screen.
  --   Static and dynamic load balancing is indidated respectively
  --   by the _Static_ and _Dynamic_ parts of the procedure name.

  -- REQUIRED : vol'range = 1..nt.

  -- ON ENTRY :
  --   nt       number of tasks in the computation;
  --   dim      dimension of the ambient space (before lifting);
  --   mcc      a regular mixed cell configuration.

  -- ON RETURN :
  --   vol      (optional) vol(t) contains the contribution to the
  --            mixed volume computed by thread t;
  --            mixvol equals the sum of vol(t) for t in 1..nt.
  --   mixvol   the mixed volume of the cells in mcc.

  procedure Silent_Static_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32 );
  procedure Silent_Dynamic_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32 );
  procedure Silent_Static_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                vol : out Standard_Integer_Vectors.Vector; 
                mixvol : out natural32 );
  procedure Silent_Dynamic_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                vol : out Standard_Integer_Vectors.Vector; 
                mixvol : out natural32 );

  -- DESCRIPTION :
  --   Given a mixed cell configuration, computes the mixed volume.
  --   No intermediate output is written to screen.
  --   Static and dynamic load balancing is indidated respectively
  --   by the _Static_ and _Dynamic_ parts of the procedure name.

  -- REQUIRED : vol'range = 1..nt.

  -- ON ENTRY :
  --   nt       number of tasks in the computation;
  --   dim      dimension of the ambient space (before lifting);
  --   mcc      a regular mixed cell configuration.

  -- ON RETURN :
  --   vol      (optional) vol(t) contains the contribution to the
  --            mixed volume computed by thread t;
  --            mixvol equals the sum of vol(t) for t in 1..nt.
  --   mixvol   the mixed volume of the cells in mcc.

end Multitasking_Volume_Computation;
