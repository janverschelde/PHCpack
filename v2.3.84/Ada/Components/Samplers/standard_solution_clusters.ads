with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;         
with Standard_Complex_Solutions;

package Standard_Solution_Clusters is

-- DESCRIPTION :
--   A solution cluster is a list of solutions close to each other.

  procedure Separate_Clusters
               ( file : in file_type; tol : in double_float;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 nbclus : out integer32 );

  -- DESCRIPTION :
  --   Separates solutions in different clusters, depending on tolerance.

  -- REQUIRED : not Is_Null(sols).

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   tol       tolerance to decide whether two solutions are clustered;
  --   sols      given list of solutions.

  -- ON RETURN :
  --   sols      multiplicity field is cluster number;
  --   nbclus    number of clusters in the list of solutions.

  function Select_Cluster
               ( sols : Standard_Complex_Solutions.Solution_List;
                 k : integer32 )
               return Standard_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Given the list with separated clusters, the solutions in the k-th
  --   cluster are returned.

  function Center_of_Cluster
               ( sols : Standard_Complex_Solutions.Solution_List )
               return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the average of the solution vectors in the cluster.

  -- REQUIRED : not Is_Null(sols).

  procedure Distances_to_Center
               ( sols : in Standard_Complex_Solutions.Solution_List;
                 center : in Standard_Complex_Vectors.Vector;
                 max,min,avg : out double_float );

  -- DESCRIPTION :
  --   Returns the maximal, minimal and average distance of the solution
  --   vectors in the solution cluster to the center vector.

end Standard_Solution_Clusters;
