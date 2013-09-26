with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Lifted_Configurations is

-- DESCRIPTION :
--   This package provides facilities for examining the lifting values
--   assigned to points of a configuration to construct a subdivision.

  type Incidence_Matrix is
    array ( integer32 range <>, integer32 range <> ) of boolean;

  procedure Collect_Points_and_Lifting
               ( n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                 mixsub : in out Mixed_Subdivision;
                 pts : out Array_of_Lists; nbp : out natural32;
                 lifting : out Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Collects the points that are used as vertices in the subdivision.

  -- ON ENTRY :
  --   n         dimension of the points in the unlifted configuration;
  --   mix       mix(i) equals the number of occurences of the ith polytope;
  --   mixsub    collection of lifted mixed cells.

  -- ON RETURN :
  --   mixsub    lifting is replaced by index of the point;
  --   pts       supports, last coordinate is the index of the point;
  --   nbp       number of points in the supports;
  --   lifting   vector of lifting values, range of lif is 1..nbp.

  procedure Get_Point ( pts : in Array_of_Lists; k : in integer32;
                        lpk : out Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Searches the point with index k in the lists pts and assigns
  --   the pointer to lpk if found, otherwise lpk is the null pointer.

  function Create ( np : integer32; mixsub : Mixed_Subdivision )
                  return Incidence_Matrix;

  -- DESCRIPTION :
  --   The (i,j)th entry of the matrix on return answers the question
  --   whether the ith point belongs to the jth cell.

  -- REQUIRED :
  --   The last coordinate of each point in the cells must correspond
  --   to its index.

  procedure Sort ( im : in Incidence_Matrix;
                   order : out Standard_Integer_Vectors.Vector;
                   fail : out boolean );

  -- DESCRIPTION :
  --   Sorts the points so that the order is placeable.
  --   If fail is false, then no such order could be found.
  --   This procedure assumes that the cells were placed as on a stack,
  --   i.e.: that the last cell was always added in front of the list.
  --   The order on return also start with the last point added.

  procedure Assign_Lifting
              ( n : in integer32; mixsub : in out Mixed_Subdivision;
                lifting : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Assigns the lifting to points in the cells of mixsub.

  -- REQUIRED :
  --   The last coordinate of each point in the cells must correspond
  --   to its index.

end Lifted_Configurations;
