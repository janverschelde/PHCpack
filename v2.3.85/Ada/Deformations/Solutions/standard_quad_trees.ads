with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Point_Lists;               use Standard_Point_Lists;

package Standard_Quad_Trees is

-- DESCRIPTION :
--   This package provides facilities to build quad trees for solutions
--   of polynomials, represented with standard floating-point numbers.

-- DATA STRUCTURES :

  type Quad_Node;
  type Link_to_Quad_Node is access Quad_Node;

  type Quad_Node ( leaf : boolean ) is record
    depth : natural32;        -- depth of the node in the quad tree
    size : natural32;         -- total number of points
    case leaf is
      when true =>
        pts : Point_List;                  -- points at leaf
      when false =>
        ne,nw,sw,se : Link_to_Quad_Node;   -- four quadrants
    end case;
  end record;

-- PARTITIONING INTO QUADRANTS :

  procedure Partition
               ( pl : in Point_List; cx,cy : in double_float;
                 ne_cnt,nw_cnt,sw_cnt,se_cnt : out natural32;
                 ne_pts,nw_pts,sw_pts,se_pts : out Point_List );

  -- DESCRIPTION :
  --   Partitions the points in the list pl in four quadrants,
  --   using (cx,cy) as center.

  -- ON ENTRY :
  --   pl        a point list;
  --   cx        x-coordinate of the center of the point list;
  --   cy        y-coordinate of the center of the point list.

  -- ON RETURN :
  --   ne_cnt    number of points in the north east quadrant;
  --   nw_cnt    number of points in the north west quadrant;
  --   sw_cnt    number of points in the south west quadrant;
  --   se_cnt    number of points in the south east quadrant;
  --   ne_pts    list of points in north east quadrant;
  --   nw_pts    list of points in north west quadrant;
  --   sw_pts    list of points in south west quadrant;
  --   se_pts    list of points in south east quadrant.

-- CREATORS :

  function Create_Root_Leaf ( pl : Point_List ) return Link_to_Quad_Node;
  function Create_Root_Leaf ( pl : Point_List; size : natural32 )
                            return Link_to_Quad_Node;

  -- DESCRIPTION :
  --   Returns a pointer to a leaf quad node with the given point list.
  --   If size is omitted, then it will be computed as Length_Of(pl).

  function Create_Leaf ( pl : Point_List; size,depth : natural32 )
                       return Link_to_Quad_Node;

  -- DESCRIPTION :
  --   Returns a leaf node with the given contents.

  procedure Split_Leaf ( lqn : in out Link_to_Quad_Node );

  -- DESCRIPTION :
  --   Applies the Quadrant Partition procedure to the node,
  --   which must be a leaf.

  procedure Create ( lqn : in out Link_to_Quad_Node;
                     max_depth,min_size : in natural32 );

  -- DESCRIPTION :
  --   Recursive creation of a quad tree of depth at most max_depth
  --   and with splits only for nodes with size > min_size.

  function Create ( pl : Point_List;
                    max_depth,min_size : natural32 )
                  return Link_to_Quad_Node;
  function Create ( pl : Point_List; size : natural32;
                    max_depth,min_size : natural32 )
                  return Link_to_Quad_Node;

  -- DESCRIPTION :
  --   Returns a quad tree for the list, bounded by depth and size.

  -- ON ENTRY :
  --   pl        list of points in the plane;
  --   size      equals Length_Of(pl), may be omitted;
  --   max_depth is the maximal depth of the quad tree, i.e.: a leaf
  --             with depth equal to max_depth will not be splitted;
  --   min_size  is the minimal size of a leaf, i.e.: a leaf
  --             with size equal to min_size will not be splitted.

  -- ON RETURN :
  --   The root of a quad tree of depth at most max_depth.

-- SELECTORS :

  function Number_of_Leaves ( lqn : Link_to_Quad_Node ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of leaves under the given quad node.

  function Number_of_Nodes ( lqn : Link_to_Quad_Node ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of nodes under the given quad node.

  generic
    with procedure Process_Leaf ( lqn : in out Link_to_Quad_Node );
  procedure Enumerate_Leaves ( root : in out Link_to_Quad_Node );

  -- DESCRIPTION :
  --   Visits all nodes in depth first search mode and applies the
  --   procedure Process_Leaf to every leaf.

-- SORTING THE LEAVES :

  procedure Sort_Leaves ( root : in out Link_to_Quad_Node );

  -- DESCRIPTION :
  --   Sorts the leaves of the quad tree.

  generic
    with procedure Report ( lp1,lp2 : in Link_to_Point );
  procedure Clusters ( root : in Link_to_Quad_Node; tol : in double_float );

  -- DESCRIPTION :
  --   Searches the entire quad tree and reports the clusters.
  --   The leaves are expected to be sorted.

-- DESTRUCTORS :

  procedure Shallow_Clear ( lqn : in out Link_to_Quad_Node );

  -- DESCRIPTION :
  --   Releases the memory occupied by the pointer, which may lead
  --   to wasted memory if the data in the node is not released.

  procedure Deep_Clear ( lqn : in out Quad_Node );
  procedure Deep_Clear ( lqn : in out Link_to_Quad_Node );

  -- DESCRIPTION :
  --   Before releasing the memory occupied by the pointer,
  --   all data in the node is released as well.

end Standard_Quad_Trees;
