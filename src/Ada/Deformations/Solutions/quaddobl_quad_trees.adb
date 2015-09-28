with unchecked_deallocation;

package body QuadDobl_Quad_Trees is

-- PARTITIONING INTO QUADRANTS :

  procedure Partition
               ( pl : in Point_List; cx,cy : in quad_double;
                 ne_cnt,nw_cnt,sw_cnt,se_cnt : out natural32;
                 ne_pts,nw_pts,sw_pts,se_pts : out Point_List ) is

    tmp : Point_List := pl;
    lp : Link_to_Point;
    ne_last,nw_last,sw_last,se_last : Point_List;

  begin
    ne_cnt := 0; nw_cnt := 0; sw_cnt := 0; se_cnt := 0;
    while not Is_Null(tmp) loop
      lp := Head_Of(tmp);
      if lp.x < cx then               -- point lies in the west
        if lp.y < cy then 
          sw_cnt := sw_cnt + 1;       -- point lies in the south west
          Append(sw_pts,sw_last,lp);
        else             
          nw_cnt := nw_cnt + 1;       -- point lies in the north west
          Append(nw_pts,nw_last,lp);
        end if;
      else                            -- point lies in the east
        if lp.y < cy then
          se_cnt := se_cnt + 1;       -- point lies in the south east 
          Append(se_pts,se_last,lp);
        else             
          ne_cnt := ne_cnt + 1;       -- point lies in the north east
          Append(ne_pts,ne_last,lp);
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Partition;

-- CREATORS :

  function Create_Root_Leaf ( pl : Point_List ) return Link_to_Quad_Node is
  begin
    return Create_Root_Leaf(pl,Length_Of(pl));
  end Create_Root_Leaf;

  function Create_Root_Leaf ( pl : Point_List; size : natural32 )
                            return Link_to_Quad_Node is

    res : Link_to_Quad_Node;
    qn : Quad_Node(true);

  begin
    qn.depth := 0;
    qn.size := size;
    qn.pts := pl;
    res := new Quad_Node'(qn);
    return res;
  end Create_Root_Leaf;

  function Create_Leaf ( pl : Point_List; size,depth : natural32 ) 
                       return Link_to_Quad_Node is

    res : Link_to_Quad_Node;
    qn : Quad_Node(true);

  begin
    qn.depth := depth;
    qn.size := size;
    qn.pts := pl;
    res := new Quad_Node'(qn);
    return res;
  end Create_Leaf;

  procedure Split_Leaf ( lqn : in out Link_to_Quad_Node ) is

    cx,cy : quad_double;
    res : Quad_Node(false);
    ne_cnt,nw_cnt,sw_cnt,se_cnt : natural32;
    ne_pts,nw_pts,sw_pts,se_pts : Point_List;

  begin
    if lqn.leaf then
      Center(lqn.pts,cx,cy);
      Partition(lqn.pts,cx,cy,ne_cnt,nw_cnt,sw_cnt,se_cnt,
                              ne_pts,nw_pts,sw_pts,se_pts);
      res.ne := Create_Leaf(ne_pts,ne_cnt,lqn.depth+1);
      res.nw := Create_Leaf(nw_pts,nw_cnt,lqn.depth+1);
      res.sw := Create_Leaf(sw_pts,sw_cnt,lqn.depth+1);
      res.se := Create_Leaf(se_pts,se_cnt,lqn.depth+1);
      res.depth := lqn.depth;
      res.size := lqn.size;
      Shallow_Clear(lqn.pts);
      Shallow_Clear(lqn);
      lqn := new Quad_Node'(res);
    end if;
  end Split_Leaf;

  procedure Create ( lqn : in out Link_to_Quad_Node;
                     max_depth,min_size : in natural32 ) is
  begin
    if lqn.depth < max_depth and lqn.size > min_size then
      Split_Leaf(lqn);
      Create(lqn.ne,max_depth+1,min_size);
      Create(lqn.nw,max_depth+1,min_size);
      Create(lqn.sw,max_depth+1,min_size);
      Create(lqn.se,max_depth+1,min_size);
    end if;
  end Create;

  function Create ( pl : Point_List;
                    max_depth,min_size : natural32 )
                  return Link_to_Quad_Node is
  begin
    return Create(pl,Length_Of(pl),max_depth,min_size);
  end Create;

  function Create ( pl : Point_List; size : natural32;
                    max_depth,min_size : natural32 )
                  return Link_to_Quad_Node is

    root : Link_to_Quad_Node := Create_Root_Leaf(pl,size);

  begin
    Create(root,max_depth,min_size);
    return root;
  end Create;

-- SELECTORS :

  function Number_of_Leaves ( lqn : Link_to_Quad_Node ) return natural32 is
  begin
    if lqn.leaf
     then return 1;
     else return Number_of_Leaves(lqn.ne) + Number_of_Leaves(lqn.nw)
               + Number_of_Leaves(lqn.sw) + Number_of_Leaves(lqn.se);
    end if;
  end Number_of_Leaves;

  function Number_of_Nodes ( lqn : Link_to_Quad_Node ) return natural32 is
  begin
    if lqn.leaf
     then return 1;
     else return Number_of_Nodes(lqn.ne) + Number_of_Nodes(lqn.nw)
               + Number_of_Nodes(lqn.sw) + Number_of_Nodes(lqn.se) + 1;
    end if;
  end Number_of_Nodes;

  procedure Enumerate_Leaves ( root : in out Link_to_Quad_Node ) is
  begin
    if root.leaf
     then Process_Leaf(root);
     else Enumerate_Leaves(root.ne); Enumerate_Leaves(root.nw);
          Enumerate_Leaves(root.sw); Enumerate_Leaves(root.se);
    end if;
  end Enumerate_Leaves;

-- SORTING :

  procedure Sort_Leaves ( root : in out Link_to_Quad_Node ) is
  
    procedure Sort_Leaf ( qn : in out Link_to_Quad_Node ) is
    begin
     -- put("Sorting "); put(qn.size,1); put_line(" points ...");
      Sort(qn.pts);
    end Sort_Leaf;
    procedure Sort_the_Leaves is new Enumerate_Leaves(Sort_Leaf);

  begin
    Sort_the_Leaves(root);
  end Sort_Leaves;

  procedure Clusters ( root : in Link_to_Quad_Node; tol : in double_float ) is

    procedure Report_Clusters is
      new QuadDobl_Point_Lists.Clusters(Report);

  begin
    if root.leaf
     then Report_Clusters(root.pts,tol);
     else Clusters(root.ne,tol); Clusters(root.nw,tol);
          Clusters(root.sw,tol); Clusters(root.se,tol);
    end if;
  end Clusters;

-- DESTRUCTORS :

  procedure free is new unchecked_deallocation(Quad_Node,Link_to_Quad_Node);

  procedure Shallow_Clear ( lqn : in out Link_to_Quad_Node ) is
  begin
    free(lqn);
  end Shallow_Clear;

  procedure Deep_Clear ( lqn : in out Quad_Node ) is
  begin
    if lqn.leaf
     then Deep_Clear(lqn.pts); lqn.size := 0;
     else Deep_Clear(lqn.ne); Deep_Clear(lqn.nw);
          Deep_Clear(lqn.sw); Deep_Clear(lqn.se);
    end if;
  end Deep_Clear;

  procedure Deep_Clear ( lqn : in out Link_to_Quad_Node ) is
  begin
    if lqn /= null
     then Deep_Clear(lqn.all); free(lqn);
    end if;
  end Deep_Clear;

end QuadDobl_Quad_Trees;
