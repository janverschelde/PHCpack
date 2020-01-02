with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Point_Lists;
with Standard_Quad_Trees;
with DoblDobl_Point_Lists;
with DoblDobl_Quad_Trees;
with QuadDobl_Point_Lists;
with QuadDobl_Quad_Trees;

procedure ts_quadpart is
 
-- DESCRIPTION :
--   Tests the partition procedure for the construction of quad trees.

  function Random_Points
             ( n : integer32 ) return Standard_Point_Lists.Point_List is

  -- DECRIPTION :
  --   Returns a list of n randomly generated points,
  --   with coordinates uniformly distributed in [-1, +1].

    use Standard_Point_Lists;
    res,res_last : Point_List;

  begin
    for i in 1..n loop
      declare
        pt : Point;
        lpt : Link_to_Point;
      begin
        pt.label := i;
        pt.x := Standard_Random_Numbers.Random;
        pt.y := Standard_Random_Numbers.Random;
        lpt := new Point'(pt);
        Append(res,res_last,lpt);
      end;
    end loop;
    return res;
  end Random_Points;

  function Random_Points
             ( n : integer32 ) return DoblDobl_Point_Lists.Point_List is

  -- DECRIPTION :
  --   Returns a list of n randomly generated points,
  --   with coordinates uniformly distributed in [-1, +1].

    use DoblDobl_Point_Lists;
    res,res_last : Point_List;

  begin
    for i in 1..n loop
      declare
        pt : Point;
        lpt : Link_to_Point;
      begin
        pt.label := i;
        pt.x := DoblDobl_Random_Numbers.Random;
        pt.y := DoblDobl_Random_Numbers.Random;
        lpt := new Point'(pt);
        Append(res,res_last,lpt);
      end;
    end loop;
    return res;
  end Random_Points;

  function Random_Points
             ( n : integer32 ) return QuadDobl_Point_Lists.Point_List is

  -- DECRIPTION :
  --   Returns a list of n randomly generated points,
  --   with coordinates uniformly distributed in [-1, +1].

    use QuadDobl_Point_Lists;
    res,res_last : Point_List;

  begin
    for i in 1..n loop
      declare
        pt : Point;
        lpt : Link_to_Point;
      begin
        pt.label := i;
        pt.x := QuadDobl_Random_Numbers.Random;
        pt.y := QuadDobl_Random_Numbers.Random;
        lpt := new Point'(pt);
        Append(res,res_last,lpt);
      end;
    end loop;
    return res;
  end Random_Points;

  procedure Quadrant_Partition ( pl : in Standard_Point_Lists.Point_List ) is

  -- DESCRIPTION :
  --   Partitions the list pl into four quadrants,
  --   where the center of the quadrant is the point in the list
  --   with average coordinates.

    use Standard_Point_Lists;
    use Standard_Quad_Trees;

    timer : Timing_Widget;
    cx,cy : double_float;
    ne_cnt,nw_cnt,sw_cnt,se_cnt : natural32;
    ne_pts,nw_pts,sw_pts,se_pts : Point_List;

  begin
    tstart(timer);
    Center(pl,cx,cy);
    tstop(timer);
    new_line;
    put("Center : ("); put(cx); put(","); put(cy); put_line(").");
    new_line;
    print_times(standard_output,timer,"computing the center");
    tstart(timer);
    Partition(pl,cx,cy,ne_cnt,nw_cnt,sw_cnt,se_cnt,
                       ne_pts,nw_pts,sw_pts,se_pts);
    tstop(timer);
    new_line;
    put_line("Distribution of the points into quadrants : ");
    put("  #points in north east quadrant : "); put(ne_cnt,1); new_line;
    put("  #points in north west quadrant : "); put(nw_cnt,1); new_line;
    put("  #points in south west quadrant : "); put(sw_cnt,1); new_line;
    put("  #points in south east quadrant : "); put(se_cnt,1); new_line;
    new_line;
    print_times(standard_output,timer,"partitioning into quadrants");
  end Quadrant_Partition;

  procedure Quadrant_Partition ( pl : in DoblDobl_Point_Lists.Point_List ) is

  -- DESCRIPTION :
  --   Partitions the list pl into four quadrants,
  --   where the center of the quadrant is the point in the list
  --   with average coordinates.

    use DoblDobl_Point_Lists;
    use DoblDobl_Quad_Trees;

    timer : Timing_Widget;
    cx,cy : double_double;
    ne_cnt,nw_cnt,sw_cnt,se_cnt : natural32;
    ne_pts,nw_pts,sw_pts,se_pts : Point_List;

  begin
    tstart(timer);
    Center(pl,cx,cy);
    tstop(timer);
    new_line;
    put("Center : ("); put(cx); put(","); put(cy); put_line(").");
    new_line;
    print_times(standard_output,timer,"computing the center");
    tstart(timer);
    Partition(pl,cx,cy,ne_cnt,nw_cnt,sw_cnt,se_cnt,
                       ne_pts,nw_pts,sw_pts,se_pts);
    tstop(timer);
    new_line;
    put_line("Distribution of the points into quadrants : ");
    put("  #points in north east quadrant : "); put(ne_cnt,1); new_line;
    put("  #points in north west quadrant : "); put(nw_cnt,1); new_line;
    put("  #points in south west quadrant : "); put(sw_cnt,1); new_line;
    put("  #points in south east quadrant : "); put(se_cnt,1); new_line;
    new_line;
    print_times(standard_output,timer,"partitioning into quadrants");
  end Quadrant_Partition;

  procedure Quadrant_Partition ( pl : in QuadDobl_Point_Lists.Point_List ) is

  -- DESCRIPTION :
  --   Partitions the list pl into four quadrants,
  --   where the center of the quadrant is the point in the list
  --   with average coordinates.

    use QuadDobl_Point_Lists;
    use QuadDobl_Quad_Trees;

    timer : Timing_Widget;
    cx,cy : quad_double;
    ne_cnt,nw_cnt,sw_cnt,se_cnt : natural32;
    ne_pts,nw_pts,sw_pts,se_pts : Point_List;

  begin
    tstart(timer);
    Center(pl,cx,cy);
    tstop(timer);
    new_line;
    put("Center : ("); put(cx); put(","); put(cy); put_line(").");
    new_line;
    print_times(standard_output,timer,"computing the center");
    tstart(timer);
    Partition(pl,cx,cy,ne_cnt,nw_cnt,sw_cnt,se_cnt,
                       ne_pts,nw_pts,sw_pts,se_pts);
    tstop(timer);
    new_line;
    put_line("Distribution of the points into quadrants : ");
    put("  #points in north east quadrant : "); put(ne_cnt,1); new_line;
    put("  #points in north west quadrant : "); put(nw_cnt,1); new_line;
    put("  #points in south west quadrant : "); put(sw_cnt,1); new_line;
    put("  #points in south east quadrant : "); put(se_cnt,1); new_line;
    new_line;
    print_times(standard_output,timer,"partitioning into quadrants");
  end Quadrant_Partition;

  procedure Standard_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Generates n random points and partitions the collection.

    use Standard_Point_Lists;

    pl : constant Point_List := Random_Points(n);

  begin
    Quadrant_Partition(pl);
  end Standard_Test;

  procedure DoblDobl_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Generates n random points and partitions the collection.

    use DoblDobl_Point_Lists;

    pl : constant Point_List := Random_Points(n);

  begin
    Quadrant_Partition(pl);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Generates n random points and partitions the collection.

    use QuadDobl_Point_Lists;

    pl : constant Point_List := Random_Points(n);

  begin
    Quadrant_Partition(pl);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of points and
  --   the precision before calling the proper test procedure.

    nb : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Creating quad trees for random point lists ...");
    new_line;
    put("Give the number of points : "); get(nb);
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test(nb);
      when '1' => DoblDobl_Test(nb);
      when '2' => QuadDobl_Test(nb);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_quadpart;
