with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Standard_Point_Lists;              use Standard_Point_Lists;

procedure ts_ptlists is

-- DESCRIPTION :
--   Tests on sorting of point lists.

  procedure put ( pl : in Point_List ) is

  -- DESCRIPTION :
  --   Writes the list of points to screen.

    tmp : Point_List := pl;
    lpt : Link_to_Point;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      put(lpt.label,1); put(" : ");
      put(lpt.x); put("  "); put(lpt.y); new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a number and then
  --   generates a lists of that many random points.

    first,last,pl : Point_List;
    n : integer32 := 0;

  begin
    put("Give the number of points : "); get(n);
    for i in 1..n loop
      declare
        pt : Point;
        lpt1,lpt2 : Link_to_Point;
      begin
        pt.label := i;
        pt.x := Standard_Random_Numbers.Random;
        pt.y := Standard_Random_Numbers.Random;
        lpt1 := new Point'(pt);
        lpt2 := new Point'(pt);
        Append(first,last,lpt1);
        Insert(pl,lpt2);
      end;
    end loop;
    put("A list of "); put(n,1); put_line(" random points :");
    put(first);
    Sort(first);
    put_line("The sorted list :"); put(first);
    put_line("The list constructed with insert:"); put(pl);
  end Main;

begin
  Main;
end ts_ptlists;
