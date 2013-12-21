with Standard_Integer64_Simplices;       use Standard_Integer64_Simplices;
with Global_Dynamic64_Triangulation;     use Global_Dynamic64_Triangulation;
with Dynamic64_Lifting_Functions;        use Dynamic64_Lifting_Functions;

package body Standard_Dynamic64_Triangulations is

-- UTILITIES :

  procedure Constant_Lifting
                ( L : in List; liftval : in integer64;
                  lifted,lifted_last : in out List ) is

    tmp : List := l;
    pt : Link_to_Vector;
 
  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      declare
        lpt : constant Link_to_Vector := new Vector(pt'first..pt'last+1);
      begin
        lpt(pt'range) := pt.all;
        lpt(lpt'last) := liftval;
        Append(lifted,lifted_last,lpt);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Constant_Lifting;

  procedure Initialize_Triangulation
                  ( L : in List; order : in boolean;
                    rest,lifted,lifted_last : in out List;
                    t : in out Triangulation ) is

  -- DESCRIPTION :
  --   Performs initialization of the dynamic lifting algorithm.

  -- ON ENTRY :
  --   L            the list of points to be processed;
  --   order        if true, then the order of the points will be considered;
  --   lifted       eventually points that already have been lifted;
  --   t            triangulation of the lifted points.

  -- ON RETURN :
  --   rest         rest of point list to be processed,
  --                if empty, then the problem is degenerate;
  --   lifted       points in the initial simplex if t was empty,
  --                otherwise left unchanged;
  --   lifted_last  pointer to the last element of lifted;
  --   t            if empty on entry, then it contains an initial simplex,
  --                if the problem was not degenerate.

    null_list : List;
    s : Simplex;

  begin
    if Is_Null(t) then
      Initial_Simplex(L,order,s,rest);    -- start from scratch
      if (s /= Null_Simplex) then
        Construct(s,t);
        lifted := Deep_Create(Vertices(s));
        lifted_last := lifted;
        while not Is_Null(Tail_Of(lifted_last)) loop
          lifted_last := Tail_Of(lifted_last);
        end loop;
      else
        rest := null_list;            -- degenerate problem
      end if;
    else
      rest := L;                          -- re-start
    end if;
  end Initialize_Triangulation;

  procedure Next_Point ( L : in out List; order : in boolean;
                         point : out Link_to_Vector ) is

  -- DESCRIPTION :
  --   Selects the next point out of the list L.

  -- ON ENTRY :
  --   L          a nonempty list of points;
  --   order      if true, then the next point in the list will be choosen,
  --              otherwise, the point will be picked randomly.

  -- ON RETURN :
  --   L          a list without the point;
  --   point      the choosen point.

    pt : Link_to_Vector;

  begin
    if order
     then pt := Head_Of(L);
     else pt := Max_Extreme(L,Head_Of(L)'last,-5,5); Swap_to_Front(L,pt);
    end if;
    L := Tail_Of(L);
    point := pt;
  end Next_Point;

-- BASIC VERSION : WITHOUT OUTPUT GENERICS :

  procedure Dynamic_Lifting
                ( L : in List; order,inter : in boolean;
                  maxli : in integer64; lifted,lifted_last : in out List;
                  t : in out Triangulation ) is

    rest,inner : List;
    pt : Link_to_Vector;
    nexli : integer64 := 1;

  begin
    Initialize_Triangulation(L,order,rest,lifted,lifted_last,t);
   -- ITERATE FOR ALL POINTS IN rest :
    while not Is_Null(rest) loop
      Next_Point(rest,order,pt);
      declare
        point : Link_to_Vector := new Vector(1..pt'last+1);
      begin
        point(1..pt'last) := pt.all;
        point(point'last) := 1; -- try to obtain an optimal lifting value !!
        if inter and then Is_In(t,point.all) then
          Clear(point); Construct(pt,inner);
        else
          nexli := Lift_to_Place(t,point.all);
          if (maxli > 0) and then (nexli > maxli)
           then Flatten(t); nexli := 1;
          end if;
          point(point'last) := nexli;
          Update(t,point);
          Append(lifted,lifted_last,point);
        end if;
      end;
    end loop;
    if inter                              -- lift out the interior points
     then Constant_Lifting(inner,nexli,lifted,lifted_last);
    end if;
  exception
    when constraint_error  -- probably due to a too high lifting
       => Flatten(t);
          Dynamic_Lifting(rest,order,inter,maxli,lifted,lifted_last,t);
  end Dynamic_Lifting;

  procedure Dynamic_Lifting_with_Flat
                ( L : in List; order,inter : in boolean;
                  maxli : in integer64; lifted,lifted_last : in out List;
                  t : in out Triangulation ) is

    rest,inner : List;
    pt : Link_to_Vector;
    nexli : integer64 := 1;

  begin
    Initialize_Triangulation(L,order,rest,lifted,lifted_last,t);
   -- ITERATE FOR ALL POINTS IN rest :
    while not Is_Null(rest) loop
      Next_Point(rest,order,pt);
      declare
        point : Link_to_Vector := new Vector(1..pt'last+1);
      begin
        point(1..pt'last) := pt.all;
        point(point'last) := 1; -- try to obtain an optimal lifting value !!
        if inter and then Is_In(t,point.all) then
          Clear(point); Construct(pt,inner);
        else
          nexli := Lift_to_Place(t,point.all);
          if (maxli > 0) and then (nexli > maxli)
           then Before_Flattening(t,lifted); Flatten(t); nexli := 1;
          end if;
          point(point'last) := nexli;
          Update(t,point);
          Append(lifted,lifted_last,point);
        end if;
      end;
    end loop;
    if inter                               -- lift out the interior points
     then Constant_Lifting(inner,nexli,lifted,lifted_last);
    end if;
  exception
    when constraint_error    -- probably due to a too high lifting
       => Before_Flattening(t,lifted); Flatten(t);
          Dynamic_Lifting_with_Flat
             (rest,order,inter,maxli,lifted,lifted_last,t);
  end Dynamic_Lifting_with_Flat;

  procedure Dynamic_Lifting_with_New
                ( L : in List; order,inter : in boolean;
                  maxli : in integer64; lifted,lifted_last : in out List;
                  t : in out Triangulation ) is

    rest,inner : List;
    pt : Link_to_Vector;
    nexli : integer64 := 1;
    newt : Triangulation;

  begin
    Initialize_Triangulation(L,order,rest,lifted,lifted_last,t);
   -- ITERATE FOR ALL POINTS IN rest :
    while not Is_Null(rest) loop
      Next_Point(rest,order,pt);
      declare
        point : Link_to_Vector := new Vector(1..pt'last+1);
      begin
        point(1..pt'last) := pt.all;
        point(point'last) := 1; -- try to obtain an optimal lifting value !!
        if inter and then Is_In(t,point.all) then
          Clear(point);  Construct(pt,inner);
        else
          nexli := Lift_to_Place(t,point.all);
          if (maxli > 0) and then (nexli > maxli)
           then Flatten(t); nexli := 1;
          end if;
          point(point'last) := nexli;
          Update(t,point,newt);
          Process_New_Simplices(newt,point.all);
          Append(lifted,lifted_last,point);
        end if;
      end;
    end loop;
    if inter                               -- lift out the interior points
     then Constant_Lifting(inner,nexli,lifted,lifted_last);
    end if;
  exception
    when constraint_error    -- probably due to a too high lifting
       => Flatten(t);
          Dynamic_Lifting_with_New(rest,order,inter,maxli,lifted,lifted_last,t);
  end Dynamic_Lifting_with_New;

  procedure Dynamic_Lifting_with_Flat_and_New
                ( L : in List; order,inter : in boolean;
                  maxli : in integer64; lifted,lifted_last : in out List;
                  t : in out Triangulation ) is

    rest,inner : List;
    pt : Link_to_Vector;
    nexli : integer64 := 1;
    newt : Triangulation;

  begin
    Initialize_Triangulation(L,order,rest,lifted,lifted_last,t);
   -- ITERATE FOR ALL POINTS IN rest :
    while not Is_Null(rest) loop
      Next_Point(rest,order,pt);
      declare
        point : Link_to_Vector := new Vector(1..pt'last+1);
      begin
        point(1..pt'last) := pt.all;
        point(point'last) := 1; -- try to obtain an optimal lifting value !!
        if inter and then Is_In(t,point.all) then
          Clear(point);  Construct(pt,inner);
        else
          nexli := Lift_to_Place(t,point.all);
          if (maxli > 0) and then (nexli > maxli)
           then Before_Flattening(t,lifted); Flatten(t); nexli := 1;
          end if;
          point(point'last) := nexli;
          Update(t,point,newt);
          Process_New_Simplices(newt,point.all);
          Append(lifted,lifted_last,point);
        end if;
      end;
    end loop;
    if inter                               -- lift out the interior points
     then Constant_Lifting(inner,nexli,lifted,lifted_last);
    end if;
  exception
    when constraint_error    -- probably due to a too high lifting
       => Before_Flattening(t,lifted); Flatten(t);
          Dynamic_Lifting_with_Flat_and_New
             (rest,order,inter,maxli,lifted,lifted_last,t);
  end Dynamic_Lifting_with_Flat_and_New;

end Standard_Dynamic64_Triangulations;
