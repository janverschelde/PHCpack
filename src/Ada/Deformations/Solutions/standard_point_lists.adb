with unchecked_deallocation;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package body Standard_Point_Lists is

-- HASH FUNCTIONS :

  function Complex_Inner_Product
             ( x,y : Standard_Complex_Vectors.Vector )
	     return Complex_Number is

    res : Complex_Number := Create(0.0);

  begin
    for i in x'range loop
      res := res + x(i)*Conjugate(y(i));
    end loop;
    return res;
  end Complex_Inner_Product;

-- CREATORS :

  function Create ( sv,h1,h2 : Standard_Complex_Vectors.Vector;
                    label : integer32 ) return Point is

    res : Point;
    ip1 : constant Complex_Number := Complex_Inner_Product(h1(sv'range),sv);
    ip2 : constant Complex_Number := Complex_Inner_Product(h2(sv'range),sv);

  begin
    res.label := label;
    res.x := REAL_PART(ip1) + IMAG_PART(ip1);
    res.y := REAL_PART(ip2) + IMAG_PART(ip2);
    return res;
  end Create;

  function Create ( ls : Link_to_Solution;
                    h1,h2 : Standard_Complex_Vectors.Vector;
                    label : integer32 ) return Point is
  begin
    return Create(ls.v,h1,h2,label);
  end Create;

  procedure Append ( first,last : in out Point_List;
                     h1,h2 : in Standard_Complex_Vectors.Vector;
                     label : in integer32;
                     sv : in Standard_Complex_Vectors.Vector ) is

    lp : constant Link_to_Point := new Point'(Create(sv,h1,h2,label));

  begin
    Append(first,last,lp);
  end Append;

  procedure Append ( first,last : in out Point_List;
                     h1,h2 : in Standard_Complex_Vectors.Vector;
                     label : in integer32; ls : in Link_to_Solution ) is

    lp : constant Link_to_Point := new Point'(Create(ls,h1,h2,label));

  begin
    Append(first,last,lp);
  end Append;

-- SELECTOR :

  procedure Center ( pl : in Point_List; cx,cy : out double_float ) is

    cnt : natural32 := 0;
    tmp : Point_List := pl;
    lp : Link_to_Point;

  begin
    cx := 0.0; cy := 0.0;
    while not Is_Null(tmp) loop
      lp := Head_Of(tmp);
      cnt := cnt + 1;
      cx := cx + lp.x;
      cy := cy + lp.y;
      tmp := Tail_Of(tmp);
    end loop;
    cx := cx/double_float(cnt);
    cy := cy/double_float(cnt);
  end Center;

-- SORTING :

  function "<" ( lp1,lp2 : Link_to_Point ) return boolean is
  begin
    if lp1.x < lp2.x then
      return true;
    elsif lp1.x > lp2.x then
      return false;
    elsif lp1.y < lp2.y then
      return true;
    else
      return false;
    end if;
  end "<";

  procedure Swap ( lp1,lp2 : in out Link_to_Point ) is

    b : integer32;
    z : double_float;

  begin
    b := lp1.label;
    lp1.label := lp2.label;
    lp2.label := b;
    z := lp1.x; lp1.x := lp2.x; lp2.x := z;
    z := lp1.y; lp1.y := lp2.y; lp2.y := z;
  end Swap;

  procedure Sort ( pl : in out Point_List ) is

    t1 : Point_List := pl;
    t2 : Point_List;
    lp1,lp2,min : Link_to_Point;

  begin
    while not Is_Null(t1) loop
      lp1 := Head_Of(t1);
      min := lp1;
      t2 := t1;
      while not Is_Null(t2) loop
        lp2 := Head_Of(t2);
        if lp2 < min
         then min := lp2;
        end if;
        t2 := Tail_Of(t2);
      end loop;
      if min /= lp1
       then Swap(min,lp1);
      end if;
      t1 := Tail_Of(t1);
    end loop;
  end Sort;

  procedure Insert ( pl : in out Point_List; pt : in Link_to_Point ) is

    first,second : Point_List;
    lpt : Link_to_Point;
    done : boolean;

  begin
    if Is_Null(pl) then                -- first special case: list is empty
      Construct(pt,pl);
    else
      first := pl;
      lpt := Head_Of(first);
      if pt < lpt then         -- second special case: first element larger
        Construct(pt,pl);
      else
        second := Tail_Of(first);    -- initialize pair of chasing pointers
        done := false;
        while not Is_Null(second) loop
          lpt := Head_Of(second);
          if lpt < pt then           -- point is larger than current lpt
            first := second;
            second := Tail_Of(second);
          else
            Construct(pt,second);    -- insert point into the list
            Swap_Tail(first,second);
            done := true;
          end if;
          exit when done;
        end loop;
        if not done then -- first points to the last element of the list
          Construct(pt,second);
          Swap_Tail(first,second);
        end if;
      end if;
    end if;
  end Insert;

  function Equal ( lp1,lp2 : Link_to_Point;
                   tol : double_float ) return boolean is
  begin
    if abs(lp1.x - lp2.x) > tol then
      return false;
    elsif abs(lp1.y - lp2.y) > tol then
      return false;
    else
      return true;
    end if;
  end Equal;

  procedure Insert_no_Duplicates
              ( pl : in out Point_List; pt : in Link_to_Point;
                tol : in double_float; lbl : out integer32 ) is

    first,second : Point_List;
    lpt : Link_to_Point;
    done : boolean;

  begin
    if Is_Null(pl) then                -- first special case: list is empty
      Construct(pt,pl);
      lbl := pt.label;
    else
      first := pl;
      lpt := Head_Of(first);
      if Equal(pt,lpt,tol) then  -- second special case: pt already in list
        lbl := lpt.label;
      elsif pt < lpt then       -- third special case: first element larger
        Construct(pt,pl);
        lbl := pt.label;
      else
        second := Tail_Of(first);    -- initialize pair of chasing pointers
        done := false;
        while not Is_Null(second) loop
          lpt := Head_Of(second);
          if Equal(pt,lpt,tol) then  -- point already belongs to the list
            lbl := lpt.label;
            done := true;
          elsif lpt < pt then        -- point is larger than current lpt
            first := second;
            second := Tail_Of(second);
          else
            Construct(pt,second);    -- insert point into the list
            Swap_Tail(first,second);
            lbl := pt.label;
            done := true;
          end if;
          exit when done;
        end loop;
        if not done then -- first points to the last element of the list
          Construct(pt,second);
          Swap_Tail(first,second);
          lbl := pt.label;
        end if;
      end if;
    end if;
  end Insert_no_Duplicates;

  procedure Insert_with_Duplicates
              ( pl : in out Point_List; pt : in Link_to_Point;
                tol : in double_float; cnt : out integer32;
                ptpl  : out Point_List ) is

    first,second : Point_List;
    lpt : Link_to_Point;
    done : boolean;

  begin
    if Is_Null(pl) then                -- first special case: list is empty
      Construct(pt,pl); ptpl := pl;
      cnt := 1;
    else
      first := pl;
      lpt := Head_Of(first);
      if Equal(pt,lpt,tol) then  -- second special case: pt already in list
        cnt := 2;                -- multiplicity is at least two
        loop
          first := Tail_Of(first);
          exit when Is_Null(first);
          cnt := cnt + 1;
        end loop;
        Construct(pt,pl); ptpl := pl;
      elsif pt < lpt then       -- third special case: first element larger
        Construct(pt,pl); ptpl := pl;
        cnt := 1;
      else
        second := Tail_Of(first);    -- initialize pair of chasing pointers
        done := false;
        while not Is_Null(second) loop
          lpt := Head_Of(second);
          if Equal(pt,lpt,tol) then  -- point already belongs to the list
            Construct(pt,second); ptpl := second;
            Swap_Tail(first,second);
            cnt := 2;                -- multiplicity is at least two
            loop
              second := Tail_Of(second);
              exit when Is_Null(second);
              cnt := cnt + 1;
            end loop;
            done := true;
          elsif lpt < pt then        -- point is larger than current lpt
            first := second;
            second := Tail_Of(second);
          else
            Construct(pt,second);    -- insert point into the list
            ptpl := second;
            Swap_Tail(first,second);
            cnt := 1;
            done := true;
          end if;
          exit when done;
        end loop;
        if not done then -- first points to the last element of the list
          Construct(pt,second); ptpl := second;
          Swap_Tail(first,second);
          cnt := 1;
        end if;
      end if;
    end if;
  end Insert_with_Duplicates;

  procedure Clusters ( pl : in Point_List; tol : in double_float ) is

    t1,t2 : Point_List;
    p1,p2 : Link_to_Point;

  begin
    if not Is_Null(pl) then
      t1 := pl; t2 := Tail_Of(t1);
      while not Is_Null(t2) loop
        p1 := Head_Of(t1); p2 := Head_Of(t2);
        if Equal(p1,p2,tol)
         then Report(p1,p2);
        end if;
        t1 := Tail_Of(t1); t2 := Tail_Of(t2);
      end loop;
    end if;
  end Clusters;

-- DESTRUCTORS :

  procedure free is new unchecked_deallocation(Point,Link_to_Point);

  procedure Clear ( lp : in out Link_to_Point ) is
  begin
    free(lp);
  end Clear;

  procedure Shallow_Clear ( pl : in out Point_List ) is
  begin
    List_of_Points.Clear(List_of_Points.List(pl));
  end Shallow_Clear;

  procedure Deep_Clear ( pl : in out Point_List ) is

    tmp : Point_List := pl;
    lp : Link_to_Point;

  begin
    while not Is_Null(tmp) loop
      lp := Head_Of(tmp);
      free(lp);
      tmp := Tail_Of(tmp);
    end loop;
    List_of_Points.Clear(List_of_Points.List(pl));
  end Deep_Clear;

end Standard_Point_Lists;
