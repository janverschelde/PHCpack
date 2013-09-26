with unchecked_deallocation;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Floating_Support_Functions;         use Floating_Support_Functions;
with Floating_Lifting_Utilities;

package body Floating_Mixed_Subdivisions is

-- CREATORS :

  function Compute_Inner_Normal
             ( n : integer32; pts : Array_of_Lists )
             return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..n) := (1..n => 0.0);
    tmp : List;
    pt,first : Standard_Floating_Vectors.Link_to_Vector;
    m1 : constant integer32 := n-1;
    A : Matrix(1..m1,1..m1);
    b : Standard_Floating_Vectors.Vector(1..m1);
    piv : Standard_Integer_Vectors.Vector(1..m1);
    ind : integer32 := 0;

  begin
    for i in pts'range loop
      first := Head_Of(pts(i));
      tmp := Tail_Of(pts(i));
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        ind := ind + 1;
        for j in 1..m1 loop
          A(ind,j) := pt(j) - first(j);
        end loop;
        b(ind) := first(n) - pt(n);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    lufac(A,m1,piv,ind);
    if ind = 0 then
      lusolve(A,m1,piv,b);
      res(b'range) := b;
      res(n) := 1.0;
    end if;
    return res;
  end Compute_Inner_Normal;

  function Recompute_Inner_Normal
              ( n : integer32; b : double_float; pts : Array_of_Lists )
              return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..n) := (1..n => 0.0);
    tmp : List;
    pt,first : Standard_Floating_Vectors.Link_to_Vector;
    m1 : constant integer32 := n-1;
    A : Matrix(1..m1,1..m1);
    w : Standard_Floating_Vectors.Vector(1..m1);
    piv : Standard_Integer_Vectors.Vector(1..m1);
    ind : integer32 := 0;

  begin
    for i in pts'range loop
      first := Head_Of(pts(i));
      tmp := Tail_Of(pts(i));
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        ind := ind + 1;
        for j in 1..m1 loop
          A(ind,j) := pt(j) - first(j);
        end loop;
        if first(n) = b then
          w(ind) := 10.0*first(n) - pt(n);
        elsif pt(n) = b then
          w(ind) := first(n) - 10.0*pt(n);
        else
          w(ind) := first(n) - pt(n);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    lufac(A,m1,piv,ind);
    if ind = 0 then
      lusolve(A,m1,piv,w);
      res(w'range) := w;
      res(n) := 1.0;
    end if;
    return res;
  end Recompute_Inner_Normal;

  function Create ( pts : Array_of_Lists;
                    nor : Standard_Floating_Vectors.Vector;
                    tol : double_float ) return Mixed_Cell is

    res : Mixed_Cell;
    sup : double_float;

  begin
    res.nor := new Standard_Floating_Vectors.Vector'(nor);
    res.pts := new Array_of_Lists(pts'range);
    for k in pts'range loop
      sup := Minimal_Support(pts(k),nor);
      res.pts(k) := Face(pts(k),nor,sup,tol);
    end loop;
    return res;
  end Create;

  function Create ( pts : Array_of_Lists; nors : List; tol : double_float )
                  return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : List := nors;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Create(pts,Head_Of(tmp).all,tol));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( pts : Array_of_Lists; mixsub : Mixed_Subdivision;
                    tol : double_float ) return Mixed_Subdivision is

    tmp,res,res_last : Mixed_Subdivision;

  begin tmp := mixsub;
    while not Is_Null(tmp) loop
      Append(res,res_last,Create(pts,Head_Of(tmp).nor.all,tol));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  procedure Update ( pts : in Array_of_Lists;
                     nor : in Standard_Floating_Vectors.Vector;
                     mixsub,mixsub_last : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Given a tuple of point sets and a normal,
  --   the mixed subdivision will be updated.

    tmp : Mixed_Subdivision := mixsub;
    done : boolean := false;

  begin
    while not Is_Null(tmp) and not done loop
      declare
        mic : constant Mixed_Cell := Head_Of(tmp);
        last : List;
      begin
        if Standard_Floating_Vectors.Equal(mic.nor.all,nor) then
          for k in mic.pts'range loop
            last := mic.pts(k);
            while not Is_Null(Tail_Of(last)) loop
              last := Tail_Of(last);
            end loop;
            Deep_Concat_Diff(mic.pts(k),last,pts(k));
          end loop;
          Set_Head(tmp,mic);
          done := true;
        else
          tmp := Tail_Of(tmp);
        end if;
      end;
    end loop;
    if not done then
      declare
        mic : Mixed_Cell;
      begin
        mic.pts := new Array_of_Lists(pts'range);
        Copy(pts,mic.pts.all);
        mic.nor := new Standard_Floating_Vectors.Vector'(nor);
        mic.sub := null;
        Append(mixsub,mixsub_last,mic);
      end;
    end if;
  end Update;

-- CONVERTORS :

  function Create ( s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Floating_VecVecs.Array_of_VecVecs is

    res : Standard_Floating_VecVecs.Array_of_VecVecs(s'range);
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in s'range loop
      res(i) := new Standard_Floating_VecVecs.Vecvec
                      (1..integer32(Length_Of(s(i))));
      tmp := s(i);
      for j in res(i)'range loop
        lpt := Head_Of(tmp);
        res(i)(j) := new Standard_Floating_Vectors.Vector'(lpt.all);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Create;

  function Position ( lifpts : Standard_Floating_VecVecs.VecVec;
                      point : Standard_Floating_Vectors.Vector )
                    return integer32 is
  begin
    for i in lifpts'range loop
      if Standard_Floating_Vectors.Equal(lifpts(i).all,point)
       then return i;
      end if;
    end loop;
    return 0;
  end Position;

  function Create_Labels
              ( pts : Standard_Floating_VecVecs.Array_of_VecVecs;
                mic : Mixed_Cell ) return Mixed_Labels is

    res : Mixed_Labels;
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    res.nor := new Standard_Floating_Vectors.Vector'(mic.nor.all);
    res.lab := new Standard_Integer_VecVecs.VecVec(mic.pts'range);
    for i in mic.pts'range loop
      res.lab(i) := new Standard_Integer_Vectors.Vector
                          (1..integer32(Length_Of(mic.pts(i))));
      tmp := mic.pts(i);
      for j in res.lab(i)'range loop
        lpt := Head_Of(tmp);
        res.lab(i)(j) := Position(pts(i).all,lpt.all);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    if mic.sub /= null then
      res.sub := new Mixed_Sublabeling'
                       (Create_Labeled_Subdivision(pts'last,mic.sub.all));
    end if;
    return res;
  end Create_Labels;

  function Create_Coordinates
              ( pts : Standard_Floating_VecVecs.Array_of_VecVecs;
                mlb : Mixed_Labels ) return Mixed_Cell is

    res : Mixed_Cell;
    last : List;

  begin
    res.nor := new Standard_Floating_Vectors.Vector'(mlb.nor.all);
    res.pts := new Arrays_of_Floating_Vector_Lists.Array_of_Lists(pts'range);
    for i in mlb.lab'range loop
      last := res.pts(i);
      for j in mlb.lab(i)'range loop
        Append(res.pts(i),last,pts(i)(mlb.lab(i)(j)).all);
      end loop;
    end loop;
    if mlb.sub /= null then
      res.sub := new Mixed_Subdivision'
                       (Create_Coordinate_Subdivision(pts'last,mlb.sub.all));
    end if;
    return res;
  end Create_Coordinates;

  function Create_Labeled_Subdivision
              ( r : integer32; sub : Mixed_Subdivision )
              return Mixed_Sublabeling is

    res : Mixed_Sublabeling;
    last : Lists_of_Mixed_Labels.List;
    ls : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r)
       := Floating_Lifting_Utilities.Lifted_Supports(r,sub);
    tmp : Mixed_Subdivision := sub;
    mic : Mixed_Cell;

  begin
    res.pts := new Standard_Floating_VecVecs.Array_of_VecVecs'(Create(ls));
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Lists_of_Mixed_Labels.Append
        (res.cells,last,Create_Labels(res.pts.all,mic));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create_Labeled_Subdivision;

  function Create_Coordinate_Subdivision
              ( r : integer32; sub : Mixed_Sublabeling )
              return Mixed_Subdivision is

    res,last : Mixed_Subdivision;
    mic : Mixed_Cell;
    tmp : Lists_of_Mixed_Labels.List := sub.cells;
    mlb : Mixed_Labels;

  begin
    while not Lists_of_Mixed_Labels.Is_Null(tmp) loop
      mlb := Lists_of_Mixed_Labels.Head_Of(tmp);
      mic := Create_Coordinates(sub.pts.all,mlb);
      Append(res,last,mic);
      tmp := Lists_of_Mixed_Labels.Tail_Of(tmp);
    end loop;
    return res;
  end Create_Coordinate_Subdivision;

-- CONSTRUCTORS :

  procedure Copy ( mic1 : in Mixed_Cell; mic2 : in out Mixed_Cell ) is

    use Standard_Floating_Vectors;

  begin
    Deep_Clear(mic2);
    if mic1.nor /= null
     then mic2.nor := new Standard_Floating_Vectors.Vector'(mic1.nor.all);
    end if;
    if mic1.pts /= null
     then mic2.pts := new Array_of_Lists(mic1.pts'range);
          Copy(mic1.pts.all,mic2.pts.all);
    end if;
    if mic1.sub /= null
     then mic2.sub := new Mixed_Subdivision;
          Copy(mic1.sub.all,mic2.sub.all);
    end if;
  end Copy;

  procedure Copy ( mixsub1 : in Mixed_Subdivision; 
                   mixsub2 : in out Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mixsub1;
    mixsub2_last : Mixed_Subdivision;

  begin
    Deep_Clear(mixsub2);
    while not Is_Null(tmp) loop
      declare
        mic1,mic2 : Mixed_Cell;
      begin
        mic1 := Head_Of(tmp);
        Copy(mic1,mic2);
        Append(mixsub2,mixsub2_last,mic2);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Copy;

  procedure Append_Diff ( first,last : in out Mixed_Subdivision;
                          mic : in Mixed_Cell ) is
  begin
    if not Is_In(first,mic)
     then Append(first,last,mic);
    end if;
  end Append_Diff;

  procedure Concat_Diff ( first,last : in out Mixed_Subdivision;
                          mixsub : in Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : constant Mixed_Cell := Head_Of(tmp);
      begin
        if not Is_In(first,mic)
         then Append_Diff(first,last,mic);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Concat_Diff;

  procedure Construct ( mixsub : in Mixed_Subdivision;
                        first : in out Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : constant Mixed_Cell := Head_Of(tmp);
      begin
        Construct(mic,first);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Construct;

  procedure Construct_Diff ( mixsub : in Mixed_Subdivision; 
                             first : in out Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : constant Mixed_Cell := Head_Of(tmp);
      begin
        if not Is_In(first,mic)
         then Construct(mic,first);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Construct_Diff;

-- SELECTORS :

  function Equal ( mic1,mic2 : Mixed_Cell ) return boolean is
  begin
    if not Standard_Floating_Vectors.Equal(mic1.nor,mic2.nor) then
      return false;
    elsif Equal(mic1.pts,mic2.pts) then
      return Equal(mic1.sub,mic2.sub);
    else 
      return false;
    end if;
  end Equal;

  function Is_Sub ( mixsub1,mixsub2 : Mixed_Subdivision ) return boolean is

  -- DESCRIPTION :
  --   Returns true when every cell in mixsub1 also belongs to mixsub2.

    tmp : Mixed_Subdivision := mixsub1;

  begin
    while not Is_Null(tmp) loop
      if not Is_In(mixsub2,Head_Of(tmp))
       then return false;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return true;
  end Is_Sub;

  function Equal ( mixsub1,mixsub2 : Mixed_Subdivision ) return boolean is
  begin
    if Is_Sub(mixsub1,mixsub2)
     then return Is_Sub(mixsub2,mixsub1);
     else return false;
    end if;
  end Equal;

  function Equal ( mixsub1,mixsub2 : Link_to_Mixed_Subdivision )
                 return boolean is
  begin
    if mixsub1 = null and then mixsub2 /= null then
      return false;
    elsif mixsub2 = null then
      return true;
    else
      return Equal(mixsub1.all,mixsub2.all);
    end if;
  end Equal;

  function Is_In ( mixsub : Mixed_Subdivision;
                   normal : Standard_Floating_Vectors.Vector )
                 return boolean is

    tmp : Mixed_Subdivision := mixsub;
    c : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      c := Head_Of(tmp);
      if Standard_Floating_Vectors.Equal(c.nor.all,normal)
       then return true;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return false;
  end Is_In;

  function Is_In ( mixsub : Mixed_Subdivision; mic : Mixed_Cell )
                 return boolean is

    tmp : Mixed_Subdivision := mixsub;
    mic1 : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic1 := Head_Of(tmp);
      if Equal(mic1,mic)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_Original ( mic : Mixed_Cell; b : double_float ) return boolean is

    n : constant integer32 := mic.nor'last;
    origin : Standard_Floating_Vectors.Vector(1..n) := (1..n => 0.0);

  begin
    origin(n) := b;
    for i in mic.pts'range loop
      if Is_In(mic.pts(i),origin)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Original;

  function Is_Stable
             ( nor : Standard_Floating_Vectors.Vector;
               b : double_float; pts : Array_of_Lists ) return boolean is

    newnor : Standard_Floating_Vectors.Vector(nor'range)
           := Recompute_Inner_Normal(nor'last,b,pts);

  begin
    for i in nor'range loop
      if nor(i) < 0.0 then
        if newnor(i) < nor(i) - 1.0
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Is_Stable;

  function Zero_Type
             ( nor : Standard_Floating_Vectors.Vector;
               b : double_float; pts : Array_of_Lists )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(nor'first..nor'last-1);
    newnor : Standard_Floating_Vectors.Vector(nor'range)
           := Recompute_Inner_Normal(nor'last,b,pts);

  begin
    for i in res'range loop
      if nor(i) = 0.0 then
        res(i) := 1;
      elsif nor(i) < 0.0 then
        if newnor(i) < nor(i) - 1.0 then
          res(i) := -1;
        else
          res(i) := +1;
        end if;
      elsif nor(i) > 0.0 then
        if newnor(i) > nor(i) + 1.0 then
          res(i) := 0;
        else
          res(i) := +1;
        end if;
      end if;
    end loop;
    return res;
  end Zero_Type;

  function Is_Stable ( mic : Mixed_Cell; b : double_float ) return boolean is
  begin
    if Is_Original(mic,b)
     then return true;
     else return Is_Stable(mic.nor.all,b,mic.pts.all);
    end if;
  end Is_Stable;

  procedure Split_Original_Cells
              ( mixsub : in Mixed_Subdivision; b : in double_float;
                orgmcc,stbmcc : out Mixed_Subdivision;
                orgcnt,stbcnt : out natural32 ) is

    orglast,stblast : Mixed_Subdivision;
    mic : Mixed_Cell;
    tmp : Mixed_Subdivision := mixsub;
    cntorg,cntstb : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      if Is_Original(mic,b) then
        Append(orgmcc,orglast,mic);
        cntorg := cntorg + 1; 
      elsif Is_Stable(mic.nor.all,b,mic.pts.all) then
        Append(stbmcc,stblast,mic);
        cntstb := cntstb + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    orgcnt := cntorg; stbcnt := cntstb;
  end Split_Original_Cells;

-- DESTRUCTORS :

  procedure free is 
    new unchecked_deallocation(Mixed_Subdivision,Link_to_Mixed_Subdivision);
  procedure free is 
    new unchecked_deallocation(Mixed_Sublabeling,Link_to_Mixed_Sublabeling);

  procedure Deep_Clear ( mic : in out Mixed_Cell ) is
  begin
    Standard_Floating_Vectors.Clear(mic.nor);
    Deep_Clear(mic.pts); Deep_Clear(mic.sub);
  end Deep_Clear;

  procedure Shallow_Clear ( mic : in out Mixed_Cell ) is
  begin
    Standard_Floating_Vectors.Clear(mic.nor);
    Shallow_Clear(mic.pts); Shallow_Clear(mic.sub);
  end Shallow_Clear;

  procedure Clear ( mlb : in out Mixed_Labels ) is
  begin
    Standard_Floating_Vectors.Clear(mlb.nor);
    Standard_Integer_VecVecs.Deep_Clear(mlb.lab);
    Clear(mlb.sub);
  end Clear;

  procedure Deep_Clear ( mixsub : in out Mixed_Subdivision ) is

    tmp : Mixed_Subdivision;

  begin
    tmp := mixsub;
    while not Is_Null(tmp) loop
      declare
	mic : Mixed_Cell := Head_Of(tmp);
      begin
	Deep_Clear(mic);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Shallow_Clear(mixsub);
  end Deep_Clear;

  procedure Deep_Clear ( mixsub : in out Link_to_Mixed_Subdivision ) is
  begin
    if mixsub /= null
     then Deep_Clear(mixsub.all); free(mixsub);
    end if;
  end Deep_Clear;

  procedure Shallow_Clear ( mixsub : in out Mixed_Subdivision ) is
  begin
    Lists_of_Mixed_Cells.Clear(Lists_of_Mixed_Cells.List(mixsub));
  end Shallow_Clear;

  procedure Shallow_Clear ( mixsub : in out Link_to_Mixed_Subdivision ) is
  begin
    if mixsub /= null
     then Shallow_Clear(mixsub.all); free(mixsub);
    end if;
  end Shallow_Clear;

  procedure Clear ( mixsub : in out Mixed_Sublabeling ) is

    tmp : Lists_of_Mixed_Labels.List := mixsub.cells;
    mlb : Mixed_Labels;

  begin
    Standard_Floating_VecVecs.Deep_Clear(mixsub.pts);
    while not Lists_of_Mixed_Labels.Is_Null(tmp) loop
      mlb := Lists_of_Mixed_Labels.Head_Of(tmp);
      Clear(mlb);
      tmp := Lists_of_Mixed_Labels.Tail_Of(tmp);
    end loop;
    Lists_of_Mixed_Labels.Clear(mixsub.cells);
  end Clear;

  procedure Clear ( mixsub : in out Link_to_Mixed_Sublabeling ) is
  begin
    if mixsub /= null
     then Clear(mixsub.all); free(mixsub);
    end if;
  end Clear;

end Floating_Mixed_Subdivisions;
