with unchecked_deallocation;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Integer_Support_Functions;          use Integer_Support_Functions;
with Standard_Integer_Norms;             use Standard_Integer_Norms;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;

package body Integer_Mixed_Subdivisions is

-- CREATORS :

  procedure Compute_Inner_Normal ( mic : in out Mixed_Cell ) is

    len : constant integer32 
        := integer32(Length_Of(mic.pts.all)) - mic.pts'length;
    im : matrix(1..len,mic.nor'range);
    tmp : List;
    pt,first : Link_to_Vector;
    cnt : integer32 := 0;

  begin
    for i in mic.pts'range loop               -- compute the inner normal
      first := Head_Of(mic.pts(i));
      tmp := Tail_Of(mic.pts(i));
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        cnt := cnt + 1;
        for j in im'range(2) loop
          im(cnt,j) := pt(j) - first(j);
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    Upper_Triangulate(im);
    Scale(im);
    Solve0(im,mic.nor.all);
    Normalize(mic.nor.all);
    if mic.nor(mic.nor'last) < 0             -- orientate the normal
     then Min(mic.nor);
    end if;
  end Compute_Inner_Normal;

  function Create ( pts : Array_of_Lists; nor : Vector ) return Mixed_Cell is

    res : Mixed_Cell;
    sup : integer32;

  begin
    res.nor := new Vector'(nor);
    res.pts := new Array_of_Lists(pts'range);
    for k in pts'range loop
      sup := Minimal_Support(pts(k),nor);
      res.pts(k) := Face(pts(k),nor,sup);
    end loop;
    return res;
  end Create;

  function Create ( pts : Array_of_Lists; nors : List )
                  return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : List := nors;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Create(pts,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( pts : Array_of_Lists; mixsub : Mixed_Subdivision )
                  return Mixed_Subdivision is

    tmp,res,res_last : Mixed_Subdivision;

  begin
    tmp := mixsub;
    while not Is_Null(tmp) loop
      Append(res,res_last,Create(pts,Head_Of(tmp).nor.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  procedure Update ( pts : in Array_of_Lists; nor : in Vector;
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
        if Equal(mic.nor.all,nor) then
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
        mic.nor := new Standard_Integer_Vectors.Vector'(nor);
        mic.sub := null;
        Append(mixsub,mixsub_last,mic);
      end;
    end if;
  end Update;

  procedure Update ( mixsub,mixsub_last : in out Mixed_Subdivision;
                     cells : in Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := cells;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Update(mic.pts.all,mic.nor.all,mixsub,mixsub_last);
      tmp := Tail_Of(tmp);
    end loop;
  end Update;

-- CONSTRUCTORS :

  procedure Copy ( mic1 : in Mixed_Cell; mic2 : in out Mixed_Cell ) is
  begin
    Deep_Clear(mic2);
    if mic1.nor /= null
     then mic2.nor := new Standard_Integer_Vectors.Vector'(mic1.nor.all);
    end if;
    if mic1.pts /= null
     then mic2.pts := new Array_of_Lists(mic1.pts'range);
          Copy(mic1.pts.all,mic2.pts.all);
    end if;
    if mic1.sub /= null then
      mic2.sub := new Mixed_Subdivision;
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
    if not Equal(mic1.nor,mic2.nor) then
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

  function Is_In ( mixsub : Mixed_Subdivision; normal : Vector )
                 return boolean is

    tmp : Mixed_Subdivision := mixsub;
    c : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      c := Head_Of(tmp);
      if Equal(c.nor.all,normal)
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

-- DESTRUCTORS :

  procedure free is new unchecked_deallocation
      (Mixed_Subdivision,Link_to_Mixed_Subdivision);

  procedure Deep_Clear ( mic : in out Mixed_Cell ) is
  begin
    Clear(mic.nor); Deep_Clear(mic.pts); Deep_Clear(mic.sub);
  end Deep_Clear;

  procedure Shallow_Clear ( mic : in out Mixed_Cell ) is
  begin
    Clear(mic.nor); Shallow_Clear(mic.pts); Shallow_Clear(mic.sub);
  end Shallow_Clear;

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

end Integer_Mixed_Subdivisions;
