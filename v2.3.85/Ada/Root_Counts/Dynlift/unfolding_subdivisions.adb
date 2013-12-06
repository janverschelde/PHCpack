with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Integer_Support_Functions;          use Integer_Support_Functions;
with Flatten_Mixed_Subdivisions;         use Flatten_Mixed_Subdivisions;

package body Unfolding_Subdivisions is

  function Different_Normals ( mixsub : Mixed_Subdivision ) return List is
  
    tmp : Mixed_Subdivision := mixsub;
    res,res_last : List;

  begin
    while not Is_Null(tmp) loop
      Append_Diff(res,res_last,Head_Of(tmp).nor.all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Different_Normals;

  function Extract ( normal : Vector; mixsub : Mixed_Subdivision )
                   return Mixed_Subdivision is

    tmp : Mixed_Subdivision := mixsub;
    res,res_last : Mixed_Subdivision;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        if mic.nor.all = normal
         then Append(res,res_last,mic);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Extract;

  function Merge_Same_Normal ( mixsub : Mixed_Subdivision )
                             return Mixed_Cell is

  -- DESCRIPTION :
  --   All cells with the same inner normal will be put in one cell,
  --   that will be contained in the mixed subdivision on return.

  -- REQUIRED :
  --   not Is_Null(mixsub) and all mixed cells have the same inner normal.

    tmp : Mixed_Subdivision;
    resmic,mic : Mixed_Cell;

  begin
    mic := Head_Of(mixsub);
    resmic.nor := new Standard_Integer_Vectors.Vector'(mic.nor.all);
    resmic.pts := new Array_of_Lists'(mic.pts.all);
    tmp := Tail_Of(mixsub);
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      declare
        last : List;
      begin
        for k in mic.pts'range loop
          last := resmic.pts(k);
          while not Is_Null(Tail_Of(last)) loop
            last := Tail_Of(last);
          end loop;
          Deep_Concat_Diff(resmic.pts(k),last,mic.pts(k));
        end loop;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return resmic;
  end Merge_Same_Normal;

  function Merge_Same_Normal ( mixsub : Mixed_Subdivision )
                             return Mixed_Subdivision is

  -- REQUIRED :
  --   not Is_Null(mixsub) and all mixed cells have the same inner normal.

    resmic : Mixed_Cell := Merge_Same_Normal(mixsub);
    ressub : Mixed_Subdivision;

  begin
    Construct(resmic,ressub);
    return ressub;
  end Merge_Same_Normal;

  function Merge ( mixsub : Mixed_Subdivision ) return Mixed_Subdivision is

  -- NOTE :
  --   Cells with an unique normal are simply taken over in the result,
  --   cells with the same normal are merged, hereby the refinement of these
  --   cells is destroyed.  Though, one could do better...

  begin
    if Is_Null(mixsub)
     then return mixsub;
     else
       declare
         tmp : Mixed_Subdivision := mixsub;
         res,res_last : Mixed_Subdivision;
         mic : Mixed_Cell;
       begin
         while not Is_Null(tmp) loop
           mic := Head_Of(tmp);
           if not Is_In(res,mic.nor.all)
            then
              if not Is_In(Tail_Of(tmp),mic.nor.all)
               then Append(res,res_last,mic);
               else declare
                      tmpmic : Mixed_Subdivision := Extract(mic.nor.all,tmp);
                      bigmic : Mixed_Cell := Merge_Same_Normal(tmpmic);
                    begin
                      Append(res,res_last,bigmic);
                    end;
              end if;
           end if;
           tmp := Tail_Of(tmp);
         end loop;
         return res;
       end;
    end if;
  end Merge;

  function Relift ( l : List; point : Vector ) return List is

    tmp,res : List;
    pt : Link_to_Vector;

  begin
    Copy(l,res);
    tmp := res;
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if pt.all = point
       then pt(pt'last) := 1;
       else pt(pt'last) := 0;
      end if;
      Set_Head(tmp,pt);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Relift;

  function Relift ( pts : Array_of_Lists; point : Vector )
                  return Array_of_Lists is

    res : Array_of_Lists(pts'range);

  begin
    for i in pts'range loop
      res(i) := Relift(pts(i),point);
    end loop;
    return res;
  end Relift;

  function Relift ( mic : Mixed_Cell; point : Vector ) return Mixed_Cell is

    res : Mixed_Cell;

  begin
    res.pts := new Array_of_Lists'(Relift(mic.pts.all,point));
    res.nor := new Standard_Integer_Vectors.Vector'(point'range => 0);
    Compute_Inner_Normal(res);
    return res;
  end Relift;

  function Relift ( mixsub : Mixed_Subdivision; point : Vector )
                  return Mixed_Subdivision is

    tmp,res,res_last : Mixed_Subdivision;

  begin
    tmp := mixsub;
    while not Is_Null(tmp) loop
      Append(res,res_last,Relift(Head_Of(tmp),point));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Relift;

  function Is_In_Point ( pt : Link_to_Vector; l : List ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the first n coordinates of pt belong to l.

    tmp : List := l;
    lpt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      if lpt(lpt'first..lpt'last-1) = pt(pt'first..pt'last-1)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In_Point;

  function Different_Points ( l1,l2 : List ) return natural32 is

  -- DESCRIPTION :
  --   Return the number of different points of the list l2 w.r.t. l1.

    res : natural32 := 0;
    tmp : List := l2;

  begin
    while not Is_Null(tmp) loop
      if not Is_In_Point(Head_Of(tmp),l1)
       then res := res + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Different_Points;

  function Different_Points ( l1,l2 : List ) return List is

  -- DESCRIPTION :
  --   Return the list of different points of the list l2 w.r.t. l1.

    res,res_last : List;
    tmp : List := l2;

  begin
    while not Is_Null(tmp) loop
      if not Is_In_Point(Head_Of(tmp),l1)
       then Append(res,res_last,Head_Of(tmp).all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Different_Points;

  function Different_Points ( pts : Array_of_Lists; mic : Mixed_Cell )
                            return natural32 is

  -- DESCRIPTION :
  --   Return the number of different points of the cell mic w.r.t. pts.

    res : natural32 := 0;

  begin
    for i in pts'range loop
      res := res + Different_Points(pts(i),mic.pts(i));
    end loop;
    return res;
  end Different_Points;

  function Different_Points ( pts : Array_of_Lists; mic : Mixed_Cell )
                            return Array_of_Lists is

  -- DESCRIPTION :
  --   Return the different points of the cell mic w.r.t. pts.

    res : Array_of_Lists(pts'range);

  begin
    for i in pts'range loop
      res(i) := Different_Points(pts(i),mic.pts(i));
    end loop;
    return res;
  end Different_Points;

  procedure Add ( l : in out List; pts : in List ) is

  -- DESCRIPTION :
  --   Adds the points in pts to l.

    tmp : List := pts;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      declare
        npt : Link_to_Vector := new Vector'(pt.all);
      begin
        Construct(npt,l);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Add;

  procedure Add ( l : in out Array_of_Lists; pts : in Array_of_Lists ) is

  -- DESCRIPTION :
  --   Adds the points in pts to l.

  begin
    for i in l'range loop
      Add(l(i),pts(i));
    end loop;
  end Add;

  procedure Put_Next_to_Front ( mixsub : in out Mixed_Subdivision;
                                pts : in Array_of_Lists ) is
 
  -- DESCRIPTION :
  --   Selects the next mixed cell to be processed, and puts in front
  --   of the list of cells mixsub.

    mic1 : Mixed_Cell := Head_Of(mixsub);
    min1 : natural32 := Different_Points(pts,mic1);
    tmp : Mixed_Subdivision := Tail_Of(mixsub);
    min : natural32;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      min := Different_Points(pts,mic);
      if min < min1 then
        min1 := min;
        Set_Head(mixsub,mic);
        Set_Head(tmp,mic1);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Put_Next_to_Front;

  procedure Relift ( l : in out List; ref : in List ) is

  -- DESCRIPTION :
  --   Gives all points in l, which belong to ref, lifting value 1.

    tmp : List := l;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if Is_In(ref,pt)
       then pt(pt'last) := 1;
       else pt(pt'last) := 0;
      end if;
      Set_Head(tmp,pt);
      tmp := Tail_Of(tmp);
    end loop;
  end Relift;

  procedure Relift ( l : in out List ) is

  -- DESCRIPTION :
  --   Gives all points lifting value 1.

    tmp : List := l;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      pt(pt'last) := 1;
      Set_Head(tmp,pt);
      tmp := Tail_Of(tmp);
    end loop;
  end Relift;

  procedure Relift ( l : in out Array_of_Lists; ref : in Array_of_Lists ) is

  -- DESCRIPTION :
  --   Gives all points in l, which belong to ref, lifting value 1.

  begin
    for i in l'range loop
      Relift(l(i),ref(i));
    end loop;
  end Relift;

  procedure Relift ( l : in out Array_of_Lists ) is

  -- DESCRIPTION :
  --   Gives all points lifting value 1.

  begin
    for i in l'range loop
      Relift(l(i));
    end loop;
  end Relift;

  procedure Relift ( mic : in out Mixed_Cell; pts : in out Array_of_Lists ) is

  -- DESCRIPTION :
  --   Gives the points in mic, which belong to pts lifting 1,
  --   and computes the new inner normal.

  begin
    Relift(mic.pts.all,pts);
    Relift(pts);
  end Relift;

  procedure Orientate_Inner_Normal 
                ( mic : in out Mixed_Cell; pts : in Array_of_Lists ) is

  -- DESCRIPTION :
  --   Orientates the normal of mic w.r.t. the points in pts.

    done : boolean := false;

  begin
    for i in pts'range loop
      if Minimal_Support(mic.pts(i),mic.nor.all) 
          > Minimal_Support(pts(i),mic.nor.all)
       then Min(mic.nor);
            done := true;
      end if;
      exit when done;
    end loop;
  end Orientate_Inner_Normal;

  procedure Unfolding ( mixsub : in out Mixed_Subdivision ) is

    tmp : Mixed_Subdivision;

  begin
    if not Is_Null(mixsub)
     then
       declare
         mic : Mixed_Cell := Head_Of(mixsub);
         pts : Array_of_Lists(mic.pts'range);
       begin
         Flatten(mic);
         Copy(mic.pts.all,pts);
         Process(mic,pts);
         tmp := Tail_Of(mixsub);
         while not Is_Null(tmp) loop
           Put_Next_to_Front(tmp,pts);
           mic := Head_Of(tmp);
           declare
             newpts : Array_of_Lists(pts'range);
           begin
             newpts := Different_Points(pts,mic);
             Relift(mic,newpts);
             Compute_Inner_Normal(mic);
            -- Orientate_Inner_Normal(mic,pts);
             Process(mic,newpts);
             Add(pts,newpts);
             Deep_Clear(newpts);
           end;
           tmp := Tail_Of(tmp);
         end loop;
       end;
    end if;
  end Unfolding;
  
end Unfolding_Subdivisions;
