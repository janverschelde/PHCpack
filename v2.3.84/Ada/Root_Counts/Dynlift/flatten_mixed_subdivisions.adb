with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;

package body Flatten_Mixed_Subdivisions is

  procedure Flatten ( L : in out List ) is

    tmp : List := L;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if pt(pt'last) /= 0
       then pt(pt'last) := 0; Set_Head(tmp,pt);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Flatten;

  procedure Flatten ( L : in out Array_of_Lists ) is
  begin
    for i in L'range loop
      Flatten(L(i));
    end loop;
  end Flatten;

  procedure Flatten ( mic : in out Mixed_Cell ) is
  begin
    Flatten(mic.pts.all);
    mic.nor.all := (mic.nor'range => 0);
    mic.nor(mic.nor'last) := 1;
  end Flatten;

  procedure Old_Flatten ( mixsub : in out Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Flatten(mic);
      Set_Head(tmp,mic);
      tmp := Tail_Of(tmp);
    end loop;
  end Old_Flatten;

-- NEW FLATTENING, USING THE RECURSIVE DATA STRUCTURE :

  function Collect_Supports ( n : integer32; mixsub : Mixed_Subdivision )
                            return Array_of_Lists is

  -- DESCRIPTION :
  --   Returns the array of list of points that occur in the cells
  --   of the mixed subdivision.

  -- REQUIRED : not Is_Null(mixsub).

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell := Head_Of(mixsub);
    tmppts : List;
    pt : Link_to_Vector;
    res,res_last : Array_of_Lists(mic.pts'range);

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      for k in mic.pts'range loop
        tmppts := mic.pts(k);
        while not Is_Null(tmppts) loop
          pt := Head_Of(tmppts);    
          if not Is_In(res(k),pt)
           then Append(res(k),res_last(k),pt.all);
          end if;
          tmppts := Tail_Of(tmppts);
        end loop;
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Collect_Supports;

  procedure Flatten ( mixsub : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Flattens the mixed subdivision, i.e., the modified mixed subdivision
  --   contains one flattened cells with all the points that occured in the
  --   subdivision.  The original mixed subdivision is stored as the 
  --   subdivision of that flattened cell.

  begin
    if not Is_Null(mixsub) then
      declare
        n : constant integer32 := Head_Of(mixsub).nor'length-1;
        mic : Mixed_Cell;
        res : Mixed_Subdivision;
      begin
        mic.nor := new Standard_Integer_Vectors.Vector(1..n+1);
        mic.pts := new Array_of_Lists'(Collect_Supports(n,mixsub));
        Flatten(mic);
        mic.sub := new Mixed_Subdivision'(mixsub);
        Construct(mic,res);
        mixsub := res;
      end;
    end if;
  end Flatten;

end Flatten_Mixed_Subdivisions;
