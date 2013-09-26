with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Permutations;                       use Permutations;
with Permute_Operations;                 use Permute_Operations;

package body Generating_Mixed_Cells is

-- FIRST TARGET ROUTINE :

  function Permute ( L : Lists_of_Integer_Vectors.List; p : Permutation )
                   return Lists_of_Integer_Vectors.List is

  -- DESCRIPTION :
  --   Applies the permutation p to all elements in the list l.

    use Lists_of_Integer_Vectors;

    tmp,res,res_last : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      declare
        plv,lv : Standard_Integer_Vectors.Link_to_Vector;
      begin
        lv := Head_Of(tmp);
        plv := new Standard_Integer_Vectors.Vector'(p*lv.all);
        plv(plv'last) := lv(lv'last);  -- same lifting !
        Append(res,res_last,plv.all);
        Standard_Integer_Vectors.Clear(plv);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Permute;

  function Permute ( p : Permutation; mix : Vector;
                     mic : Integer_Mixed_Subdivisions.Mixed_Cell )
                   return Integer_Mixed_Subdivisions.Mixed_Cell is

  -- DESCRIPTION :
  --   Permutes the components of mic.pts according to the permutation p.

    use Lists_of_Integer_Vectors;
    use Arrays_of_Integer_Vector_Lists;
    use Integer_Mixed_Subdivisions;

    res : Mixed_Cell;
    index : integer32;
 
  begin
    res.nor := new Standard_Integer_Vectors.Vector'(mic.nor.all);
    res.pts := new Array_of_Lists(mic.pts'range);
    for k in res.pts'range loop
      index := Compute_Index(p(k),mix);
      Copy(mic.pts(index),res.pts(k));
    end loop;
    return res;
  end Permute;

  function Permute ( mic : Integer_Mixed_Subdivisions.Mixed_Cell;
                     p : Permutation )
                   return Integer_Mixed_Subdivisions.Mixed_Cell is
 
  -- DESCRIPTION :
  --   Applies permutation p on the mixed cell mic.

    use Arrays_of_Integer_Vector_Lists;
    use Integer_Mixed_Subdivisions;

    res : Mixed_Cell;

  begin
    res.nor := new Standard_Integer_Vectors.Vector'(p*mic.nor.all);
    res.nor(res.nor'last) := mic.nor(mic.nor'last);
    res.pts := new Array_of_Lists(mic.pts'range);
    for k in mic.pts'range loop
      res.pts(k) := Permute(mic.pts(k),p);
    end loop;
    return res;
  end Permute;

  procedure Permute_and_Append
               ( v,w : in List_of_Permutations;
                 mic : in Integer_Mixed_Subdivisions.Mixed_Cell;
                 mix : in Vector;
                 mixsub,mixsub_last
                     : in out Integer_Mixed_Subdivisions.Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Applies all permutations to the mixed cell and appends the results.

    use Integer_Mixed_Subdivisions;

    lv,lw : List_of_Permutations;

  begin
    lv := v; lw := w;
    while not Is_Null(lv) loop
      declare
        vmic,wmic : Mixed_Cell;
      begin
        vmic := Permute(mic,Permutation(Head_Of(lv).all));
        wmic := Permute(Permutation(Head_Of(lw).all),mix,vmic);
        Deep_Clear(vmic);
        if not Is_In(mixsub,wmic.nor.all)
         then Append(mixsub,mixsub_last,wmic);
        end if;
      end;
      lv := Tail_Of(lv);
      lw := Tail_Of(lw);
    end loop;
  end Permute_and_Append;

  function Generating_Cells
              ( v,w : List_of_Permutations; mix : Vector;
                mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision )
              return Integer_Mixed_Subdivisions.Mixed_Subdivision is

    use Integer_Mixed_Subdivisions;

    tmp,res,res_last,done,done_last : Mixed_Subdivision;

  begin
    tmp := mixsub;
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        if not Is_In(done,mic.nor.all)
         then Append(res,res_last,mic);
              Permute_and_Append(v,w,mic,mix,done,done_last);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(done);
    return res;
  end Generating_Cells;

-- SECOND TARGET ROUTINE :

  function Permute ( l : Lists_of_Floating_Vectors.List; p : Permutation )
                   return Lists_of_Floating_Vectors.List is

  -- DESCRIPTION :
  --   Applies the permutation p to all elements in the list l.

    use Lists_of_Floating_Vectors;

    tmp,res,res_last : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      declare
        plv,lv : Standard_Floating_Vectors.Link_to_Vector;
      begin
        lv := Head_Of(tmp);
        plv := new Standard_Floating_Vectors.Vector'(p*lv.all);
        plv(plv'last) := lv(lv'last);  -- same lifting !
        Append(res,res_last,plv.all);
        Standard_Floating_Vectors.Clear(plv);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Permute;

  function Permute ( p : Permutation; mix : Vector;
                     mic : Floating_Mixed_Subdivisions.Mixed_Cell )
                   return Floating_Mixed_Subdivisions.Mixed_Cell is

  -- DESCRIPTION :
  --   Permutes the components of mic.pts according to the permutation p.

    use Lists_of_Floating_Vectors;
    use Arrays_of_Floating_Vector_Lists;
    use Floating_Mixed_Subdivisions;

    res : Mixed_Cell;
    index : integer32;
 
  begin
    res.nor := new Standard_Floating_Vectors.Vector'(mic.nor.all);
    res.pts := new Array_of_Lists(mic.pts'range);
    for k in res.pts'range loop
      index := Compute_Index(p(k),mix);
      Copy(mic.pts(index),res.pts(k));
    end loop;
    return res;
  end Permute;

  function Permute ( mic : Floating_Mixed_Subdivisions.Mixed_Cell;
                     p : Permutation )
                   return Floating_Mixed_Subdivisions.Mixed_Cell is
 
  -- DESCRIPTION :
  --   Applies permutation p on the mixed cell mic.

    use Arrays_of_Floating_Vector_Lists;
    use Floating_Mixed_Subdivisions;

    res : Mixed_Cell;

  begin
    res.nor := new Standard_Floating_Vectors.Vector'(p*mic.nor.all);
    res.nor(res.nor'last) := mic.nor(mic.nor'last);
    res.pts := new Array_of_Lists(mic.pts'range);
    for k in mic.pts'range loop
      res.pts(k) := Permute(mic.pts(k),p);
    end loop;
    return res;
  end Permute;

  procedure Permute_and_Append
               ( v,w : in List_of_Permutations;
                 mic : in Floating_Mixed_Subdivisions.Mixed_Cell;
                 mix : in Vector;
                 mixsub,mixsub_last
                     : in out Floating_Mixed_Subdivisions.Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Applies all permutations to the mixed cell and appends the results.

    use Floating_Mixed_Subdivisions;

    lv,lw : List_of_Permutations;

  begin
    lv := v; lw := w;
    while not Is_Null(lv) loop
      declare
        vmic,wmic : Mixed_Cell;
      begin
        vmic := Permute(mic,Permutation(Head_Of(lv).all));
        wmic := Permute(Permutation(Head_Of(lw).all),mix,vmic);
        Deep_Clear(vmic);
        if not Is_In(mixsub,wmic.nor.all)
         then Append(mixsub,mixsub_last,wmic);
        end if;
      end;
      lv := Tail_Of(lv);
      lw := Tail_Of(lw);
    end loop;
  end Permute_and_Append;

  function Generating_Cells
              ( v,w : List_of_Permutations; mix : Vector;
                mixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision )
              return Floating_Mixed_Subdivisions.Mixed_Subdivision is

    use Floating_Mixed_Subdivisions;

    tmp,res,res_last,done,done_last : Mixed_Subdivision;

  begin
    tmp := mixsub;
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        if not Is_In(done,mic.nor.all) then
          Append(res,res_last,mic);
          Permute_and_Append(v,w,mic,mix,done,done_last);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(done);
    return res;
  end Generating_Cells;

-- THIRD TARGET ROUTINE :

  function Generate_Cells
              ( v,w : List_of_Permutations; mix : Vector;
                mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision )
              return Integer_Mixed_Subdivisions.Mixed_Subdivision is

    use Integer_Mixed_Subdivisions;

    tmp,res,res_last : Mixed_Subdivision;

  begin
    tmp := mixsub;
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        Permute_and_Append(v,w,mic,mix,res,res_last);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Generate_Cells;

  function Exists_Permutation
               ( v1,v2 : Standard_Integer_Vectors.Link_to_Vector )
               return boolean is

  -- DESCRIPTION :
  --   Returns true if there exists a permutation p: v1 = p*v1.

  begin
    if v1(v1'last) /= v2(v2'last)
     then return false;  -- they must have the same lifted component !
     else return Permutable(v1(v1'first..v1'last-1),v2(v2'first..v2'last-1));
    end if;
  end Exists_Permutation;

  function Permutable ( mic : Integer_Mixed_Subdivisions.Mixed_Cell;
                        mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision )
                      return boolean is

    use Integer_Mixed_Subdivisions;
    tmp : Mixed_Subdivision := mixsub;
    mic2 : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic2 := Head_Of(tmp);
      if Exists_Permutation(mic.nor,mic2.nor)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Permutable;

  function Generating_Cells
             ( mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision )
             return Integer_Mixed_Subdivisions.Mixed_Subdivision is

    use Integer_Mixed_Subdivisions;

    tmp,res,res_last : Mixed_Subdivision;
    mic : Mixed_Cell;
  
  begin
    tmp := mixsub;
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      if not Permutable(mic,res)
       then Append(res,res_last,mic);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Generating_Cells;

-- FOURTH TARGET ROUTINE :

  function Exists_Permutation
              ( v1,v2 : Standard_Floating_Vectors.Link_to_Vector )
              return boolean is

  -- DESCRIPTION :
  --   Returns true if there exists a permutation p: v1 = p*v1.

  begin
    if v1(v1'last) /= v2(v2'last)
     then return false;  -- they must have the same lifted component !
     else return Permutable(v1(v1'first..v1'last-1),v2(v2'first..v2'last-1));
    end if;
  end Exists_Permutation;

  function Permutable ( mic : Floating_Mixed_Subdivisions.Mixed_Cell;
                        mixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision )
                      return boolean is

    use Floating_Mixed_Subdivisions;

    tmp : Mixed_Subdivision := mixsub;
    mic2 : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic2 := Head_Of(tmp);
      if Exists_Permutation(mic.nor,mic2.nor)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Permutable;

  function Generating_Cells
             ( mixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision )
             return Floating_Mixed_Subdivisions.Mixed_Subdivision is

    use Floating_Mixed_Subdivisions;

    tmp,res,res_last : Mixed_Subdivision;
    mic : Mixed_Cell;
  
  begin
    tmp := mixsub;
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      if not Permutable(mic,res)
       then Append(res,res_last,mic);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Generating_Cells;

end Generating_Mixed_Cells;
