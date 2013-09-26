with Permutations;                       use Permutations;
with Permutations_of_Faces;              use Permutations_of_Faces;

package body Faces_of_Symmetric_Polytopes is

-- ON FACES : group * faces -> invariant subgroup

  function Stabilizer ( v : List_of_Permutations; f : Face ) 
                      return List_of_Permutations is

    tmp,res,res_last : List_of_Permutations;

  begin
    tmp := v;
    while not Is_Null(tmp) loop
      declare
        p : constant Permutation := Permutation(Head_Of(tmp).all);
        pf : Face := Permute(f,p);
      begin
        if Is_Equal(f,pf)
         then Append(res,res_last,p);
        end if;
        Deep_Clear(pf);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Stabilizer;

  function Stabilizer_Lifted ( v : List_of_Permutations; f : Face ) 
                             return List_of_Permutations is

    tmp,res,res_last : List_of_Permutations;

  begin
    tmp := v;
    while not Is_Null(tmp) loop
      declare
        p : constant Permutation := Permutation(Head_Of(tmp).all);
        pf : Face := Permute_Lifted(f,p);
      begin
        if Is_Equal(f,pf)
         then Append(res,res_last,p);
        end if;
        Deep_Clear(pf);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Stabilizer_Lifted;

-- ON FACES : group * faces -> invariant faces

  function Invariant_Faces ( v : List_of_Permutations;
                             f : Faces ) return Faces is

    tmpf,res,res_last : Faces;
    tmpv : List_of_Permutations;
    inva : boolean;

  begin
    tmpf := f;
    while not Is_Null(tmpf) loop
      inva := false;
      declare
        ff : constant Face := Head_Of(tmpf);
        cf : Face;
      begin
        tmpv := v;
        while not Is_Null(tmpv) loop
          inva := Invariant(ff,Permutation(Head_Of(tmpv).all));
          exit when not inva;
          tmpv := Tail_Of(tmpv);
        end loop;
        if inva
         then Copy(ff,cf); Append(res,res_last,cf);
        end if;
      end;
      tmpf := Tail_Of(tmpf);
    end loop;
    return res;
  end Invariant_Faces;

  function Invariant_Lifted_Faces ( v : List_of_Permutations;
                                    f : Faces ) return Faces is

    tmpf,res,res_last : Faces;
    tmpv : List_of_Permutations;
    inva : boolean;

  begin
    tmpf := f;
    while not Is_Null(tmpf) loop
      inva := false;
      declare
        ff : constant Face := Head_Of(tmpf);
        cf : Face;
      begin
        tmpv := v;
        while not Is_Null(tmpv) loop
          inva := Invariant_Lifted(ff,Permutation(Head_Of(tmpv).all));
          exit when not inva;
          tmpv := Tail_Of(tmpv);
        end loop;
        if inva
         then Copy(ff,cf); Append(res,res_last,cf);
        end if;
      end;
      tmpf := Tail_Of(tmpf);
    end loop;
    return res;
  end Invariant_Lifted_Faces;

-- ON FACES : group * faces -> generated faces

  function Generated_Faces ( v : List_of_Permutations; f : Faces )
                           return Faces is

    res : Faces;

  begin
    return res;
  end Generated_Faces;

  function Generated_Lifted_Faces
                          ( v : List_of_Permutations; f : Faces )
                          return Faces is

    res : Faces;

  begin
    return res;
  end Generated_Lifted_Faces;

-- ON FACES : group * faces -> generating faces

  function Generating_Faces ( f : Faces ) return Faces is

    tmp,res,res_last : Faces;
    lf : Face;

  begin
    tmp := f;
    while not Is_Null(tmp) loop
      lf := Head_Of(tmp);
      if not Permutable(lf,res)
       then Append(res,res_last,lf);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Generating_Faces;

  function Generating_Lifted_Faces ( f : Faces ) return Faces is

    tmp,res,res_last : Faces;
    lf : Face;

  begin
    tmp := f;
    while not Is_Null(tmp) loop
      lf := Head_Of(tmp);
      if not Permutable_Lifted(lf,res)
       then Append(res,res_last,lf);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Generating_Lifted_Faces;

  function Generating_Faces ( v : List_of_Permutations; f : Faces )
                            return Faces is

    tmp,res,res_last : Faces;
    tmpv : List_of_Permutations;
    found : boolean;
    lf : Face;

  begin
    tmp := f;
    while not Is_Null(tmp) loop
       lf := Head_Of(tmp);
       tmpv := v;
       while not Is_Null(tmpv) loop
          declare
            pv : constant Permutation := Permutation(Head_Of(tmpv).all);
          begin
            found := Is_In(res,Permute(lf,pv));
          end;
          exit when found;
          tmpv := Tail_Of(tmpv);
        end loop;
        if not found
         then Append(res,res_last,lf);
        end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Generating_Faces;

  function Generating_Lifted_Faces ( v : List_of_Permutations; f : Faces )
                                   return Faces is

    tmp,res,res_last : Faces;
    tmpv : List_of_Permutations;
    found : boolean;
    lf : Face;

  begin
    tmp := f;
    while not Is_Null(tmp) loop
       lf := Head_Of(tmp);
       tmpv := v;
       while not Is_Null(tmpv) loop
          declare
            pv : constant Permutation := Permutation(Head_Of(tmpv).all);
          begin
            found := Is_In(res,Permute_Lifted(lf,pv));
          end;
          exit when found;
          tmpv := Tail_Of(tmpv);
        end loop;
        if not found
         then Append(res,res_last,lf);
        end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Generating_Lifted_Faces;

-- ON TUPLES OF FACES : group * faces -> invariant faces

  function Invariant_Faces ( v : List_of_Permutations;
                             af : Array_of_Faces ) return Array_of_Faces is

    res : Array_of_Faces(af'range);

  begin
    for i in res'range loop
      res(i) := Invariant_Faces(v,af(i));
    end loop;
    return res;
  end Invariant_Faces;

  function Invariant_Lifted_Faces
                           ( v : List_of_Permutations; af : Array_of_Faces )
                           return Array_of_Faces is

    res : Array_of_Faces(af'range);

  begin
    for i in res'range loop
      res(i) := Invariant_Lifted_Faces(v,af(i));
    end loop;
    return res;
  end Invariant_Lifted_Faces;

-- ON TUPLES OF FACES : group * faces -> generators of faces

  function Generating_Faces ( af : Array_of_Faces ) return Array_of_Faces is

    res,res_last : Array_of_Faces(af'range);
    tmp : Faces;
    lf : Face;
    found : boolean;

  begin
    for i in af'range loop
      tmp := af(i);
      while not Is_Null(tmp) loop
        lf := Head_Of(tmp);
        for j in res'range loop
          found := Permutable(lf,res(j));
          exit when found;
        end loop;
        if not found
         then Append(res(i),res_last(i),lf);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Generating_Faces;

  function Generating_Lifted_Faces
             ( af : Array_of_Faces ) return Array_of_Faces is

    res,res_last : Array_of_Faces(af'range);
    tmp : Faces;
    lf : Face;
    found : boolean;

  begin
    for i in af'range loop
      tmp := af(i);
      while not Is_Null(tmp) loop
        lf := Head_Of(tmp);
        for j in res'range loop
          found := Permutable_Lifted(lf,res(j));
          exit when found;
        end loop;
        if not found
         then Append(res(i),res_last(i),lf);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Generating_Lifted_Faces;

  function Generating_Faces 
                 ( v : List_of_Permutations; af : Array_of_Faces )
                 return Array_of_Faces is

    res : Array_of_Faces(af'range);

  begin
    for i in res'range loop
      res(i) := Generating_Faces(v,af(i));
    end loop;
    return res;
  end Generating_Faces;

  function Generating_Lifted_Faces 
                 ( v : List_of_Permutations; af : Array_of_Faces )
                 return Array_of_Faces is

    res : Array_of_Faces(af'range);

  begin
    for i in res'range loop
      res(i) := Generating_Lifted_Faces(v,af(i));
    end loop;
    return res;
  end Generating_Lifted_Faces;

  function Generating_Faces
                 ( v,w : List_of_Permutations; af : Array_of_Faces )
                 return Array_of_Faces is

    res,res_last : Array_of_Faces(af'range);
    tmp : Faces;
    lf : Face;
    found : boolean;
    tmpv,tmpw : List_of_Permutations;

  begin
    for i in af'range loop
      tmp := af(i);
      while not Is_Null(tmp) loop
        lf := Head_Of(tmp);
        tmpv := v; tmpw := w;
        while not Is_Null(tmpv) loop
          declare
            pv : constant Permutation := Permutation(Head_Of(tmpv).all);
            pw : constant Permutation := Permutation(Head_Of(tmpw).all);
          begin
            found := Is_In(res(pw(i)),Permute(lf,pv));
          end;
          exit when found;
          tmpv := Tail_Of(tmpv); tmpw := Tail_Of(tmpw);
        end loop;
        if not found
         then Append(res(i),res_last(i),lf);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Generating_Faces;

  function Generating_Lifted_Faces
                 ( v,w : List_of_Permutations; af : Array_of_Faces )
                 return Array_of_Faces is

    res,res_last : Array_of_Faces(af'range);
    tmp : Faces;
    lf : Face;
    found : boolean;
    tmpv,tmpw : List_of_Permutations;

  begin
    for i in af'range loop
      tmp := af(i);
      while not Is_Null(tmp) loop
        lf := Head_Of(tmp);
        tmpv := v; tmpw := w;
        while not Is_Null(tmpv) loop
          declare
            pv : constant Permutation := Permutation(Head_Of(tmpv).all);
            pw : constant Permutation := Permutation(Head_Of(tmpw).all);
          begin
            found := Is_In(res(pw(i)),Permute_Lifted(lf,pv));
          end;
          exit when found;
          tmpv := Tail_Of(tmpv); tmpw := Tail_Of(tmpw);
        end loop;
        if not found
         then Append(res(i),res_last(i),lf);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Generating_Lifted_Faces;

end Faces_of_Symmetric_Polytopes;
