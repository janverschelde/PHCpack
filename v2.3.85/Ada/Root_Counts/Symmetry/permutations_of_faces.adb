with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;
with Permute_Operations;                 use Permute_Operations;

package body Permutations_of_Faces is

  function Invariant ( f : Face; p : Permutation ) return boolean is

    fp : Face := Permute(f,p);
    res : constant boolean := Is_Equal(f,fp);

  begin
    Deep_Clear(fp);
    return res;
  end Invariant;

  function Invariant_Lifted ( f : Face; p : Permutation ) return boolean is

    fp : Face := Permute_Lifted(f,p);
    res : constant boolean := Is_Equal(f,fp);

  begin
    Deep_Clear(fp);
    return res;
  end Invariant_Lifted;

  function Permute ( f : Face; p : Permutation ) return Face is

    res : constant Face := new VecVec(f'range);

  begin
    for i in f'range loop
      res(i) := new Vector'(p*f(i).all);
    end loop;
    return res;
  end Permute;

  function Permute_Lifted ( f : Face; p : Permutation ) return Face is

    res : constant Face := new VecVec(f'range);

  begin
    for i in f'range loop
      declare
        pt : constant Vector := f(i)(f(i)'first..f(i)'last-1);
      begin
        res(i) := new Vector(f(i)'range);
        res(i)(pt'range) := p*pt;
        res(i)(res(i)'last) := f(i)(f(i)'last);
      end;
    end loop;
    return res;
  end Permute_Lifted;

  function Permutable ( f1,f2 : Face ) return boolean is
    
    res : boolean;

  begin
    for i in f1'range loop
      res := false;
      for j in f2'range loop
        res := Permutable(f1(i).all,f2(j).all);
        exit when res;
      end loop;
      exit when not res;
    end loop;
    return res;
  end Permutable;

  function Permutable_Lifted ( f1,f2 : Face ) return boolean is
    
    res : boolean;

  begin
    for i in f1'range loop
      res := false;
      for j in f2'range loop
        if f1(i)(f1(i)'last) = f2(j)(f2(j)'last)  -- same lifting
         then res := Permutable(f1(i)(f1(i)'first..f1(i)'last-1),
                                f2(j)(f2(j)'first..f2(j)'last-1));
        end if;
        exit when res;
      end loop;
      exit when not res;
    end loop;
    return res;
  end Permutable_Lifted;

  function Permutable ( f1 : Face; f2 : Faces ) return boolean is

    tmp : Faces := f2;

  begin
    while not Is_Null(tmp) loop
      if Permutable(f1,Head_Of(tmp))
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Permutable;

  function Permutable_Lifted ( f1 : Face; f2 : Faces ) return boolean is

    tmp : Faces := f2;

  begin
    while not Is_Null(tmp) loop
      if Permutable_Lifted(f1,Head_Of(tmp))
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Permutable_Lifted;

end Permutations_of_Faces;
