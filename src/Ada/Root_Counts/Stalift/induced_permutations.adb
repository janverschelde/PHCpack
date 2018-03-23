with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Vectors;

package body Induced_Permutations is

  function Remove_Artificial_Origin
             ( L : Lists_of_Floating_Vectors.List;
               b : double_float )
             return Lists_of_Floating_Vectors.List is

    use Lists_of_Floating_Vectors;
    res,res_last : List;
    lv : Standard_Floating_Vectors.Link_to_Vector;
    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      if lv(lv'last) /= b
       then Append(res,res_last,lv.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Artificial_Origin;

  procedure Remove_Artificial_Origin
              ( L : in out Lists_of_Floating_Vectors.List;
                b : in double_float ) is

    use Lists_of_Floating_Vectors;
    res : List := Remove_Artificial_Origin(L,b);

  begin
    Deep_Clear(L);
    L := res;
  end Remove_Artificial_Origin;

  function Remove_Artificial_Origin
              ( L : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                b : double_float )
              return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Floating_Vector_Lists.Array_of_Lists(L'range);

  begin
    for i in L'range loop
      res(i) := Remove_Artificial_Origin(L(i),b);
    end loop;
    return res;
  end Remove_Artificial_Origin;

  procedure Remove_Artificial_Origin
              ( L : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                b : in double_float ) is
  begin
    for i in L'range loop
      Remove_Artificial_Origin(L(i),b);
    end loop;
  end Remove_Artificial_Origin;

  function Is_Subset
             ( lifted,original : Lists_of_Floating_Vectors.List )
             return boolean is

    use Lists_of_Floating_Vectors;
    tmp : List := lifted;
    lv : Standard_Floating_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      if not Is_In(original,lv(lv'first..lv'last-1))
       then return false;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return true;
  end Is_Subset;

  function Permutation
             ( s,ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               mix : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(s'range) := (s'range => 0);
    offset : integer32 := 0;
    min_index : integer32;
    found : boolean;
    use Lists_of_Floating_Vectors;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        if j > 1
         then offset := offset + 1;
        end if;
        for k in s'range loop
          found := false;
          if res(k) = 0 then
            if Is_Subset(ls(i),s(k)) then
              res(k) := i + offset; found := true;
              min_index := k;
              if found then
                for kk in k+1..s'last loop
                  if res(kk) = 0 then
                    if Is_Subset(ls(i),s(kk)) then
                      if Length_Of(s(kk)) < Length_Of(s(min_index)) then
                        res(min_index) := 0;
                        min_index := kk;
                        res(min_index) := i + offset;
                      end if;
                    end if;
                  end if;
                end loop;
              end if;
            end if;
          end if;
          exit when found;
        end loop;
      end loop;
    end loop;
    return res;
  end Permutation;

  function Shift_Indices 
             ( p : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is
  begin
    if p'first = 1 then
      return p;
    else
      declare
        res : Standard_Integer_Vectors.Vector(1..p'last+1);
      begin
        for i in p'range loop
          res(i+1) := p(i);
        end loop;
        return res;
      end;
    end if;
  end Shift_Indices;

  function Relabel_for_Zero
             ( p : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

     inzero : boolean := false;

  begin
    for i in p'range loop
      if p(i) = 0
       then inzero := true; exit;
      end if;
    end loop;
    if not inzero then
      return Shift_Indices(p);
    else
      declare
        rlp : Standard_Integer_Vectors.Vector(p'range);
      begin
        for i in p'range loop
          rlp(i) := p(i) + 1;
        end loop;
        return Shift_Indices(rlp);
      end;
    end if;
  end Relabel_for_Zero;
 
  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

    pf : Standard_Complex_Poly_Systems.Poly_Sys(f'range);
    rp : constant Standard_Integer_Vectors.Vector := Relabel_for_Zero(p);

  begin
    for i in p'range loop
      pf(rp(i)) := f(i);
    end loop;
    f := pf;
  end Permute;

  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out Standard_Complex_Laur_Systems.Laur_Sys ) is

    pf : Standard_Complex_Laur_Systems.Laur_Sys(f'range);
    rp : constant Standard_Integer_Vectors.Vector := Relabel_for_Zero(p);
  
  begin
    for i in p'range loop
      pf(rp(i)) := f(i);
    end loop;
    f := pf;
  end Permute;

  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    pf : DoblDobl_Complex_Poly_Systems.Poly_Sys(f'range);
    rp : constant Standard_Integer_Vectors.Vector := Relabel_for_Zero(p);

  begin
    for i in pf'range loop
      pf(rp(i)) := f(i);
    end loop;
    f := pf;
  end Permute;

  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    pf : DoblDobl_Complex_Laur_Systems.Laur_Sys(f'range);
    rp : constant Standard_Integer_Vectors.Vector := Relabel_for_Zero(p);

  begin
    for i in p'range loop
      pf(rp(i)) := f(i);
    end loop;
    f := pf;
  end Permute;

  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    pf : QuadDobl_Complex_Poly_Systems.Poly_Sys(f'range);
    rp : constant Standard_Integer_Vectors.Vector := Relabel_for_Zero(p);

  begin
    for i in p'range loop
      pf(rp(i)) := f(i);
    end loop;
    f := pf;
  end Permute;

  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    pf : QuadDobl_Complex_Laur_Systems.Laur_Sys(f'range);
    rp : constant Standard_Integer_Vectors.Vector := Relabel_for_Zero(p);

  begin
    for i in p'range loop
      pf(rp(i)) := f(i);
    end loop;
    f := pf;
  end Permute;

end Induced_Permutations;
