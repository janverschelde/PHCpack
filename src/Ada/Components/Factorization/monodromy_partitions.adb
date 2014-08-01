with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors_io;         use Standard_Natural_Vectors_io;

package body Monodromy_Partitions is 

  function Init_Factors ( d : natural32 ) return Link_to_VecVec is

    res : constant Link_to_VecVec := new VecVec(1..integer32(d));

  begin
    for i in 1..d loop
      res(integer32(i)) := new Standard_Natural_Vectors.Vector'(1..1 => i);
    end loop;
    return res;
  end Init_Factors;

  function Map ( t1,t2 : Standard_Complex_Vectors.Vector; tol : double_float )
               return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(t1'range) := (t1'range => 0);
    use Standard_Complex_Numbers;

  begin
    for i in t1'range loop
      for j in t2'range loop
        if (abs(REAL_PART(t1(i)) - REAL_PART(t2(j))) < tol)
            and then (abs(IMAG_PART(t1(i)) - IMAG_PART(t2(j))) < tol)
         then res(i) := natural32(j); exit;
        end if;
      end loop;
    end loop;
    return res;
  end Map;

  function Map ( t1,t2 : DoblDobl_Complex_Vectors.Vector; tol : double_float )
               return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(t1'range) := (t1'range => 0);
    use Double_Double_Numbers,DoblDobl_Complex_Numbers;

  begin
    for i in t1'range loop
      for j in t2'range loop
        if (abs(REAL_PART(t1(i)) - REAL_PART(t2(j))) < tol)
            and then (abs(IMAG_PART(t1(i)) - IMAG_PART(t2(j))) < tol)
         then res(i) := natural32(j); exit;
        end if;
      end loop;
    end loop;
    return res;
  end Map;

  function Map ( t1,t2 : QuadDobl_Complex_Vectors.Vector; tol : double_float )
               return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(t1'range) := (t1'range => 0);
    use Quad_Double_Numbers,QuadDobl_Complex_Numbers;

  begin
    for i in t1'range loop
      for j in t2'range loop
        if (abs(REAL_PART(t1(i)) - REAL_PART(t2(j))) < tol)
            and then (abs(IMAG_PART(t1(i)) - IMAG_PART(t2(j))) < tol)
         then res(i) := natural32(j); exit;
        end if;
      end loop;
    end loop;
    return res;
  end Map;

  procedure Write_Map ( file : in file_type;
                        map : in Standard_Natural_Vectors.Vector ) is
  begin
    put_line(file,"The map of monodromy loops : ");
    for i in map'range loop
      put(file,i,3);
      put(file," -> ");
      put(file,map(i),1);
      new_line(file);
    end loop;
  end Write_Map;

  function Is_In ( v : Standard_Natural_Vectors.Vector;
                   i : natural32 ) return boolean is
  begin
    for k in v'range loop
      if v(k) = i
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_Connected ( deco : Link_to_VecVec; i,j : in natural32 ) 
                        return boolean is

    use Standard_Natural_Vectors;

  begin
    for k in deco'range loop
      if deco(k) /= null then
        if Is_In(deco(k).all,i) then
          if Is_In(deco(k).all,j)
           then return true;
           else return false;
          end if;
        elsif Is_In(deco(k).all,j) then
          return false;
        end if;
      end if;
    end loop;
    return false;
  end Is_Connected;

  function Merge_Sort ( v1,v2 : Standard_Natural_Vectors.Vector )
                      return Standard_Natural_Vectors.Vector is 

  -- DESCRIPTION :
  --   Given two vectors v1 and v2, sorted in increasing order,
  --   the vector on return contains all entries of v1 and v2,
  --   and is also sorted in increasing order.

    res : Standard_Natural_Vectors.Vector(v1'first..v1'last+v2'length);
    i1 : integer32 := v1'first;
    i2 : integer32 := v2'first;
    ind : integer32 := res'first-1;

  begin
    loop
      ind := ind+1;
      if v1(i1) < v2(i2)
       then res(ind) := v1(i1); i1 := i1 + 1;
       else res(ind) := v2(i2); i2 := i2 + 1;
      end if;
      exit when ((i1 > v1'last) or (i2 > v2'last));
    end loop;
    if i1 > v1'last then
      for i in i2..v2'last loop
        ind := ind + 1;
        res(ind) := v2(i);
      end loop;
    else
      for i in i1..v1'last loop
        ind := ind + 1;
        res(ind) := v1(i);
      end loop;
    end if;
    return res;
  end Merge_Sort;

  procedure Connect ( deco : in Link_to_VecVec;
                      i,j : in natural32 ) is

    use Standard_Natural_Vectors;

  begin
    for k1 in deco'range loop
      if deco(k1) /= null and then Is_In(deco(k1).all,i) then
        for k2 in k1+1..deco'last loop
          if deco(k2) /= null and then Is_In(deco(k2).all,j) then
            declare
              -- l : constant integer32 := deco(k1)'last+deco(k2)'last;
              -- v : Standard_Natural_Vectors.Vector(1..l);
              -- l1 : constant natural := deco(k1)'last;
              v : constant Standard_Natural_Vectors.Vector
                := Merge_Sort(deco(k1).all,deco(k2).all);
            begin
              -- v(deco(k1)'range) := deco(k1).all;
              -- for k3 in l1+1..l loop
              --   v(k3) := deco(k2)(k3-l1);
              -- end loop;
              Standard_Natural_Vectors.Clear(deco(k1));
              Standard_Natural_Vectors.Clear(deco(k2));
              deco(k1) := new Standard_Natural_Vectors.Vector'(v);
            end;
          end if;
        end loop;
      end if;
    end loop;
  end Connect;

  procedure Add_Map ( deco : in Link_to_VecVec; nb : in out natural32;
                      map : in Standard_Natural_Vectors.Vector ) is
  begin
    for i in map'range loop
      if (map(i) > 0) and then not Is_Connected(deco,natural32(i),map(i))
       then Connect(deco,natural32(i),map(i));
            nb := nb - 1;
      end if;
    end loop;
    if deco = null
     then nb := 0;
     else nb := Number_of_Factors(deco.all);
    end if;
  end Add_Map;

  function Number_of_Factors ( deco : in VecVec ) return natural32 is

    use Standard_Natural_Vectors;
    res : natural32 := 0;

  begin
    for i in deco'range loop
      if deco(i) /= null
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Factors;

  procedure Write_Factors ( file : in file_type; deco : in VecVec ) is

    use Standard_Natural_Vectors;

  begin
    put(file,Number_of_Factors(deco),1); new_line(file);
    for i in deco'range loop
      if deco(i) /= null then
        put(file,natural32(deco(i)'length),1); put(file," : ");
        put(file,deco(i).all); new_line(file);
      end if;
    end loop;
  end Write_Factors;

  procedure Write_Factors ( file : in file_type; deco : in VecVec;
                            m : in Standard_Natural_Vectors.Vector ) is

    use Standard_Natural_Vectors;

  begin
    put(file,Number_of_Factors(deco),1); new_line(file);
    for i in deco'range loop
      if deco(i) /= null then
        put(file,natural32(deco(i)'length),1); put(file," : ");
        put(file,deco(i).all); put(file," : m = ");
        put(file,m(integer32(deco(i)(deco(i)'first))),1); new_line(file);
      end if;
    end loop;
  end Write_Factors;

  procedure Assign_Legend ( deco : in Link_to_VecVec;
                            legend : in Standard_Natural_Vectors.Vector ) is

    use Standard_Natural_Vectors;
    d : integer32;

  begin
    if deco /= null then
      for i in deco'range loop
        if deco(i) /= null then
          for j in deco(i)'range loop
            d := integer32(deco(i)(j));
            if d >= legend'first and d <= legend'last
             then deco(i)(j) := legend(d);
            end if;
          end loop;
        end if;
      end loop;
    end if;
  end Assign_Legend;

  procedure Merge ( deco : in Link_to_VecVec;
                    subdeco : in Link_to_VecVec ) is

    use Standard_Natural_Vectors;
    ind1,ind2 : natural32;

  begin
    if subdeco /= null then
      for i in subdeco'range loop
        if subdeco(i) /= null then
          ind1 := subdeco(i)(subdeco(i)'first);
          for j in subdeco(i)'first+1..subdeco(i)'last loop
            ind2 := subdeco(i)(j);
            if not Is_Connected(deco,ind1,ind2)
             then Connect(deco,ind1,ind2);
            end if;
          end loop;
        end if;
      end loop;
    end if;
  end Merge;

  procedure Remove_Empty_Entries ( deco : in out Link_to_VecVec ) is

    use Standard_Natural_Vectors;

  begin
    if deco /= null then
      declare
        res : VecVec(deco'range);
        ind : integer32 := res'first-1;
      begin
        for i in deco'range loop
          if deco(i) /= null then
            ind := ind + 1;
            res(ind) := new Standard_Natural_Vectors.Vector'(deco(i).all);
          end if;
        end loop;
        Deep_Clear(deco);
        deco := new VecVec'(res(res'first..ind));
      end;
    end if;
  end Remove_Empty_Entries;

end Monodromy_Partitions;
