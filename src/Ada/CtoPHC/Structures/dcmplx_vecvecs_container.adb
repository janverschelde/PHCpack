package body DCMPLX_VecVecs_Container is

  data : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;

  procedure Initialize ( n : in integer32 ) is
  begin
    data := new Standard_Complex_VecVecs.Array_of_VecVecs(1..n);
  end Initialize;

  procedure Initialize ( m : in Standard_Integer_Vectors.Vector ) is

    use Standard_Complex_VecVecs;

  begin
    if data = null
     then data := new Standard_Complex_VecVecs.Array_of_VecVecs(m'range);
    end if;
    for k in m'range loop
      data(k) := new Standard_Complex_VecVecs.VecVec(1..m(k));
    end loop;
  end Initialize;

  procedure Store_Link
              ( v : in Standard_Complex_VecVecs.Link_to_Array_of_VecVecs ) is
  begin
    data := v;
  end Store_Link;

  procedure Store_Copy
              ( v : in Standard_Complex_VecVecs.Array_of_VecVecs ) is

    cpv : constant Standard_Complex_VecVecs.Array_of_VecVecs(v'range)
        := Standard_Complex_VecVecs.Create_Copy(v);

  begin
    data := new Standard_Complex_VecVecs.Array_of_VecVecs'(cpv);
  end Store_Copy;

  procedure Store_Link ( k : in integer32;
                         v : in Standard_Complex_VecVecs.Link_to_VecVec ) is
  begin
    data(k) := v;
  end Store_Link;

  procedure Store_Copy ( k : in integer32;
                         v : in Standard_Complex_VecVecs.VecVec ) is

    cpv : constant Standard_Complex_VecVecs.VecVec(v'range )
        := Standard_Complex_VecVecs.Create_Copy(v);

  begin
    data(k) := new Standard_Complex_VecVecs.VecVec'(cpv);
  end Store_Copy;

  procedure Store_Link ( k,i : in integer32;
                         v : in Standard_Complex_Vectors.Link_to_Vector ) is
  begin
    data(k)(i) := v;
  end Store_Link;

  procedure Store_Copy ( k,i : in integer32;
                         v : in Standard_Complex_Vectors.Vector ) is
  begin
    data(k)(i) := new Standard_Complex_Vectors.Vector'(v);
  end Store_Copy;

  function Size return integer32 is

    use Standard_Complex_VecVecs;

  begin
    if data = null
     then return 0;
     else return data'last;
    end if;
  end Size;

  function Size ( k : integer32 ) return integer32 is

    use Standard_Complex_VecVecs;

  begin
    if data = null then
      return 0;
    else
      if k > 0 and k <= data'last
       then return data(k)'last;
       else return 0;
      end if;
    end if;
  end Size;

  function Get return Standard_Complex_VecVecs.Link_to_Array_of_VecVecs is
  begin
    return data;
  end Get;

  function Get ( k : integer32 )
               return Standard_Complex_VecVecs.Link_to_VecVec is

    res : Standard_Complex_VecVecs.Link_to_VecVec := null;

    use Standard_Complex_VecVecs;

  begin
    if data /= null then
      if k > 0 and k <= data'last
       then res := data(k);
      end if;
    end if;
    return res;
  end Get;

  function Get ( k,i : integer32 )
               return Standard_Complex_Vectors.Link_to_Vector is

    res : Standard_Complex_Vectors.Link_to_Vector;
    dk : Standard_Complex_VecVecs.Link_to_VecVec;

    use Standard_Complex_VecVecs;

  begin
    if data /= null then
      if k > 0 and k <= data'last
       then dk := data(k);
      end if;
    end if;
    if dk /= null then
      if i > 0 and i <= dk'last
       then res := dk(i);
      end if;
    end if;
    return res;
  end Get;

  procedure Clear is

    use Standard_Complex_VecVecs;

  begin
    if data /= null
     then Deep_Clear(data);
    end if;
  end Clear;

begin
  data := null;
end DCMPLX_VecVecs_Container;
