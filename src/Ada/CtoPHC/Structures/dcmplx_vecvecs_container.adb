with Ada.Text_IO;                       use Ada.Text_IO;

package body DCMPLX_VecVecs_Container is

  data : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;

  procedure Initialize ( n : in integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Initialize 0 ...");
    end if;
    data := new Standard_Complex_VecVecs.Array_of_VecVecs(1..n);
  end Initialize;

  procedure Initialize ( m : in Standard_Integer_Vectors.Vector;
                         vrblvl : in integer32 := 0 ) is

    use Standard_Complex_VecVecs;

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Initialize 1 ...");
    end if;
    if data = null
     then data := new Standard_Complex_VecVecs.Array_of_VecVecs(m'range);
    end if;
    for k in m'range loop
      data(k) := new Standard_Complex_VecVecs.VecVec(1..m(k));
    end loop;
  end Initialize;

  procedure Store_Link
              ( v : in Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Store_Link 0 ...");
    end if;
    data := v;
  end Store_Link;

  procedure Store_Copy
              ( v : in Standard_Complex_VecVecs.Array_of_VecVecs;
                vrblvl : in integer32 := 0 ) is

    cpv : constant Standard_Complex_VecVecs.Array_of_VecVecs(v'range)
        := Standard_Complex_VecVecs.Create_Copy(v);

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Store_Copy 0 ...");
    end if;
    data := new Standard_Complex_VecVecs.Array_of_VecVecs'(cpv);
  end Store_Copy;

  procedure Store_Link ( k : in integer32;
                         v : in Standard_Complex_VecVecs.Link_to_VecVec;
                         vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Store_Link 1 ...");
    end if;
    data(k) := v;
  end Store_Link;

  procedure Store_Copy ( k : in integer32;
                         v : in Standard_Complex_VecVecs.VecVec;
                         vrblvl : in integer32 := 0 ) is

    cpv : constant Standard_Complex_VecVecs.VecVec(v'range )
        := Standard_Complex_VecVecs.Create_Copy(v);

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Store_Copy 1 ...");
    end if;
    data(k) := new Standard_Complex_VecVecs.VecVec'(cpv);
  end Store_Copy;

  procedure Store_Link ( k,i : in integer32;
                         v : in Standard_Complex_Vectors.Link_to_Vector;
                         vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Store_Link 2 ...");
    end if;
    data(k)(i) := v;
  end Store_Link;

  procedure Store_Copy ( k,i : in integer32;
                         v : in Standard_Complex_Vectors.Vector;
                         vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Store_Copy 2 ...");
    end if;
    data(k)(i) := new Standard_Complex_Vectors.Vector'(v);
  end Store_Copy;

  function Size ( vrblvl : in integer32 := 0 ) return integer32 is

    use Standard_Complex_VecVecs;

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Size 0 ...");
    end if;
    if data = null
     then return 0;
     else return data'last;
    end if;
  end Size;

  function Size ( k : integer32; vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_VecVecs;

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Size 1 ...");
    end if;
    if data = null then
      return 0;
    else
      if k > 0 and k <= data'last
       then return data(k)'last;
       else return 0;
      end if;
    end if;
  end Size;

  function Size ( k,i : integer32;
                  vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_VecVecs;

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Size 2 ...");
    end if;
    if data = null then
      return 0;
    else
      if k < 1 or k > data'last then
        return 0;
      else
        declare
          vk : constant Link_to_VecVec := data(k);
        begin
          if i < vk'first or i > vk'last
           then return 0;
           else return vk(i)'last;
          end if;
        end;
      end if;
    end if;
  end Size;

  function Get ( vrblvl : integer32 := 0 )
               return Standard_Complex_VecVecs.Link_to_Array_of_VecVecs is
  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Get 0 ...");
    end if;
    return data;
  end Get;

  function Get ( k : integer32; vrblvl : integer32 := 0 )
               return Standard_Complex_VecVecs.Link_to_VecVec is

    res : Standard_Complex_VecVecs.Link_to_VecVec := null;

    use Standard_Complex_VecVecs;

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Get 1 ...");
    end if;
    if data /= null then
      if k > 0 and k <= data'last
       then res := data(k);
      end if;
    end if;
    return res;
  end Get;

  function Get ( k,i : integer32; vrblvl : integer32 := 0 )
               return Standard_Complex_Vectors.Link_to_Vector is

    res : Standard_Complex_Vectors.Link_to_Vector;
    dk : Standard_Complex_VecVecs.Link_to_VecVec;

    use Standard_Complex_VecVecs;

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Get 2 ...");
    end if;
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

  procedure Clear ( vrblvl : in integer32 := 0 ) is

    use Standard_Complex_VecVecs;

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_container.Clear ...");
    end if;
    if data /= null
     then Deep_Clear(data);
    end if;
  end Clear;

begin
  data := null;
end DCMPLX_VecVecs_Container;
