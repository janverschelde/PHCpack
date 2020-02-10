with unchecked_deallocation;

package body Generic_VecVecs is

  function Create_Copy ( v : VecVec ) return VecVec is

    res : VecVec(v'range);

  begin
    for i in v'range loop
      if v(i) /= null
       then res(i) := new Vector'(v(i).all);
      end if;
    end loop;
    return res;
  end Create_Copy;

  function Create_Copy ( v : Array_of_VecVecs ) return Array_of_VecVecs is

    res : Array_of_VecVecs(v'range);

  begin
    for i in v'range loop
      if v(i) /= null
       then res(i) := new VecVec'(Create_Copy(v(i).all));
      end if;
    end loop;
    return res;
  end Create_Copy;

  procedure Copy ( v : in VecVec; w : in out VecVec ) is
  begin
    Clear(w);
    for i in v'range loop
      if v(i) /= null then
        declare
          vec : constant Vector := v(i).all;
        begin
          w(i) := new Vector'(vec);
        end;
      end if;
    end loop;
  end Copy;

  procedure Clear ( v : in out VecVec ) is
  begin
    for i in v'range loop
      Vectors.Clear(v(i));
    end loop;
  end Clear;

  procedure Shallow_Clear ( v : in out Link_to_VecVec ) is

    procedure free is new unchecked_deallocation(VecVec,Link_to_VecVec);

  begin
    free(v);
  end Shallow_Clear;

  procedure Shallow_Clear ( v : in out Link_to_Array_of_VecVecs ) is

    procedure free is
      new unchecked_deallocation(Array_of_VecVecs,Link_to_Array_of_VecVecs);

  begin
    free(v);
  end Shallow_Clear;
          
  procedure Clear ( v : in out Array_of_VecVecs ) is
  begin
    for i in v'range loop
      Deep_Clear(v(i));
    end loop;
  end Clear;

  procedure Deep_Clear ( v : in out Link_to_VecVec ) is
  begin
    if v /= null
     then Clear(v.all);
          Shallow_Clear(v);
    end if;
  end Deep_Clear;

  procedure Deep_Clear ( v : in out Link_to_Array_of_VecVecs ) is
  begin
    if v /= null
     then Clear(v.all);
          Shallow_Clear(v);
    end if;
  end Deep_Clear;

end Generic_VecVecs;
