with unchecked_deallocation;

package body Generic_VecMats is

  procedure Clear ( v : in out VecMat ) is
  begin
    for i in v'range loop
      Matrices.Clear(v(i));
    end loop;
  end Clear;

  procedure Shallow_Clear ( v : in out Link_to_VecMat ) is

    procedure free is new unchecked_deallocation(VecMat,Link_to_VecMat);

  begin
    free(v);
  end Shallow_Clear;

  procedure Deep_Clear ( v : in out Link_to_VecMat ) is
  begin
    if v /= null
     then Clear(v.all); Shallow_Clear(v);
    end if;
  end Deep_Clear;

end Generic_VecMats;
