with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;

package body Multprec_Floating_Vector_Tools is

  function Create ( v : Standard_Floating_Vectors.Vector )
                  return Multprec_Floating_Vectors.Vector is

    res : Multprec_Floating_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Create(v(i));
    end loop;
    return res;
  end Create;

  function Create ( v : Standard_Floating_VecVecs.VecVec )
                  return Multprec_Floating_VecVecs.VecVec is

    res : Multprec_Floating_VecVecs.VecVec(v'range);
    use Standard_Floating_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then res(i) := new Multprec_Floating_Vectors.Vector'(Create(v(i).all));
      end if;
    end loop;
    return res;
  end Create;

  function Create ( v : Standard_Floating_VecVecs.Array_of_VecVecs )
                  return Multprec_Floating_VecVecs.Array_of_VecVecs is

    res : Multprec_Floating_VecVecs.Array_of_VecVecs(v'range);
    use Standard_Floating_VecVecs;

  begin
    for i in v'range loop
      if v(i) /= null
       then res(i) := new Multprec_Floating_VecVecs.VecVec'(Create(v(i).all));
      end if;
    end loop;
    return res;
  end Create;

  procedure Set_Size ( v : in out Multprec_Floating_Vectors.Vector;
                       size : in natural32 ) is
  begin
    for i in v'range loop
      Set_Size(v(i),size);
    end loop;
  end Set_Size;

  procedure Set_Size ( v : in out Multprec_Floating_VecVecs.VecVec;
                       size : in natural32 ) is

    use Multprec_Floating_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Set_Size(v(i).all,size);
      end if;
    end loop;
  end Set_Size;

  procedure Set_Size ( v : in out Multprec_Floating_VecVecs.Array_of_VecVecs;
                       size : in natural32 ) is

    use Multprec_Floating_VecVecs;

  begin
    for i in v'range loop
      if v(i) /= null
       then Set_Size(v(i).all,size);
      end if;
    end loop;
  end Set_Size;

  function Round ( v : Multprec_Floating_Vectors.Vector )
                 return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Round(v(i));
    end loop;
    return res;
  end Round;

end Multprec_Floating_Vector_Tools;
