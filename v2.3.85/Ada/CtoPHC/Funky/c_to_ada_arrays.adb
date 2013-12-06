with Interfaces.C;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;

package body C_to_Ada_Arrays is

  function Convert ( v : C_Integer_Array )
                   return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(0..integer32(v'length)-1);

  begin
    for i in v'range loop
      res(integer32(i)) := integer32(v(i));
    end loop;
    return res;
  end Convert;

  function Convert ( v : Standard_Integer_Vectors.Vector )
                   return C_Integer_Array is

    res : C_Integer_Array(0..Interfaces.C.size_T(v'length-1));
    ind : Interfaces.C.size_T := 0;

    use Interfaces.C;

  begin
    for i in v'range loop
      res(ind) := Interfaces.C.int(v(i));
      ind := ind + 1;
    end loop;
    return res;
  end Convert;

  function Convert ( v : C_Double_Array )
                   return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(0..integer32(v'length)-1);

  begin
    for i in v'range loop
      res(integer32(i)) := double_float(v(i));
    end loop;
    return res;
  end Convert;

  function Convert ( v : Standard_Floating_Vectors.Vector )
                   return C_Double_Array is

    res : C_Double_Array(0..Interfaces.C.size_T(v'length-1));
    ind : Interfaces.C.size_T := 0;

    use Interfaces.C;

  begin
    for i in v'range loop
      res(ind) := Interfaces.C.double(v(i));
      ind := ind + 1;
    end loop;
    return res;
  end Convert;

  function Convert ( v : C_Double_Array )
                   return Standard_Complex_Vectors.Vector is

    vlen : constant integer32 := integer32(v'length);
    res : Standard_Complex_Vectors.Vector(1..vlen/2);
    ind : Interfaces.C.size_T := 0;

    use Interfaces.C;

  begin
    for i in res'range loop
      res(i) := Create(double_float(v(ind)),double_float(v(ind+1)));
      ind := ind + 2;
    end loop;
    return res;
  end Convert;

  function Convert ( v : Standard_Complex_Vectors.Vector )
                   return C_Double_Array is

    res : C_Double_Array(0..Interfaces.C.size_T(2*v'length-1));
    ind : Interfaces.C.size_T := 0;

    use Interfaces.C;

  begin
    for i in v'range loop
      res(ind) := Interfaces.C.double(REAL_PART(v(i)));
      ind := ind + 1;
      res(ind) := Interfaces.C.double(IMAG_PART(v(i)));
      ind := ind + 1;
    end loop;
    return res;
  end Convert;

  function Convert ( m : Standard_Complex_Matrices.Matrix )
                   return C_Double_Array is

    len : constant integer32 := 2*m'length(1)*m'length(2);
    res : C_Double_Array(0..Interfaces.C.size_T(len-1));
    ind : Interfaces.C.size_T := 0;

    use Interfaces.C;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        res(ind) := Interfaces.C.double(REAL_PART(m(i,j)));
        ind := ind + 1;
        res(ind) := Interfaces.C.double(IMAG_PART(m(i,j)));
        ind := ind + 1;
      end loop;
    end loop;
    return res;
  end Convert;

end C_to_Ada_Arrays;
