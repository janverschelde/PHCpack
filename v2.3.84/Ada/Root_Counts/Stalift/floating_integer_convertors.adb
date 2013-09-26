with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package body Floating_Integer_Convertors is

  function Convert ( v : Standard_Integer_Vectors.Vector )
                   return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(v'range);

  begin
    for i in res'range loop
      res(i) := double_float(v(i));
    end loop;
    return res;
  end Convert;

  function Convert ( v : Standard_Floating_Vectors.Vector )
                   return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(v'range);

  begin
    for i in res'range loop
      res(i) := integer32(v(i));
    end loop;
    return res;
  end Convert;

  function Convert ( L : Lists_of_Integer_Vectors.List )
                   return Lists_of_Floating_Vectors.List is

    res,res_last : Lists_of_Floating_Vectors.List;
    tmp : Lists_of_Integer_Vectors.List := L;

    use Lists_of_Integer_Vectors;

  begin
    while not Is_Null(tmp) loop
      Lists_of_Floating_Vectors.Append(res,res_last,Convert(Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

  function Convert ( L : Lists_of_Floating_Vectors.List )
                   return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    tmp : Lists_of_Floating_Vectors.List := L;

    use Lists_of_Floating_Vectors;

  begin
    while not Is_Null(tmp) loop
      Lists_of_Integer_Vectors.Append(res,res_last,Convert(Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

  function Convert ( L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                   return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Floating_Vector_Lists.Array_of_Lists(L'range);

  begin
    for i in l'range loop
      res(i) := Convert(L(i));
    end loop;
    return res;
  end Convert;

  function Convert ( L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                   return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(L'range);

  begin
    for i in L'range loop
      res(i) := Convert(L(i));
    end loop;
    return res;
  end Convert;

  function Convert ( m : Integer_Mixed_Subdivisions.Mixed_Cell )
                   return Floating_Mixed_Subdivisions.Mixed_Cell is

    res : Floating_Mixed_Subdivisions.Mixed_Cell;
    use Integer_Mixed_Subdivisions;

  begin
    res.nor := new Standard_Floating_Vectors.Vector'(Convert(m.nor.all));
    res.pts := 
      new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(Convert(m.pts.all));
    if m.sub = null
     then res.sub := null;
     else res.sub := 
      new Floating_Mixed_Subdivisions.Mixed_Subdivision'(Convert(m.sub.all));
    end if;
    return res;
  end Convert;

  function Convert ( m : Floating_Mixed_Subdivisions.Mixed_Cell )
                   return Integer_Mixed_Subdivisions.Mixed_Cell is

    res : Integer_Mixed_Subdivisions.Mixed_Cell;
    use Floating_Mixed_Subdivisions;

  begin
    res.nor := new Standard_Integer_Vectors.Vector'(Convert(m.nor.all));
    res.pts :=
      new Arrays_of_Integer_Vector_Lists.Array_of_Lists'(Convert(m.pts.all));
    if m.sub = null
     then res.sub := null;
     else res.sub := 
       new Integer_Mixed_Subdivisions.Mixed_Subdivision'(Convert(m.sub.all));
    end if;
    return res;
  end Convert;

  function Convert ( s : Integer_Mixed_Subdivisions.Mixed_Subdivision )
                   return Floating_Mixed_Subdivisions.Mixed_Subdivision is

    res,res_last : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    tmp : Integer_Mixed_Subdivisions.Mixed_Subdivision := s;

    use Integer_Mixed_Subdivisions;

  begin
    while not Is_Null(tmp) loop
      Floating_Mixed_Subdivisions.Append(res,res_last,Convert(Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

  function Convert ( s : Floating_Mixed_Subdivisions.Mixed_Subdivision )
                   return Integer_Mixed_Subdivisions.Mixed_Subdivision is

    res,res_last : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    tmp : Floating_Mixed_Subdivisions.Mixed_Subdivision := s;

    use Floating_Mixed_Subdivisions;

  begin
    while not Is_Null(tmp) loop
      Integer_Mixed_Subdivisions.Append(res,res_last,Convert(Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

end Floating_Integer_Convertors;
