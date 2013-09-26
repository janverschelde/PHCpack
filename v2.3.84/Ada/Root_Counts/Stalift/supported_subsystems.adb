with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;

package body Supported_Subsystems is

  function Is_Equal ( v1 : Standard_Floating_Vectors.Vector;
                      v2 : Standard_Natural_Vectors.Vector )
                    return boolean is
  begin
    for i in v2'range loop
      if double_float(v2(i)) /= v1(i)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Equal;

  function Is_Equal ( v1 : Standard_Floating_Vectors.Vector;
                      v2 : Standard_Integer_Vectors.Vector )
                    return boolean is
  begin
    for i in v2'range loop
      if double_float(v2(i)) /= v1(i)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Equal;

  function Is_In ( s : Lists_of_Floating_Vectors.List;
                   v : Standard_Natural_Vectors.Vector ) return boolean is

    use Lists_of_Floating_Vectors;
    tmp : List := s;
    l2v : Standard_Floating_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      if Is_Equal(l2v.all,v)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( s : Lists_of_Floating_Vectors.List;
                   v : Standard_Integer_Vectors.Vector ) return boolean is

    use Lists_of_Floating_Vectors;
    tmp : List := s;
    l2v : Standard_Floating_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      if Is_Equal(l2v.all,v)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Select_Terms ( p : Standard_Complex_Polynomials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      v : constant Standard_Natural_Vectors.Vector(t.dg'range) := t.dg.all;

    begin
      if Is_In(s,v)
       then Add(res,t);
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Select_Terms;

  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Poly_Systems;
    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector; 
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Poly_Systems;
    res : Poly_Sys(p'range);
    ind : integer32 := 0;

  begin
    for i in m'range loop
      for j in 1..m(i) loop
        ind := ind + 1;
        res(ind) := Select_Terms(p(ind),s(i));
      end loop;
    end loop;
    return res;
  end Select_Terms;

end Supported_Subsystems;
