with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Floating_Integer_Convertors;
with Transforming_Integer32_Vector_Lists;

package body Supports_of_Polynomial_Systems is

  function Create ( p : Standard_Complex_Polynomials.Poly )
                  return Lists_of_Integer_Vectors.List is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    res,res_last : List;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           cont : out boolean ) is

      h : Link_to_Vector;

    begin
      h := new Standard_Integer_Vectors.Vector(t.dg'range);
      for j in h'range loop
        h(j) := integer32(t.dg(j));
      end loop;
      Append(res,res_last,h);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Create;

  function Create ( p : Standard_Complex_Laurentials.Poly )
                  return Lists_of_Integer_Vectors.List is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    res,res_last : List;

    procedure Visit_Term ( t : in Standard_Complex_Laurentials.Term;
                           cont : out boolean ) is

      h : Link_to_Vector;

    begin
      h := new Standard_Integer_Vectors.Vector(t.dg'range);
      for j in h'range loop
        h(j) := t.dg(j);
      end loop;
      Append(res,res_last,h);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Laurentials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Polynomials.Poly )
                  return Lists_of_Integer_Vectors.List is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    res,res_last : List;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           cont : out boolean ) is

      h : Link_to_Vector;

    begin
      h := new Standard_Integer_Vectors.Vector(t.dg'range);
      for j in h'range loop
        h(j) := integer32(t.dg(j));
      end loop;
      Append(res,res_last,h);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Laurentials.Poly )
                  return Lists_of_Integer_Vectors.List is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    res,res_last : List;

    procedure Visit_Term ( t : in DoblDobl_Complex_Laurentials.Term;
                           cont : out boolean ) is

      h : Link_to_Vector;

    begin
      h := new Standard_Integer_Vectors.Vector(t.dg'range);
      for j in h'range loop
        h(j) := t.dg(j);
      end loop;
      Append(res,res_last,h);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Laurentials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Create;


  function Create ( p : QuadDobl_Complex_Polynomials.Poly )
                  return Lists_of_Integer_Vectors.List is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    res,res_last : List;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           cont : out boolean ) is

      h : Link_to_Vector;

    begin
      h := new Standard_Integer_Vectors.Vector(t.dg'range);
      for j in h'range loop
        h(j) := integer32(t.dg(j));
      end loop;
      Append(res,res_last,h);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Create;

  function Create ( p : QuadDobl_Complex_Laurentials.Poly )
                  return Lists_of_Integer_Vectors.List is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    res,res_last : List;

    procedure Visit_Term ( t : in QuadDobl_Complex_Laurentials.Term;
                           cont : out boolean ) is

      h : Link_to_Vector;

    begin
      h := new Standard_Integer_Vectors.Vector(t.dg'range);
      for j in h'range loop
        h(j) := t.dg(j);
      end loop;
      Append(res,res_last,h);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Laurentials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Create;

  function Random_Complex_Polynomial
             ( s : Lists_of_Integer_Vectors.List )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Integer_Vectors,Lists_of_Integer_Vectors;
    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    tmp : List := s;

  begin
    while not Is_Null(tmp) loop
      declare
        v : Link_to_Vector := Head_Of(tmp);
        t : Term;
      begin
        t.cf := Standard_Random_Numbers.Random1;
        t.dg := new Standard_Natural_Vectors.Vector'(v'range => 0);
        for i in v'range loop
          if v(i) > 0
           then t.dg(i) := natural32(v(i));
          end if;
        end loop;
        Add(res,t); Clear(t);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Random_Complex_Polynomial;

  function Random_Complex_Laurential
             ( s : Lists_of_Integer_Vectors.List )
             return Standard_Complex_Laurentials.Poly is

    use Standard_Integer_Vectors,Lists_of_Integer_Vectors;
    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;
    tmp : List := s;

  begin
    while not Is_Null(tmp) loop
      declare
        v : Link_to_Vector := Head_Of(tmp);
        t : Term;
      begin
        t.cf := Standard_Random_Numbers.Random1;
        t.dg := new Standard_Integer_Vectors.Vector'(v.all);
        Add(res,t); Clear(t);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Random_Complex_Laurential;

  function Random_Complex_Laurent_System
             ( dim : integer32; s : Lists_of_Integer_Vectors.List )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(1..dim);

  begin
    for i in res'range loop
      res(i) := Random_Complex_Laurential(s);
    end loop;
    return res;
  end Random_Complex_Laurent_System;

  function Random_Complex_Laurent_System
             ( dim : integer32; mix : Standard_Integer_Vectors.Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(1..dim);
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop  -- use the i-th support mix(i) times
        ind := ind + 1;
        res(ind) := Random_Complex_Laurential(sup(i));
      end loop;
    end loop;
    return res;
  end Random_Complex_Laurent_System;

  function Is_Equal ( v1 : Standard_Floating_Vectors.Vector;
                      v2 : Standard_Natural_Vectors.Vector ) return boolean is
  begin
    for i in v2'range loop
      if double_float(v2(i)) /= v1(i)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Equal;

  function Is_Equal ( v1 : Standard_Floating_Vectors.Vector;
                      v2 : Standard_Integer_Vectors.Vector ) return boolean is
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
                          s : Lists_of_Integer_Vectors.List )
                        return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Select_Term ( t : in Standard_Complex_Polynomials.Term;
                            cont : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      for i in v'range loop
        v(i) := integer32(t.dg(i));
      end loop;
      if Lists_of_Integer_Vectors.Is_In(s,v)
       then Standard_Complex_Polynomials.Add(res,t);
      end if;
      cont := true;
    end Select_Term;
    procedure Select_Poly is
      new Standard_Complex_Polynomials.Visiting_Iterator(Select_Term);

  begin
    Select_Poly(p);
    return res;
  end Select_Terms;

  function Select_Terms ( p : Standard_Complex_Laurentials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return Standard_Complex_Laurentials.Poly is

    res : Standard_Complex_Laurentials.Poly
        := Standard_Complex_Laurentials.Null_Poly;

    procedure Select_Term ( t : in Standard_Complex_Laurentials.Term;
                            cont : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      for i in v'range loop
        v(i) := t.dg(i);
      end loop;
      if Lists_of_Integer_Vectors.Is_In(s,v)
       then Standard_Complex_Laurentials.Add(res,t);
      end if;
      cont := true;
    end Select_Term;
    procedure Select_Poly is
      new Standard_Complex_Laurentials.Visiting_Iterator(Select_Term);

  begin
    Select_Poly(p);
    return res;
  end Select_Terms;

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

  function Select_Terms ( p : Standard_Complex_Laurentials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      v : constant Standard_Integer_Vectors.Vector(t.dg'range) := t.dg.all;

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

  function Select_Terms ( p : DoblDobl_Complex_Polynomials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Select_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                            cont : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      for i in v'range loop
        v(i) := integer32(t.dg(i));
      end loop;
      if Lists_of_Integer_Vectors.Is_In(s,v)
       then DoblDobl_Complex_Polynomials.Add(res,t);
      end if;
      cont := true;
    end Select_Term;
    procedure Select_Poly is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Select_Term);

  begin
    Select_Poly(p);
    return res;
  end Select_Terms;

  function Select_Terms ( p : DoblDobl_Complex_Polynomials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;
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

  function Select_Terms ( p : DoblDobl_Complex_Laurentials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return DoblDobl_Complex_Laurentials.Poly is

    res : DoblDobl_Complex_Laurentials.Poly
        := DoblDobl_Complex_Laurentials.Null_Poly;

    procedure Select_Term ( t : in DoblDobl_Complex_Laurentials.Term;
                            cont : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      for i in v'range loop
        v(i) := integer32(t.dg(i));
      end loop;
      if Lists_of_Integer_Vectors.Is_In(s,v)
       then DoblDobl_Complex_Laurentials.Add(res,t);
      end if;
      cont := true;
    end Select_Term;
    procedure Select_Poly is
      new DoblDobl_Complex_Laurentials.Visiting_Iterator(Select_Term);

  begin
    Select_Poly(p);
    return res;
  end Select_Terms;

  function Select_Terms ( p : DoblDobl_Complex_Laurentials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return DoblDobl_Complex_Laurentials.Poly is

    use DoblDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      v : constant Standard_Integer_Vectors.Vector(t.dg'range) := t.dg.all;

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

  function Select_Terms ( p : QuadDobl_Complex_Polynomials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Select_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                            cont : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      for i in v'range loop
        v(i) := integer32(t.dg(i));
      end loop;
      if Lists_of_Integer_Vectors.Is_In(s,v)
       then QuadDobl_Complex_Polynomials.Add(res,t);
      end if;
      cont := true;
    end Select_Term;
    procedure Select_Poly is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Select_Term);

  begin
    Select_Poly(p);
    return res;
  end Select_Terms;

  function Select_Terms ( p : QuadDobl_Complex_Polynomials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;
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

  function Select_Terms ( p : QuadDobl_Complex_Laurentials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return QuadDobl_Complex_Laurentials.Poly is

    res : QuadDobl_Complex_Laurentials.Poly
        := QuadDobl_Complex_Laurentials.Null_Poly;

    procedure Select_Term ( t : in QuadDobl_Complex_Laurentials.Term;
                            cont : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      for i in v'range loop
        v(i) := integer32(t.dg(i));
      end loop;
      if Lists_of_Integer_Vectors.Is_In(s,v)
       then QuadDobl_Complex_Laurentials.Add(res,t);
      end if;
      cont := true;
    end Select_Term;
    procedure Select_Poly is
      new QuadDobl_Complex_Laurentials.Visiting_Iterator(Select_Term);

  begin
    Select_Poly(p);
    return res;
  end Select_Terms;

  function Select_Terms ( p : QuadDobl_Complex_Laurentials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      v : constant Standard_Integer_Vectors.Vector(t.dg'range) := t.dg.all;

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

  function Create ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    for i in p'range loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Create;

  function Create ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    for i in p'range loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    for i in p'range loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    for i in p'range loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Create;

  function Create ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    for i in p'range loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Create;

  function Create ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    for i in p'range loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Create;

  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
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

  function Select_Terms ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
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

  function Select_Terms ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Laur_Systems.Laur_Sys is

    use Standard_Complex_Laur_Systems;
    res : Laur_Sys(p'range);

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

  function Select_Terms ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Laur_Systems.Laur_Sys is

    use Standard_Complex_Laur_Systems;
    res : Laur_Sys(p'range);
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

  function Select_Terms ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
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

  function Select_Terms ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Poly_Systems;
    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    use DoblDobl_Complex_Laur_Systems;
    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector; 
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Poly_Systems;
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

  function Select_Terms ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    use DoblDobl_Complex_Laur_Systems;
    res : Laur_Sys(p'range);
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

  function Select_Terms ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
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

  function Select_Terms ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Poly_Systems;
    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists ) 
                        return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    use QuadDobl_Complex_Laur_Systems;
    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),s(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector; 
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Poly_Systems;
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

  function Select_Terms ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    use QuadDobl_Complex_Laur_Systems;
    res : Laur_Sys(p'range);
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

  function Select_Lifted
               ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                 mix : Standard_Integer_Vectors.Vector;
                 lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
               return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    intsup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(lifsup'range);

  begin
    intsup := Floating_Integer_Convertors.Convert(lifsup);
    for i in intsup'range loop
      Transforming_Integer32_Vector_Lists.Reduce(intsup(i),p'last+1);
    end loop;
    res := Select_Terms(p,mix,intsup);
    Arrays_of_Integer_Vector_Lists.Deep_Clear(intsup);
    return res;
  end Select_Lifted;

  function Select_Lifted
               ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                 mix : Standard_Integer_Vectors.Vector;
                 lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
               return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    intsup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(lifsup'range);

  begin
    intsup := Floating_Integer_Convertors.Convert(lifsup);
    for i in intsup'range loop
      Transforming_Integer32_Vector_Lists.Reduce(intsup(i),p'last+1);
    end loop;
    res := Select_Terms(p,mix,intsup);
    Arrays_of_Integer_Vector_Lists.Deep_Clear(intsup);
    return res;
  end Select_Lifted;

-- PROCEDURES without extra memory allocation :

  procedure Select_Coefficients
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Lists_of_Integer_Vectors.List;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : out Standard_Complex_Vectors.Vector ) is

    ind : integer32 := 0;
    tmp : Lists_of_Integer_Vectors.List := s;
    ls : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      ls := Lists_of_Integer_Vectors.Head_Of(tmp);
      for i in 1..integer32(dim) loop
        deg(i) := natural32(ls(i));
      end loop;
      ind := ind + 1;
      cff(ind) := Standard_Complex_Polynomials.Coeff(p,deg);
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Laurentials.Poly;
                s : in Lists_of_Integer_Vectors.List;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : out Standard_Complex_Vectors.Vector ) is

    ind : integer32 := 0;
    tmp : Lists_of_Integer_Vectors.List := s;
    ls : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      ls := Lists_of_Integer_Vectors.Head_Of(tmp);
      for i in 1..integer32(dim) loop
        deg(i) := ls(i);
      end loop;
      ind := ind + 1;
      cff(ind) := Standard_Complex_Laurentials.Coeff(p,deg);
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Lists_of_Floating_Vectors.List;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : out Standard_Complex_Vectors.Vector ) is

    ind : integer32 := 0;
    tmp : Lists_of_Floating_Vectors.List := s;
    ls : Standard_Floating_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
      ls := Lists_of_Floating_Vectors.Head_Of(tmp);
      for i in 1..integer32(dim) loop
        deg(i) := natural32(ls(i));
      end loop;
      ind := ind + 1;
      cff(ind) := Standard_Complex_Polynomials.Coeff(p,deg);
      tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Laurentials.Poly;
                s : in Lists_of_Floating_Vectors.List;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : out Standard_Complex_Vectors.Vector ) is

    ind : integer32 := 0;
    tmp : Lists_of_Floating_Vectors.List := s;
    ls : Standard_Floating_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
      ls := Lists_of_Floating_Vectors.Head_Of(tmp);
      for i in 1..integer32(dim) loop
        deg(i) := integer32(ls(i));
      end loop;
      ind := ind + 1;
      cff(ind) := Standard_Complex_Laurentials.Coeff(p,deg);
      tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec ) is
  begin
    for i in p'range loop
      Select_Coefficients(p(i),s(i),dim,deg,cff(i).all);
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec ) is
  begin
    for i in p'range loop
      Select_Coefficients(p(i),s(i),dim,deg,cff(i).all);
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec ) is
  begin
    for i in p'range loop
      Select_Coefficients(p(i),s(i),dim,deg,cff(i).all);
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec ) is
  begin
    for i in p'range loop
      Select_Coefficients(p(i),s(i),dim,deg,cff(i).all);
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                m : in Standard_Integer_Vectors.Vector;
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec ) is

    ind : integer32 := 0;

  begin
    for i in s'range loop
      for j in 1..m(i) loop
        ind := ind + 1;
        Select_Coefficients(p(ind),s(i),dim,deg,cff(ind).all);
      end loop;
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                m : in Standard_Integer_Vectors.Vector;
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec ) is

    ind : integer32 := 0;

  begin
    for i in s'range loop
      for j in 1..m(i) loop
        ind := ind + 1;
        Select_Coefficients(p(ind),s(i),dim,deg,cff(ind).all);
      end loop;
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                m : in Standard_Integer_Vectors.Vector;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec ) is

    ind : integer32 := 0;

  begin
    for i in s'range loop
      for j in 1..m(i) loop
        ind := ind + 1;
        Select_Coefficients(p(ind),s(i),dim,deg,cff(ind).all);
      end loop;
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                m : in Standard_Integer_Vectors.Vector;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec ) is

    ind : integer32 := 0;

  begin
    for i in s'range loop
      for j in 1..m(i) loop
        ind := ind + 1;
        Select_Coefficients(p(ind),s(i),dim,deg,cff(ind).all);
      end loop;
    end loop;
  end Select_Coefficients;

end Supports_of_Polynomial_Systems;
