with Standard_Complex_Series;
with DoblDobl_Complex_Series;
with TripDobl_Complex_Series;
with QuadDobl_Complex_Series;
with PentDobl_Complex_Series;
with OctoDobl_Complex_Series;
with DecaDobl_Complex_Series;
with HexaDobl_Complex_Series;

package body Series_and_Solutions is

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return Standard_Complex_Series_Vectors.Vector is

    nsl : constant integer32 := sol'last;
    dim : constant integer32 := (if idx = 0 then nsl else nsl-1);
    res : Standard_Complex_Series_Vectors.Vector(1..dim);

  begin
    if idx = 0 then
      for k in res'range loop
        res(k) := Standard_Complex_Series.Create(sol(k));
      end loop;
    else
      for k in 1..(idx-1) loop
        res(k) := Standard_Complex_Series.Create(sol(k));
      end loop;
      for k in (idx+1)..nsl loop
        res(k-1) := Standard_Complex_Series.Create(sol(k));
      end loop;
    end if;
    return res;
  end Create;

  function Create ( sol : DoblDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return DoblDobl_Complex_Series_Vectors.Vector is

    nsl : constant integer32 := sol'last;
    dim : constant integer32 := (if idx = 0 then nsl else nsl-1);
    res : DoblDobl_Complex_Series_Vectors.Vector(1..dim);

  begin
    if idx = 0 then
      for k in res'range loop
        res(k) := DoblDobl_Complex_Series.Create(sol(k));
      end loop;
    else
      for k in 1..(idx-1) loop
        res(k) := DoblDobl_Complex_Series.Create(sol(k));
      end loop;
      for k in (idx+1)..nsl loop
        res(k-1) := DoblDobl_Complex_Series.Create(sol(k));
      end loop;
    end if;
    return res;
  end Create;

  function Create ( sol : TripDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return TripDobl_Complex_Series_Vectors.Vector is

    nsl : constant integer32 := sol'last;
    dim : constant integer32 := (if idx = 0 then nsl else nsl-1);
    res : TripDobl_Complex_Series_Vectors.Vector(1..dim);

  begin
    if idx = 0 then
      for k in res'range loop
        res(k) := TripDobl_Complex_Series.Create(sol(k));
      end loop;
    else
      for k in 1..(idx-1) loop
        res(k) := TripDobl_Complex_Series.Create(sol(k));
      end loop;
      for k in (idx+1)..nsl loop
        res(k-1) := TripDobl_Complex_Series.Create(sol(k));
      end loop;
    end if;
    return res;
  end Create;

  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return QuadDobl_Complex_Series_Vectors.Vector is

    nsl : constant integer32 := sol'last;
    dim : constant integer32 := (if idx = 0 then nsl else nsl-1);
    res : QuadDobl_Complex_Series_Vectors.Vector(1..dim);

  begin
    if idx = 0 then
      for k in res'range loop
        res(k) := QuadDobl_Complex_Series.Create(sol(k));
      end loop;
    else
      for k in 1..(idx-1) loop
        res(k) := QuadDobl_Complex_Series.Create(sol(k));
      end loop;
      for k in (idx+1)..nsl loop
        res(k-1) := QuadDobl_Complex_Series.Create(sol(k));
      end loop;
    end if;
    return res;
  end Create;

  function Create ( sol : PentDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return PentDobl_Complex_Series_Vectors.Vector is

    nsl : constant integer32 := sol'last;
    dim : constant integer32 := (if idx = 0 then nsl else nsl-1);
    res : PentDobl_Complex_Series_Vectors.Vector(1..dim);

  begin
    if idx = 0 then
      for k in res'range loop
        res(k) := PentDobl_Complex_Series.Create(sol(k));
      end loop;
    else
      for k in 1..(idx-1) loop
        res(k) := PentDobl_Complex_Series.Create(sol(k));
      end loop;
      for k in (idx+1)..nsl loop
        res(k-1) := PentDobl_Complex_Series.Create(sol(k));
      end loop;
    end if;
    return res;
  end Create;

  function Create ( sol : OctoDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return OctoDobl_Complex_Series_Vectors.Vector is

    nsl : constant integer32 := sol'last;
    dim : constant integer32 := (if idx = 0 then nsl else nsl-1);
    res : OctoDobl_Complex_Series_Vectors.Vector(1..dim);

  begin
    if idx = 0 then
      for k in res'range loop
        res(k) := OctoDobl_Complex_Series.Create(sol(k));
      end loop;
    else
      for k in 1..(idx-1) loop
        res(k) := OctoDobl_Complex_Series.Create(sol(k));
      end loop;
      for k in (idx+1)..nsl loop
        res(k-1) := OctoDobl_Complex_Series.Create(sol(k));
      end loop;
    end if;
    return res;
  end Create;

  function Create ( sol : DecaDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return DecaDobl_Complex_Series_Vectors.Vector is

    nsl : constant integer32 := sol'last;
    dim : constant integer32 := (if idx = 0 then nsl else nsl-1);
    res : DecaDobl_Complex_Series_Vectors.Vector(1..dim);

  begin
    if idx = 0 then
      for k in res'range loop
        res(k) := DecaDobl_Complex_Series.Create(sol(k));
      end loop;
    else
      for k in 1..(idx-1) loop
        res(k) := DecaDobl_Complex_Series.Create(sol(k));
      end loop;
      for k in (idx+1)..nsl loop
        res(k-1) := DecaDobl_Complex_Series.Create(sol(k));
      end loop;
    end if;
    return res;
  end Create;

  function Create ( sol : HexaDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return HexaDobl_Complex_Series_Vectors.Vector is

    nsl : constant integer32 := sol'last;
    dim : constant integer32 := (if idx = 0 then nsl else nsl-1);
    res : HexaDobl_Complex_Series_Vectors.Vector(1..dim);

  begin
    if idx = 0 then
      for k in res'range loop
        res(k) := HexaDobl_Complex_Series.Create(sol(k));
      end loop;
    else
      for k in 1..(idx-1) loop
        res(k) := HexaDobl_Complex_Series.Create(sol(k));
      end loop;
      for k in (idx+1)..nsl loop
        res(k-1) := HexaDobl_Complex_Series.Create(sol(k));
      end loop;
    end if;
    return res;
  end Create;

  function Create ( sol : Standard_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return Standard_Complex_Series_Vectors.Vector is
  begin
    return Create(sol.v,idx);
  end Create;

  function Create ( sol : DoblDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return DoblDobl_Complex_Series_Vectors.Vector is
  begin
    return Create(sol.v,idx);
  end Create;

  function Create ( sol : TripDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return TripDobl_Complex_Series_Vectors.Vector is
  begin
    return Create(sol.v,idx);
  end Create;

  function Create ( sol : QuadDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return QuadDobl_Complex_Series_Vectors.Vector is
  begin
    return Create(sol.v,idx);
  end Create;

  function Create ( sol : PentDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return PentDobl_Complex_Series_Vectors.Vector is
  begin
    return Create(sol.v,idx);
  end Create;

  function Create ( sol : OctoDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return OctoDobl_Complex_Series_Vectors.Vector is
  begin
    return Create(sol.v,idx);
  end Create;

  function Create ( sol : DecaDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return DecaDobl_Complex_Series_Vectors.Vector is
  begin
    return Create(sol.v,idx);
  end Create;

  function Create ( sol : HexaDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return HexaDobl_Complex_Series_Vectors.Vector is
  begin
    return Create(sol.v,idx);
  end Create;

  function Create ( sols : Standard_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return Standard_Complex_Series_VecVecs.VecVec is

    dim : constant integer32
        := integer32(Standard_Complex_Solutions.Length_Of(sols));
    res : Standard_Complex_Series_VecVecs.VecVec(1..dim);
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    for i in res'range loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      declare
        sersol : constant Standard_Complex_Series_Vectors.Vector
               := Create(ls.all,idx);
      begin
        res(i) := new Standard_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( sols : DoblDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return DoblDobl_Complex_Series_VecVecs.VecVec is

    dim : constant integer32
        := integer32(Dobldobl_Complex_Solutions.Length_Of(sols));
    res : DoblDobl_Complex_Series_VecVecs.VecVec(1..dim);
    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    for i in res'range loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      declare
        sersol : constant DoblDobl_Complex_Series_Vectors.Vector
               := Create(ls.all,idx);
      begin
        res(i) := new DoblDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( sols : TripDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return TripDobl_Complex_Series_VecVecs.VecVec is

    dim : constant integer32
        := integer32(TripDobl_Complex_Solutions.Length_Of(sols));
    res : TripDobl_Complex_Series_VecVecs.VecVec(1..dim);
    tmp : TripDobl_Complex_Solutions.Solution_List := sols;
    ls : TripDobl_Complex_Solutions.Link_to_Solution;

  begin
    for i in res'range loop
      ls := TripDobl_Complex_Solutions.Head_Of(tmp);
      declare
        sersol : constant TripDobl_Complex_Series_Vectors.Vector
               := Create(ls.all,idx);
      begin
        res(i) := new TripDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := TripDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( sols : QuadDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return QuadDobl_Complex_Series_VecVecs.VecVec is

    dim : constant integer32
        := integer32(QuadDobl_Complex_Solutions.Length_Of(sols));
    res : QuadDobl_Complex_Series_VecVecs.VecVec(1..dim);
    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    for i in res'range loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      declare
        sersol : constant QuadDobl_Complex_Series_Vectors.Vector
               := Create(ls.all,idx);
      begin
        res(i) := new QuadDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( sols : PentDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return PentDobl_Complex_Series_VecVecs.VecVec is

    dim : constant integer32
        := integer32(PentDobl_Complex_Solutions.Length_Of(sols));
    res : PentDobl_Complex_Series_VecVecs.VecVec(1..dim);
    tmp : PentDobl_Complex_Solutions.Solution_List := sols;
    ls : PentDobl_Complex_Solutions.Link_to_Solution;

  begin
    for i in res'range loop
      ls := PentDobl_Complex_Solutions.Head_Of(tmp);
      declare
        sersol : constant PentDobl_Complex_Series_Vectors.Vector
               := Create(ls.all,idx);
      begin
        res(i) := new PentDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := PentDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( sols : OctoDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return OctoDobl_Complex_Series_VecVecs.VecVec is

    dim : constant integer32
        := integer32(OctoDobl_Complex_Solutions.Length_Of(sols));
    res : OctoDobl_Complex_Series_VecVecs.VecVec(1..dim);
    tmp : OctoDobl_Complex_Solutions.Solution_List := sols;
    ls : OctoDobl_Complex_Solutions.Link_to_Solution;

  begin
    for i in res'range loop
      ls := OctoDobl_Complex_Solutions.Head_Of(tmp);
      declare
        sersol : constant OctoDobl_Complex_Series_Vectors.Vector
               := Create(ls.all,idx);
      begin
        res(i) := new OctoDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := OctoDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( sols : DecaDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return DecaDobl_Complex_Series_VecVecs.VecVec is

    dim : constant integer32
        := integer32(DecaDobl_Complex_Solutions.Length_Of(sols));
    res : DecaDobl_Complex_Series_VecVecs.VecVec(1..dim);
    tmp : DecaDobl_Complex_Solutions.Solution_List := sols;
    ls : DecaDobl_Complex_Solutions.Link_to_Solution;

  begin
    for i in res'range loop
      ls := DecaDobl_Complex_Solutions.Head_Of(tmp);
      declare
        sersol : constant DecaDobl_Complex_Series_Vectors.Vector
               := Create(ls.all,idx);
      begin
        res(i) := new DecaDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := DecaDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( sols : HexaDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return HexaDobl_Complex_Series_VecVecs.VecVec is

    dim : constant integer32
        := integer32(HexaDobl_Complex_Solutions.Length_Of(sols));
    res : HexaDobl_Complex_Series_VecVecs.VecVec(1..dim);
    tmp : HexaDobl_Complex_Solutions.Solution_List := sols;
    ls : HexaDobl_Complex_Solutions.Link_to_Solution;

  begin
    for i in res'range loop
      ls := HexaDobl_Complex_Solutions.Head_Of(tmp);
      declare
        sersol : constant HexaDobl_Complex_Series_Vectors.Vector
               := Create(ls.all,idx);
      begin
        res(i) := new HexaDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := HexaDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Create;

end Series_and_Solutions;
