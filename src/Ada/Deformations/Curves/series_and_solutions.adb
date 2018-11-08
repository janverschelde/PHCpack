with Standard_Complex_Series;
with DoblDobl_Complex_Series;
with QuadDobl_Complex_Series;

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

  function Create ( sol : QuadDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return QuadDobl_Complex_Series_Vectors.Vector is
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

end Series_and_Solutions;
