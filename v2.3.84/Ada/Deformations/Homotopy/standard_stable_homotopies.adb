with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;

package body Standard_Stable_Homotopies is

  function Number_of_Zeroes
             ( z : Standard_Integer_Vectors.Vector ) return integer32 is

    res : integer32 := 0;

  begin
    for i in z'range loop
      if z(i) < 0 then
        return -1;
      elsif z(i) = 0 then
        res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Zeroes;

  procedure Zero_Type
              ( v : in Standard_Complex_Vectors.Vector; nz : out integer32;
                z : out Standard_Integer_Vectors.Vector ) is
  begin
    nz := 0;
    for i in v'range loop
      if ((REAL_PART(v(i)) = 0.0) and then (IMAG_PART(v(i)) = 0.0))
       then z(i) := 0; nz := nz + 1;
       else z(i) := 1;
      end if;
    end loop;
  end Zero_Type;

  function Remove_Zeroes
             ( s : in Solution; nz : in integer32;
               z : in Standard_Integer_Vectors.Vector ) return Solution is

    res : Solution(s.n-nz);
    ind : integer32 := 0;

  begin
    res.t := s.t;
    res.m := s.m;
    for i in s.v'range loop
      if z(i) /= 0 then
        ind := ind + 1;
        res.v(ind) := s.v(i);
      end if;
    end loop;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Remove_Zeroes;

  function Origin ( n,m : integer32 ) return Solution is

    res : Solution(n);

  begin
    res.t := Create(1.0);
    res.m := m;
    res.v := (1..n => Create(0.0));
    res.err := 0.0;
    res.rco := 1.0;
    res.res := 0.0;
    return res;
  end Origin;

  function Insert_Zeroes
              ( v : Standard_Complex_Vectors.Vector;
                z : Standard_Integer_Vectors.Vector )
              return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(z'range);
    ind : integer32 := v'first;

  begin
    for i in z'range loop
      if z(i) = 0 then
        res(i) := Create(0.0);
      else
        res(i) := v(ind);
        ind := ind + 1;
      end if;
    end loop;
    return res;
  end Insert_Zeroes;

  function Insert_Zeroes
              ( s : Solution; z : Standard_Integer_Vectors.Vector )
              return Solution is

    res : Solution(z'last);

  begin
    res.t := s.t;
    res.m := s.m;
    res.v := Insert_Zeroes(s.v,z);
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Insert_Zeroes;

  function Insert_Zeroes
              ( s : Solution_List; z : Standard_Integer_Vectors.Vector )
              return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Insert_Zeroes(ls.all,z));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Insert_Zeroes;

  function Is_Same ( s1,s2 : Link_to_Solution ) return boolean is

  -- DESCRIPTION :
  --   Returns true if both solutions are the same.

    nz1,nz2 : integer32;
    z1 : Standard_Integer_Vectors.Vector(s1.v'range);
    z2 : Standard_Integer_Vectors.Vector(s2.v'range);

  begin
    Zero_Type(s1.v,nz1,z1);
    Zero_Type(s2.v,nz2,z2);
    if nz1 /= nz2 then
      return false;
    elsif not Standard_Integer_Vectors.Equal(z1,z2) then
      return false;
    else
      return (nz1 = s1.n);
    end if;
  end Is_Same;

  procedure Merge_and_Concat
              ( first,last : in out Solution_List;
                sols : in Solution_List ) is

    tmp : Solution_List := sols;
    ptr : Solution_List;
    ls1,ls2 : Link_to_Solution;
    found : boolean;

  begin
    while not Is_Null(tmp) loop
      ls1 := Head_Of(tmp);
      ptr := first;
      found := false;
      while not Is_Null(ptr) loop
        ls2 := Head_Of(ptr);
        found := Is_Same(ls1,ls2);
        exit when found;
        ptr := Tail_Of(ptr);
      end loop;
      if found then
        ls2.m := ls2.m + ls1.m;
        Set_Head(ptr,ls2);
      else
        Append(first,last,ls1.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Merge_and_Concat;

  function Vanish_by_Zeroes
             ( t : Standard_Complex_Polynomials.Term;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return boolean is
  begin
    if nbz <= 0 then
      return false;
    else
      for i in z'range loop
        if ((z(i) = 0) and then (t.dg(i) /= 0))
         then return true;
        end if;
      end loop;
      return false;
    end if;
  end Vanish_by_Zeroes;

  function Vanish_by_Zeroes
             ( t : Standard_Complex_Laurentials.Term;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return boolean is
  begin
    if nbz <= 0 then
      return false;
    else
      for i in z'range loop
        if ((z(i) = 0) and then (t.dg(i) /= 0))
         then return true;
        end if;
      end loop;
      return false;
    end if;
  end Vanish_by_Zeroes;

  function Substitute_Zeroes
             ( t : Standard_Complex_Polynomials.Term;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;
    ind : integer32;

  begin
    if nbz <= 0 then
      return t;
    elsif Vanish_by_Zeroes(t,z,nbz) then
      res.cf := Create(0.0);
      return res;
    else
      res.cf := t.cf;
      res.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-nbz);
      ind := t.dg'first - 1;
      for i in z'range loop
        if z(i) /= 0 then
          ind := ind + 1;
          res.dg(ind) := t.dg(i);
        end if;
      end loop;
      return res;
    end if;
  end Substitute_Zeroes;

  function Substitute_Zeroes
             ( t : Standard_Complex_Laurentials.Term;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return Standard_Complex_Laurentials.Term is

    res : Standard_Complex_Laurentials.Term;
    ind : integer32;

  begin
    if nbz <= 0 then
      return t;
    elsif Vanish_by_Zeroes(t,z,nbz) then
      res.cf := Create(0.0);
      return res;
    else
      res.cf := t.cf;
      res.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last-nbz);
      ind := t.dg'first - 1;
      for i in z'range loop
        if z(i) /= 0 then
          ind := ind + 1;
          res.dg(ind) := t.dg(i);
        end if;
      end loop;
      return res;
    end if;
  end Substitute_Zeroes;

  function Substitute_Zeroes
             ( p : Standard_Complex_Polynomials.Poly;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Subs ( t : in Term; continue : out boolean ) is

      nt : Term := Substitute_Zeroes(t,z,nbz);

    begin
      if REAL_PART(t.cf) /= 0.0 or IMAG_PART(t.cf) /= 0.0 then
        Add(res,nt);
        Clear(nt);
      end if;
      continue := true;
    end Subs;
    procedure Substitute is new Visiting_Iterator(Subs);

  begin
    Substitute(p);
    return res;
  end Substitute_Zeroes;

  function Substitute_Zeroes
             ( p : Standard_Complex_Laurentials.Poly;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Subs ( t : in Term; continue : out boolean ) is

      nt : Term := Substitute_Zeroes(t,z,nbz);

    begin
      if REAL_PART(t.cf) /= 0.0 or IMAG_PART(t.cf) /= 0.0 then
        Add(res,nt);
        Clear(nt);
      end if;
      continue := true;
    end Subs;
    procedure Substitute is new Visiting_Iterator(Subs);

  begin
    Substitute(p);
    return res;
  end Substitute_Zeroes;

  function Substitute_Zeroes
             ( p : Poly_Sys; z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Substitute_Zeroes(p(i),z,nbz);
    end loop;
    return res;
  end Substitute_Zeroes;

  function Substitute_Zeroes
             ( p : Laur_Sys; z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Substitute_Zeroes(p(i),z,nbz);
    end loop;
    return res;
  end Substitute_Zeroes;

  function Filter ( p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);
    cnt : integer32 := p'first-1;

  begin
    for i in p'range loop
      if Standard_Complex_Laurentials.Degree(p(i)) > 0 then
        cnt := cnt + 1;
        res(cnt) := p(i);
      end if;
    end loop;
    return res(res'first..cnt);
  end Filter;

  function Filter ( p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);
    cnt : integer32 := p'first-1;

  begin
    for i in p'range loop
      if Standard_Complex_Polynomials.Degree(p(i)) > 0 then
        cnt := cnt + 1;
        res(cnt) := p(i);
      end if;
    end loop;
    return res(res'first..cnt);
  end Filter;

end Standard_Stable_Homotopies;
