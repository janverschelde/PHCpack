with Interfaces.C;                      use Interfaces.C;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;

package body Coefficient_Solution_Vectors is

  function Multiplicity ( s : Solution ) return integer32 is
  begin
    return s.m;
  end Multiplicity;

  function Multiplicities ( s : Solution_List ) return C_Integer_Array is

    len : constant Interfaces.C.size_T := Interfaces.C.size_T(Length_Of(s));
    res : C_Integer_Array(1..len);
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      res(i) := Interfaces.C.int(Multiplicity(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Multiplicities;

  function Coefficients ( s : Solution ) return C_Double_Array is

    len : constant natural32 := 2*natural32(s.n) + 5;
    res : C_Double_Array(1..Interfaces.C.size_T(len));
    ind : Interfaces.C.size_T := 2;

  begin
    res(1) := Interfaces.C.double(REAL_PART(s.t));
    res(2) := Interfaces.C.double(IMAG_PART(s.t));
    for i in s.v'range loop
      ind := ind + 1;
      res(ind) := Interfaces.C.double(REAL_PART(s.v(i)));
      ind := ind + 1;
      res(ind) := Interfaces.C.double(IMAG_PART(s.v(i)));
    end loop;
    ind := ind + 1;
    res(ind) := Interfaces.C.double(s.err);
    ind := ind + 1;
    res(ind) := Interfaces.C.double(s.rco);
    ind := ind + 1;
    res(ind) := Interfaces.C.double(s.res);
    return res;
  end Coefficients;

  function Coefficients ( s : Solution_List ) return C_Double_Array is

    n : constant natural32 := natural32(Head_Of(s).n);
    len : constant natural32 := Length_Of(s)*(2*n+5);
    res : C_Double_Array(1..Interfaces.C.size_T(len));
    tmp : Solution_List := s;
    ls : Link_to_Solution;
    ind : Interfaces.C.size_T := 0;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        cff : constant C_Double_Array := Coefficients(ls.all);
      begin
        for i in cff'range loop
          ind := ind + 1;
          res(ind) := cff(i);
        end loop;
      end;
      tmp := Tail_Of(tmp);
    end loop; 
    return res;
  end Coefficients; 

  function Concat ( n : natural32; m : C_Integer_Array; c : C_Double_Array )
                  return C_Double_Array is

    d : constant size_T := 8 + m'last-m'first+1 + c'last-c'first+1;
    res : C_Double_Array(0..d);
    ind : size_T;

  begin
    res(0) := Interfaces.C.double(d);
    res(1) := Interfaces.C.double(n);
    res(2) := Interfaces.C.double(m'first);
    res(3) := Interfaces.C.double(m'last);
    res(4) := Interfaces.C.double(c'first);
    res(5) := Interfaces.C.double(c'last);
    ind := 8;
    res(6) := Interfaces.C.double(ind);   -- start of multiplicities
    for i in m'range loop
      res(ind) := Interfaces.C.double(m(i));
      ind := ind + 1;
    end loop;
    res(7) := Interfaces.C.double(ind);   -- start of coefficients
    for i in c'range loop
      res(ind) := c(i);
      ind := ind + 1;
    end loop;
    return res;
  end Concat;

  function Dimension ( x : C_Double_Array ) return natural32 is
  begin
    return natural32(x(1));
  end Dimension;

  function Multiplicities ( x : C_Double_Array ) return C_Integer_Array is

    i1 : constant size_T := Interfaces.C.size_t(x(2));
    i2 : constant size_T := Interfaces.C.size_t(x(3));
    res : C_Integer_Array(i1..i2);
    ind : size_T := Interfaces.C.size_T(x(6));

  begin
    for i in res'range loop
      res(i) := Interfaces.C.int(x(ind));
      ind := ind + 1;
    end loop;
    return res;
  end Multiplicities;

  function Coefficients ( x : C_Double_Array ) return C_Double_Array is

    i1 : constant size_T := Interfaces.C.size_t(x(4));
    i2 : constant size_T := Interfaces.C.size_t(x(5));
    res : C_Double_Array(i1..i2);
    ind : size_T := Interfaces.C.size_T(x(7));

  begin
    for i in res'range loop
      res(i) := x(ind);
      ind := ind + 1;
    end loop;
    return res;
  end Coefficients;

  function Create ( n,m : natural32; c : C_Double_Array ) return Solution is

    res : Solution(integer32(n));
    ind : Interfaces.C.size_T := c'first;

  begin
    res.t := Create(double_float(c(ind)),double_float(c(ind+1)));
    res.m := integer32(m);
    for i in 1..integer32(n) loop
      ind := ind + 2;
      res.v(i) := Create(double_float(c(ind)),double_float(c(ind+1)));
    end loop;
    ind := ind + 2;
    res.err := double_float(c(ind));
    res.rco := double_float(c(ind+1));
    res.res := double_float(c(ind+2));
    return res;
  end Create;

  function Create ( n : natural32; m : C_Integer_Array; c : C_Double_Array )
                  return Solution_List is

    res,res_last : Solution_List;
    ind : Interfaces.C.size_T := c'first;

  begin
    for i in m'range loop
      declare
        s : Solution(integer32(n));
      begin
        s.m := integer32(m(i));
        s.t := Create(double_float(c(ind)),double_float(c(ind+1)));
        for i in 1..integer32(n) loop
          ind := ind + 2;
          s.v(i) := Create(double_float(c(ind)),double_float(c(ind+1)));
        end loop;
        ind := ind + 2;
        s.err := double_float(c(ind));
        s.rco := double_float(c(ind+1));
        s.res := double_float(c(ind+2));
        ind := ind + 3;
        Append(res,res_last,s);
      end;
    end loop;
    return res;
  end Create;

end Coefficient_Solution_Vectors;
