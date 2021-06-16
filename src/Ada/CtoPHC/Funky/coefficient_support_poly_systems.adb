with Interfaces.C;                      use Interfaces.C;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Coefficient_Support_Polynomials;   use Coefficient_Support_Polynomials;

package body Coefficient_Support_Poly_Systems is

  function Monomial_Count ( p : Poly_Sys ) return C_Integer_Array is

    res : C_Integer_Array(0..Interfaces.C.size_T(p'length-1));
    ind : integer32 := p'first;

  begin
    for i in res'range loop
      res(i) := Interfaces.C.int(Number_of_Terms(p(ind)));
      ind := ind+1;
    end loop;
    return res;
  end Monomial_Count;

  function Sum ( a : C_Integer_Array ) return integer32 is

    res : integer32 := 0;

  begin
    for i in a'range loop
      res := res + integer32(a(i));
    end loop;
    return res;
  end Sum;

  function Support ( n,m : natural32; moncnt : C_Integer_Array;
                     p : Poly_Sys ) return C_Integer_Array is

    res : C_Integer_Array(0..Interfaces.C.size_T(n*m-1));
    indmon : Interfaces.C.size_T := 0;
    indres : Interfaces.C.size_T := 0;

  begin
    for i in p'range loop
      declare
        dimsup : constant natural32 := natural32(moncnt(indmon))*n;
        sup : constant C_Integer_Array(0..Interfaces.C.size_T(dimsup-1))
            := Support(p(i));
      begin
        for j in sup'range loop
          res(indres) := sup(j);
          indres := indres + 1;
        end loop;
      end;
      indmon := indmon + 1;
    end loop;
    return res;
  end Support;

  function Coefficients ( m : natural32; moncnt : C_Integer_Array;
                          p : Poly_Sys ) return C_Double_Array is

    res : C_Double_Array(0..Interfaces.C.size_T(2*m-1));
    indmon : Interfaces.C.size_T := 0;
    indres : Interfaces.C.size_T := 0;

  begin
    for i in p'range loop
      declare
        dimcff : constant natural32 := 2*natural32(moncnt(indmon));
        cff : constant C_Double_Array(0..Interfaces.C.size_T(dimcff-1))
            := Coefficients(p(i));
      begin
        for j in cff'range loop
          res(indres) := cff(j);
          indres := indres + 1;
        end loop;
      end;
      indmon := indmon + 1;
    end loop;
    return res;
  end Coefficients;

  function Create ( n : natural32; m : C_Integer_Array;
                    c : C_Double_Array; s : C_Integer_Array )
                  return Poly_Sys is

    res : Poly_Sys(1..m'length);
    indmon : Interfaces.C.size_T := m'first;
    indcff : Interfaces.C.size_T := c'first;
    indsup : Interfaces.C.size_T := s'first;

  begin
    for i in res'range loop
      declare
        dimsup : constant natural32 := n*natural32(m(indmon));
        sup : C_Integer_Array(0..Interfaces.C.size_T(dimsup-1));
        dimcff : constant natural32 := 2*natural32(m(indmon));
        cff : C_Double_Array(0..Interfaces.C.size_T(dimcff-1));
      begin
        for j in cff'range loop
          cff(j) := c(indcff);
          indcff := indcff + 1;
        end loop;
        for j in sup'range loop
          sup(j) := s(indsup);
          indsup := indsup + 1;
        end loop;
        res(i) := Create(n,cff,sup);
      end;
      indmon := indmon + 1;
    end loop;
    return res;
  end Create;

  function Concat ( n : natural32; m : C_Integer_Array;
                    c : C_Double_Array; s : C_Integer_Array )
                  return C_Double_Array is

    d : constant size_T
      := 11 + m'last-m'first+1 + s'last-s'first+1 + c'last-c'first+1;
    res : C_Double_Array(0..d);
    ind : size_T;

  begin
    res(0) := Interfaces.C.double(d);
    res(1) := Interfaces.C.double(n);
    res(2) := Interfaces.C.double(m'first);
    res(3) := Interfaces.C.double(m'last);
    res(4) := Interfaces.C.double(s'first);
    res(5) := Interfaces.C.double(s'last);
    res(6) := Interfaces.C.double(c'first);
    res(7) := Interfaces.C.double(c'last);
    ind := 11;
    res(8) := Interfaces.C.double(ind);    -- start of monomial count 
    for i in m'range loop
      res(ind) := Interfaces.C.double(m(i));
      ind := ind + 1;
    end loop;
    res(9) := Interfaces.C.double(ind);    -- start of support
    for i in s'range loop
      res(ind) := Interfaces.C.double(s(i));
      ind := ind + 1;
    end loop;
    res(10) := Interfaces.C.double(ind);   -- start of coefficients
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

  function Monomial_Count ( x : C_Double_Array ) return C_Integer_Array is

    i1 : constant size_T := Interfaces.C.size_T(x(2));
    i2 : constant size_T := Interfaces.C.size_T(x(3));
    res : C_Integer_Array(i1..i2);
    ind : size_T := Interfaces.C.size_T(x(8));

  begin
    for i in res'range loop
      res(i) := Interfaces.C.int(x(ind));
      ind := ind + 1;
    end loop;
    return res;
  end Monomial_Count;

  function Support ( x : C_Double_Array ) return C_Integer_Array is

    i1 : constant size_T := Interfaces.C.size_T(x(4));
    i2 : constant size_T := Interfaces.C.size_T(x(5));
    res : C_Integer_Array(i1..i2);
    ind : size_T := Interfaces.C.size_T(x(9));

  begin
    for i in res'range loop
      res(i) := Interfaces.C.int(x(ind));
      ind := ind + 1;
    end loop;
    return res;
  end Support;

  function Coefficients ( x : C_Double_Array ) return C_Double_Array is

    i1 : constant size_T := Interfaces.C.size_T(x(6));
    i2 : constant size_T := Interfaces.C.size_T(x(7));
    res : C_Double_Array(i1..i2);
    ind : size_T := Interfaces.C.size_T(x(10));

  begin
    for i in res'range loop
      res(i) := x(ind);
      ind := ind + 1;
    end loop;
    return res;
  end Coefficients;

end Coefficient_Support_Poly_Systems;
