with Standard_Complex_Poly_Lists;        use Standard_Complex_Poly_Lists;
with Standard_Linear_Product_System;

package body Standard_Complex_Prod_Planes is

  function Create return Prod_Sys is

    dim : constant integer32 
        := integer32(Standard_Linear_Product_System.Dimension);
    res,res_last : Prod_Sys(1..dim);
    ind : integer32;

    procedure Store ( h : in Standard_Complex_Vectors.Vector;
                      cont : out boolean ) is

      p : constant Poly := Standard_Linear_Product_System.Polynomial(h);

    begin
      Append(res(ind),res_last(ind),p);
      cont := true;
    end Store;
    procedure Enumerate_and_Store is
      new Standard_Linear_Product_System.Enumerate_Hyperplanes(Store);

  begin
    for i in res'range loop
      ind := i;
      Enumerate_and_Store(natural32(i));
    end loop;
    return res;
  end Create;

  function Linear_Index
             ( d : Standard_Natural_Vectors.Vector ) return integer32 is

    res : integer32 := 0;

  begin
    for i in d'range loop
      if d(i) > 1 then
        return -1;
      elsif d(i) = 1 then
        if res = 0 then
          res := i;
        else
          return -1;
        end if;
      end if;
    end loop;
    return res;
  end Linear_Index;

  procedure Hyperplane_Coefficients
              ( p : in Poly; fail : out boolean;
                h : out Standard_Complex_Vectors.Vector ) is

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      ind : constant integer32 := Linear_Index(t.dg.all);

    begin
      if ind < 0 or ind > h'last then
        fail := true;
        continue := false;
        return;
      else
        h(ind) := t.cf;
        continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    for i in h'range loop
      h(i) := Create(0.0);
    end loop;
    fail := false;
    Scan_Terms(p);
  end Hyperplane_Coefficients;

  procedure Store ( p : in Prod_Sys; fail : out boolean ) is

    q : Poly;
    tmp : Prod_Poly;
    h : Standard_Complex_Vectors.Vector(0..p'last);

  begin
    Standard_Linear_Product_System.Init(natural32(p'last));
    for i in p'range loop
      tmp := p(i);
      while not Is_Null(tmp) loop
        q := Head_Of(tmp);
        if Degree(q) > 1 then
          fail := true; return;
        else
          Hyperplane_Coefficients(q,fail,h);
          if fail then
            return;
          else 
            Standard_Linear_Product_System.Add_Hyperplane(natural32(i),h);
            tmp := Tail_Of(tmp);
          end if;
        end if;
      end loop;
    end loop;
  end Store;

  function Degrees ( p : Prod_Sys ) return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Number_of_Factors(p(i));
    end loop;
    return res;
  end Degrees;

  function Values_at_Hyperplanes
              ( i : natural32; x : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector is

    nhi : constant natural32
        := Standard_Linear_Product_System.Number_of_Hyperplanes(i);
    res : Standard_Complex_Vectors.Vector(1..integer32(nhi));
    ind : integer32 := 0;

    procedure Value ( h : in Standard_Complex_Vectors.Vector;
                      cont : out boolean ) is

      v : Complex_Number := h(0);

    begin
      for k in x'range loop
        v := v + h(k)*x(k);
      end loop;
      ind := ind + 1;
      res(ind) := v;
      cont := true;
    end Value;
    procedure Value_at_Hyperplanes is
      new Standard_Linear_Product_System.Enumerate_Hyperplanes(Value);

  begin
    Value_at_Hyperplanes(i);
    return res;
  end Values_at_Hyperplanes;

  function Eval ( i : natural32; x : Standard_Complex_Vectors.Vector )
                return Complex_Number is

    res : Complex_Number := Create(1.0);

    procedure Value ( h : in Standard_Complex_Vectors.Vector;
                      cont : out boolean ) is

      v : Complex_Number := h(0);

    begin
      for k in x'range loop
        v := v + h(k)*x(k);
      end loop;
      res := res*v;
      cont := true;
    end Value;
    procedure Value_at_Hyperplanes is
      new Standard_Linear_Product_System.Enumerate_Hyperplanes(Value);

  begin
    Value_at_Hyperplanes(i);
    return res;
  end Eval;

  function Eval ( x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := Standard_Complex_Prod_Planes.Eval(natural32(i),x);
    end loop;
    return res;
  end Eval;

  function Gradient ( i : natural32; x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(x'range) := (x'range => Create(0.0));
    nhi : constant natural32
        := Standard_Linear_Product_System.Number_of_Hyperplanes(i);
    val : constant Standard_Complex_Vectors.Vector(1..integer32(nhi))
        := Values_at_Hyperplanes(i,x);
    p : Standard_Complex_Vectors.Vector(val'range);
    ind : integer32 := 0;

    procedure Multiply ( h : in Standard_Complex_Vectors.Vector;
                         cont : out boolean ) is

    -- DESCRIPTION :
    --   Multiplies the product vector p with the coefficients
    --   of a factor in the i-th linear product polynomial.

    begin
      ind := ind + 1;
      for i in res'range loop
        res(i) := res(i) + h(i)*p(ind);
      end loop;
      cont := true;
    end Multiply;
    procedure Multiply_with_Coefficients is 
      new Standard_Linear_Product_System.Enumerate_Hyperplanes(Multiply);

  begin
    for j in p'range loop
      p(j) := Create(1.0);
      for k in val'range loop
        if k /= j 
         then p(j) := p(j)*val(k);
        end if;
      end loop;
    end loop;
    Multiply_with_Coefficients(i);
    return res;
  end Gradient;

  function Jacobian ( x : Standard_Complex_Vectors.Vector ) return Matrix is

    res : Matrix(x'range,x'range);
    gx : Standard_Complex_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      gx := Gradient(natural32(i),x);
      for j in x'range loop
        res(i,j) := gx(j);
      end loop;
    end loop;
    return res;
  end Jacobian;

end Standard_Complex_Prod_Planes;
