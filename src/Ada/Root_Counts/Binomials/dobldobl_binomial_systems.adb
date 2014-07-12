with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Double_Double_Numbers;            use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;         use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Exponentiation;  use DoblDobl_Complex_Exponentiation;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;

package body DoblDobl_Binomial_Systems is

-- FORMAT of a BINOMIAL SYSTEM :  p(x) = 0 => x^A = c

  procedure Parse ( p : in Poly_Sys; nq : in integer32;
                    A : out Standard_Integer64_Matrices.Matrix;
                    c : out Vector; fail : out boolean ) is

    use DoblDobl_Complex_Polynomials;

    equ,mon : integer32;

    procedure Store ( t : in Term; continue : out boolean ) is
    begin
      if mon = 1 then
        for i in t.dg'range loop     -- store monomial in column equ
          A(i,equ) := integer64(t.dg(i));
        end loop;
        c(equ) := t.cf;              -- store coefficient at equ
        mon := mon + 1;
      else
        for i in t.dg'range loop     -- divide by new monomial
          A(i,equ) := A(i,equ) - integer64(t.dg(i));
        end loop;
        c(equ) := -t.cf/c(equ);      -- coefficient to righthand side
      end if;
      continue := true;
    end Store;
    procedure Store_Terms is new Visiting_Iterator(Store);

  begin
    for i in 1..nq loop                   -- check if binomial system
      if Number_of_Terms(p(i)) /= 2 then
        fail := true;
        return;
      end if;
    end loop;
    fail := false;
    for i in 1..nq loop
      equ := i; mon := 1;
      Store_Terms(p(i));
    end loop;
  end Parse;

  procedure Parse ( p : in Laur_Sys; nq : in integer32;
                    A : out Standard_Integer64_Matrices.Matrix;
                    c : out Vector; fail : out boolean ) is

    use DoblDobl_Complex_Laurentials;

    equ,mon : integer32;

    procedure Store ( t : in Term; continue : out boolean ) is
    begin
      if mon = 1 then
        for i in t.dg'range loop     -- store monomial in column equ
          A(i,equ) := integer64(t.dg(i));
        end loop;
        c(equ) := t.cf;              -- store coefficient at equ
        mon := mon + 1;
      else
        for i in t.dg'range loop     -- divide by new monomial
          A(i,equ) := A(i,equ) - integer64(t.dg(i));
        end loop;
        c(equ) := -t.cf/c(equ);      -- coefficient to righthand side
      end if;
      continue := true;
    end Store;
    procedure Store_Terms is new Visiting_Iterator(Store);

  begin
    for i in 1..nq loop                   -- check if binomial system
      if Number_of_Terms(p(i)) /= 2 then
        fail := true;
        return;
      end if;
    end loop;
    fail := false;
    for i in 1..nq loop
      equ := i; mon := 1;
      Store_Terms(p(i));
    end loop;
  end Parse;

  function Create ( A : Standard_Integer64_Matrices.Matrix;
                    c : Vector ) return Poly_Sys is

    use DoblDobl_Complex_Polynomials;

    res : Poly_Sys(A'range(2));
    t1,t2 : Term;
    one : constant double_double := Double_Double_Numbers.create(1.0);

  begin
    t1.cf := DoblDobl_Complex_Numbers.create(one);
    t1.dg := new Standard_Natural_Vectors.Vector(A'range(1));
    t2.dg := new Standard_Natural_Vectors.Vector(A'range(1));
    for j in res'range loop
      for i in A'range(1) loop
        if A(i,j) >= 0
         then t1.dg(i) := natural32(A(i,j)); t2.dg(i) := 0;
         else t1.dg(i) := 0; t2.dg(i) := natural32(-A(i,j));
        end if;
      end loop;
      res(j) := Create(t1);
      t2.cf := -c(j);
      Add(res(j),t2);
    end loop;
    Clear(t1); Clear(t2);
    return res;
  end Create;

  function Create ( A : Standard_Integer64_Matrices.Matrix;
                    c : Vector ) return Laur_Sys is

    use DoblDobl_Complex_Laurentials;

    res : Laur_Sys(A'range(2));
    t1,t2 : Term;
    one : constant double_double := Double_Double_Numbers.Create(1.0);

  begin
    t1.cf := DoblDobl_Complex_Numbers.Create(one);
    t1.dg := new Standard_Integer_Vectors.Vector(A'range(1));
    t2.dg := new Standard_Integer_Vectors.Vector(A'range(1));
    for j in res'range loop
      for i in A'range(1) loop
        t1.dg(i) := integer32(A(i,j));
        t2.dg(i) := 0;
      end loop;
      res(j) := Create(t1);
      t2.cf := -c(j);
      Add(res(j),t2);
    end loop;
    Clear(t1); Clear(t2);
    return res;
  end Create;

-- EVALUATION of a BINOMIAL SYSTEM :

  function Eval ( A : Standard_Integer64_Matrices.Matrix;
                  x : Vector ) return Vector is

    res : Vector(A'range(2));
    one : constant double_double := Double_Double_Numbers.create(1.0);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Complex_Numbers.create(one);
    end loop;
    for j in A'range(2) loop
      for i in A'range(1) loop
        res(j) := res(j)*(x(i)**integer(A(i,j)));
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( A : Multprec_Integer_Matrices.Matrix; x : Vector )
                return Vector is

    res : Vector(A'range(2));
    pwr : Complex_Number;
    one : constant double_double := create(1.0);

  begin
    for i in res'range loop
      res(i) := create(one);
    end loop;
    for j in A'range(2) loop
      for i in A'range(1) loop
        pwr := Polar_Exponentiation_ModTwoPi_of_Unit(x(i),A(i,j));
        res(j) := res(j)*pwr; -- res(j)*(x(i)**A(i,j));
      end loop;
    end loop;
    return res;
  end Eval;

  procedure Eval ( A : in Standard_Integer64_Matrices.Matrix;
                   x : in Vector; y : out Vector ) is

    one : constant double_double := create(1.0);

  begin
    y := (y'range => Create(one));
    for j in A'range(2) loop
      for i in A'range(1) loop
        y(j) := y(j)*(x(i)**integer(A(i,j)));
      end loop;
    end loop;
  end Eval;

  function Eval ( A : Standard_Integer64_Matrices.Matrix;
                  c,x : Vector ) return Vector is

    res : constant Vector(A'range(2)) := Eval(A,x) - c;

  begin
    return res;
  end Eval;

  function Eval ( A : Standard_Integer64_Matrices.Matrix;
                  s : Solution_List ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        res_ls : Solution(ls.n);
      begin
        res_ls.t := ls.t;
        res_ls.m := ls.m;
        res_ls.v := Eval(A,ls.v);
        res_ls.err := ls.err;
        res_ls.rco := ls.rco;
        res_ls.res := ls.res;
        Append(res,res_last,res_ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Eval;

  function Eval ( A : Multprec_Integer_Matrices.Matrix;
                  s : Solution_List ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        res_ls : Solution(ls.n);
      begin
        res_ls.t := ls.t;
        res_ls.m := ls.m;
        res_ls.v := Eval(A,ls.v);
        res_ls.err := ls.err;
        res_ls.rco := ls.rco;
        res_ls.res := ls.res;
        Append(res,res_last,res_ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Eval;

  procedure Eval ( A : in Standard_Integer64_Matrices.Matrix;
                   s : in Solution_List; w : in out Vector ) is

    tmp : Solution_List := s;
    ls : Link_to_Solution;
    one : constant double_double := create(1.0);
    wj : Complex_Number;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      for j in w'range loop
        w(j) := Create(one);
      end loop;
      for j in A'range(2) loop
        wj := w(j);
        for i in A'range(1) loop
         --  wj := wj*(ls.v(i)**(A(i,j)));
          if A(i,j) > 0 then
            for k in 1..A(i,j) loop
              wj := wj*ls.v(i);
            end loop;
          elsif A(i,j) < 0 then
            for k in 1..(-A(i,j)) loop
              wj := wj/ls.v(i);
            end loop;
          end if;
        end loop;
        w(j) := wj;
      end loop;
      for i in ls.v'range loop
        ls.v(i) := w(i);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
  end Eval;

end DoblDobl_Binomial_Systems;
