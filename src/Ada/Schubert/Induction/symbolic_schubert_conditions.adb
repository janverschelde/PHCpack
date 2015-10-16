with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Matrix_Indeterminates;
with Bracket_Monomials;                 use Bracket_Monomials;
with Symbolic_Minor_Equations;          use Symbolic_Minor_Equations;

package body Symbolic_Schubert_Conditions is

  function General_Localization_Map
             ( n,k : in integer32 ) return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(1..n,1..k);

  begin
    for i in 1..n loop
      for j in 1..k loop
        if i = j then
          res(i,j) := 1;
        elsif (i < j) or (i > n - k + j) then
          res(i,j) := 0;
        else
          res(i,j) := 2;
        end if;
      end loop;
    end loop;
    return res;
  end General_Localization_Map;

  function Symbolic_Form_of_Plane
             ( n,k : integer32; locmap : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Matrix_Indeterminates.Dimension(locmap));
    ind : integer32 := 0;
    t : Term;

  begin
    for i in 1..n loop
      for j in 1..k loop
        if locmap(i,j) = 0 then
          res(i,j) := Null_Poly;
        elsif locmap(i,j) = 1 then
          t.cf := Create(integer(1));
          t.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
          res(i,j) := Create(t);
          Clear(t);
        else
          ind := ind + 1;
          t.cf := Create(integer(1));
          t.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
          t.dg(ind) := 1;
          res(i,j) := Create(t);
          Clear(t);
        end if;
      end loop;
    end loop;
    return res;
  end Symbolic_Form_of_Plane;

  function Symbolic_Form_of_Plane
             ( n,k : integer32; locmap : Standard_Natural_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Matrix_Indeterminates.Dimension(locmap));
    ind : integer32 := 0;
    t : Term;

  begin
    for i in 1..n loop
      for j in 1..k loop
        if locmap(i,j) = 0 then
          res(i,j) := Null_Poly;
        elsif locmap(i,j) = 1 then
          t.cf := Create(integer(1));
          t.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
          res(i,j) := Create(t);
          Clear(t);
        else
          ind := ind + 1;
          t.cf := Create(integer(1));
          t.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
          t.dg(ind) := 1;
          res(i,j) := Create(t);
          Clear(t);
        end if;
      end loop;
    end loop;
    return res;
  end Symbolic_Form_of_Plane;

  function Symbolic_Form_of_Plane
             ( n,k : integer32; locmap : Standard_Natural_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Matrix_Indeterminates.Dimension(locmap));
    ind : integer32 := 0;
    t : Term;

  begin
    for i in 1..n loop
      for j in 1..k loop
        if locmap(i,j) = 0 then
          res(i,j) := Null_Poly;
        elsif locmap(i,j) = 1 then
          t.cf := Create(integer(1));
          t.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
          res(i,j) := Create(t);
          Clear(t);
        else
          ind := ind + 1;
          t.cf := Create(integer(1));
          t.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
          t.dg(ind) := 1;
          res(i,j) := Create(t);
          Clear(t);
        end if;
      end loop;
    end loop;
    return res;
  end Symbolic_Form_of_Plane;

  function Number_of_Equations ( n,k,f,i : natural32 ) return natural32 is

    m : constant natural32 := k + f;  -- number of columns
    r : constant natural32 := m - i;  -- rank condition

  begin
    if ((r+1 <= n) and (r+1 <= m)) then
      if ((m <= n) and (r+1 = m)) then
        return Number_of_Maximal_Minors(n,m);
      else
        return Number_of_Minors(n,m,r+1);
      end if;
    else
      return 0;                     -- trivial condition
    end if;
  end Number_of_Equations;

  function Number_of_Equations
             ( n : natural32; b : Bracket ) return natural32 is 

    res : natural32 := 0;
    k : constant natural32 := natural32(b'last);

  begin
    for i in b'range loop
      res := res + Number_of_Equations(n,k,b(i),natural32(i));
    end loop;
    return res;
  end Number_of_Equations;

  procedure Enumerate_NotAbove ( n : in natural32; b : in Bracket ) is

    k : constant integer32 := b'last;
    accu : Bracket(1..k);
    start : natural32;
    cont : boolean := true;

    procedure Enumerate ( i : in integer32 ) is
    begin
      if i > k then
        if not (accu <= b)
         then Process(accu,cont);
        end if;
      else
        if i = 1
         then start := 1;
         else start := accu(i-1) + 1;
        end if;
        for j in start..n loop
          accu(i) := j;
          Enumerate(i+1);
          exit when not cont;
        end loop;
      end if;
    end Enumerate;
 
  begin
    Enumerate(1);
  end Enumerate_NotAbove;

  function Number_of_NotAbove
             ( n : natural32; b : Bracket ) return natural32 is

    res : natural32 := 0;

    procedure Count ( b : in Bracket; cont : out boolean ) is
    begin
      res := res + 1;
      cont := true;
    end Count;
    procedure Enum is new Enumerate_NotAbove(Count);

  begin
    Enum(n,b);
    return res;
  end Number_of_NotAbove;

  procedure Explain_Equations
               ( n : in natural32; b : in bracket; nq : out natural32 ) is

    k : constant integer32 := b'last;
    m,r,d : natural32;
    cnt : natural32 := 0;

  begin
    put("Our "); put(k,1); put("-plane X in "); put(n,1);
    put_line("-space is subject to the following conditions : ");
    for i in 1..k loop
      m := natural32(k) + b(i);
      r := m - natural32(i);
      put("  X meets F("); put(b(i),1); put(") in a "); put(i,1); 
      put("-plane : Rank([ X | F("); put(b(i),1); put(") ]) = ");
      put(r,1); new_line; 
      put("  => all "); put(r+1,1); put("-by-"); put(r+1,1);
      put(" minors of a "); put(n,1); put("-by-"); put(m,1);
      put_line(" matrix must be zero");
      if ((r+1 > n) or (r+1 > m)) then
        put_line("  trivial condition, no minor equations");
      else
        if ((m <= n) and (r+1 = m)) then
          d := Number_of_Maximal_Minors(n,m);
        else
          d := Number_of_Minors(n,m,r+1);
        end if;
        put("  add "); put(d,1); put_line(" minor equations");
        cnt := cnt + d;
      end if;
    end loop;
    put("The Schubert conditions consist of ");
    put(cnt,1); put_line(" minor equations.");
    nq := cnt;
  end Explain_Equations;

  function Flag_Minors ( n,k,f,i : natural32 ) return Bracket_Polynomial is

    res : Bracket_Polynomial := Null_Bracket_Poly;
    m : constant natural32 := k + f;
    r : constant natural32 := m - i;

  begin
    if ((r+1 <= n) and (r+1 <= m))
     then res := Minors(n,m,r+1);
    end if;
    return res;
  end Flag_Minors;

  function Flag_Minors ( n : natural32; b : bracket ) return Bracket_System is

    res : Bracket_System(b'range);
    k : constant natural32 := natural32(b'last);

  begin
    for i in b'range loop
      res(i) := Flag_Minors(n,k,b(i),natural32(i));
    end loop;
    return res;
  end Flag_Minors;

  procedure Flag_Minor_Polynomials
              ( p : in Bracket_Polynomial; bs : in out Bracket_System;
                ind : in out integer32 ) is

    procedure Store_Term ( t : in Bracket_Term; cont : out boolean ) is

      row : Bracket_Monomial;
      first : boolean := true;

      procedure Store_Bracket ( b : in Bracket; c : out boolean ) is

        t : Bracket_Term;

      begin
        if first then
          first := false;
          row := Create(b);
        else
          t.coeff := Standard_Complex_Numbers.Create(integer(1));
          Copy_Append(row,t.monom);
          Append(t.monom,b);
          ind := ind + 1;
          bs(ind) := Create(t);
        end if;
        c := true;
      end Store_Bracket;
      procedure Store_Brackets is new Enumerate_Brackets(Store_Bracket);

    begin
      if Number_of_Brackets(t.monom) < 3 then
        ind := ind + 1;
        bs(ind) := Create(t);
      else
        Store_Brackets(t.monom);
      end if;
      cont := true;
    end Store_Term;
    procedure Store_Terms is new Enumerate_Terms(Store_Term);

  begin
    Store_Terms(p);
  end Flag_Minor_Polynomials;

  function Flag_Minor_System
              ( m : natural32; fm : Bracket_Polynomial )
              return Bracket_System is

    res : Bracket_System(1..integer32(m));
    ind : integer32 := 0;

  begin
    Flag_Minor_Polynomials(fm,res,ind);
    return res;
  end Flag_Minor_System;

  function Flag_Minor_System
              ( m : natural32; bs : Bracket_System ) return Bracket_System is

    res : Bracket_System(1..integer32(m));
    ind : integer32 := 0;

  begin
    for i in bs'range loop
      if bs(i) /= Null_Bracket_Poly
       then Flag_Minor_Polynomials(bs(i),res,ind);
      end if;
    end loop;
    return res;
  end Flag_Minor_System;

end Symbolic_Schubert_Conditions;
