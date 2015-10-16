with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;

package body Matrix_Indeterminates is

  procedure Initialize_Symbols ( n,d : in natural32 ) is
  begin
    Symbol_Table.Init(n*d);
    for i in 1..n loop
      for j in 1..d loop
        declare
          s : constant Symbol := X_ij(i,j);
        begin
          Symbol_Table.Add(s);
        end;
      end loop;
    end loop;
  end Initialize_Symbols;

  function Dimension ( locmap : Matrix ) return natural32 is

    res : natural32 := 0;

  begin
    for i in locmap'range(1) loop
      for j in locmap'range(2) loop
        if locmap(i,j) = 2
         then res := res + 1;
        end if;
      end loop;
    end loop;
    return res;
  end Dimension;

  procedure Initialize_Symbols
               ( d : in natural32; locmap : in Matrix ) is
  begin
    Symbol_Table.Init(d);
    for i in locmap'range(1) loop
      for j in locmap'range(2) loop
        if locmap(i,j) = 2 then
          declare
            s : constant Symbol := X_ij(natural32(i),natural32(j));
          begin
            Symbol_Table.Add(s);
          end;
        end if;
      end loop;
    end loop;
  end Initialize_Symbols;

  procedure Initialize_Symbols ( locmap : in Matrix ) is

    d : constant natural32 := Dimension(locmap);

  begin
    Initialize_Symbols(d,locmap);
  end Initialize_Symbols;

  function X_ij ( i,j : natural32 ) return Symbol is

    res : Symbol;   

  begin
    res(1) := 'x';
    res(2) := Convert_Hexadecimal(i);
    res(3) := Convert_Hexadecimal(j);
    for i in 4..res'last loop
      res(i) := ' ';
    end loop;
    return res;
  end X_ij;

  function Monomial ( n,d,i,j : natural32 )
                    return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.cf := Create(1.0);
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n*d) => 0);
    t.dg(integer32((i-1)*d+j)) := 1;
    res := Create(t);
    Clear(t.dg);
    return res;
  end Monomial;

  function Monomial ( n,d,i,j : natural32 )
                    return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : Poly;
    one : constant double_double := create(1.0);
    t : Term;

  begin
    t.cf := Create(one);
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n*d) => 0);
    t.dg(integer32((i-1)*d+j)) := 1;
    res := Create(t);
    Clear(t.dg);
    return res;
  end Monomial;

  function Monomial ( n,d,i,j : natural32 )
                    return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : Poly;
    one : constant quad_double := create(1.0);
    t : Term;

  begin
    t.cf := Create(one);
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n*d) => 0);
    t.dg(integer32((i-1)*d+j)) := 1;
    res := Create(t);
    Clear(t.dg);
    return res;
  end Monomial;

  procedure Reduce_Symbols ( locmap : in Matrix ) is
  begin
    for i in reverse locmap'range(1) loop
      for j in reverse locmap'range(2) loop
        if locmap(i,j) /= 2
         then Symbol_Table.Remove(X_ij(natural32(i),natural32(j)));
        end if;
      end loop;
    end loop;
  end Reduce_Symbols;

  procedure Clear_Symbols is
  begin
    Symbol_Table.Clear;
  end Clear_Symbols;

end Matrix_Indeterminates;
