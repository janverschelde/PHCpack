with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Matrix_Indeterminates;
with Straightening_Syzygies;             use Straightening_Syzygies;

package body Symbolic_Minor_Equations is

-- AUXILIARIES TO Minor_Equations :

  function Substitute ( b,minor : Bracket ) return Bracket is

  -- DESCRIPTION :
  --   The ith entry in the bracket b is replaced by minor(b(i)).

    res : Bracket(b'range);
    start : integer32;

  begin
    if b(b'first) = 0
     then start := b'first+1;         -- bracket is a coefficient bracket
          res(res'first) := 0;        -- result is also a coefficient bracket
     else start := b'first;           -- second bracket in expansion
    end if;
    for i in start..b'last loop
      res(i) := minor(integer32(b(i)));
    end loop;
    return res;
  end Substitute;

  function Substitute ( bm : Bracket_Monomial; minor : Bracket )
                      return Bracket_Monomial is

  -- DESCRIPTION :
  --   Substitutes the entries in the brackets according to the minor.

  -- REQUIRED : bm is a quadratic monomial.

    res : Bracket_Monomial;
    lb : Link_to_Bracket;
    first : boolean := true;

    procedure Substitute_Bracket ( b : in Bracket; continue : out boolean ) is

      sb : constant Bracket(b'range) := Substitute(b,minor);

    begin
      if first then
        lb := new Bracket'(sb);  -- save to preserve order
        first := false;
      else
        res := sb*lb.all;
      end if;
      continue := true;
    end Substitute_Bracket;
    procedure Substitute_Brackets is
      new Enumerate_Brackets(Substitute_Bracket);

  begin
    Substitute_Brackets(bm);
    Clear(lb);
    return res;
  end Substitute;

  function Substitute ( bt : Bracket_Term; minor : Bracket )
                      return Bracket_Term is

  -- DESCRIPTION :
  --   The ith entry in every bracket of the term according to the minor.

    res : Bracket_Term;

  begin
    res.coeff := bt.coeff;
    res.monom := Substitute(bt.monom,minor);
    return res;
  end Substitute;

  function Substitute ( bp : Bracket_Polynomial; minor : Bracket )
                      return Bracket_Polynomial is

  -- DESCRIPTION :
  --   The labels in the brackets of bp are replaced by those in the minor.
  --   The bracket polynomial represents the Laplace expansion of that minor.

    res : Bracket_Polynomial;

    procedure Substitute_Term ( t : in Bracket_Term; continue : out boolean ) is

      st : Bracket_Term := Substitute(t,minor);

    begin
      Add(res,st);
      Clear(st);
      continue := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Enumerate_Terms(Substitute_Term);

  begin
    Substitute_Terms(bp);
    return res;
  end Substitute;

-- AUXILIARIES TO Expanded_Minor :

  function Subtract ( b : Bracket; i : integer32 ) return Bracket is

  -- DESCRIPTION :
  --   Returns a smaller bracket with the ith entry removed.

    res : Bracket(b'first..b'last-1);

  begin
    res(b'first..(i-1)) := b(b'first..(i-1));
    res(i..res'last) := b((i+1)..b'last);
    return res;
  end Subtract;

  function General_Expanded_Minor
             ( m : Standard_Complex_Poly_Matrices.Matrix;
               b : Bracket )
             return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   This function treats the case for Expanded_Minor when b'length > 2.

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    sig : integer32;
    submin,acc : Poly;

  begin
    if b'last mod 2 = 0
     then sig := -1;
     else sig := +1;
    end if;
    for i in b'range loop
      if m(integer32(b(i)),b'last) /= Null_Poly then
        submin := Expanded_Minor(m,Subtract(b,i));
        if submin /= Null_Poly then
          acc := m(integer32(b(i)),b'last)*submin;
          if sig > 0
           then Add(res,acc);
           else Sub(res,acc);
          end if;
          Clear(acc);
        end if;
        Clear(submin);
      end if;
      sig := -sig;
    end loop;
    return res;
  end General_Expanded_Minor;

  function General_Expanded_Minor
             ( m : DoblDobl_Complex_Poly_Matrices.Matrix;
               b : Bracket )
             return DoblDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   This function treats the case for Expanded_Minor when b'length > 2.

    use DoblDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    sig : integer32;
    submin,acc : Poly;

  begin
    if b'last mod 2 = 0
     then sig := -1;
     else sig := +1;
    end if;
    for i in b'range loop
      if m(integer32(b(i)),b'last) /= Null_Poly then
        submin := Expanded_Minor(m,Subtract(b,i));
        if submin /= Null_Poly then
          acc := m(integer32(b(i)),b'last)*submin;
          if sig > 0
           then Add(res,acc);
           else Sub(res,acc);
          end if;
          Clear(acc);
        end if;
        Clear(submin);
      end if;
      sig := -sig;
    end loop;
    return res;
  end General_Expanded_Minor;

  function General_Expanded_Minor
             ( m : QuadDobl_Complex_Poly_Matrices.Matrix;
               b : Bracket )
             return QuadDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   This function treats the case for Expanded_Minor when b'length > 2.

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    sig : integer32;
    submin,acc : Poly;

  begin
    if b'last mod 2 = 0
     then sig := -1;
     else sig := +1;
    end if;
    for i in b'range loop
      if m(integer32(b(i)),b'last) /= Null_Poly then
        submin := Expanded_Minor(m,Subtract(b,i));
        if submin /= Null_Poly then
          acc := m(integer32(b(i)),b'last)*submin;
          if sig > 0
           then Add(res,acc);
           else Sub(res,acc);
          end if;
          Clear(acc);
        end if;
        Clear(submin);
      end if;
      sig := -sig;
    end loop;
    return res;
  end General_Expanded_Minor;

  procedure Purify ( p : in out Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Eliminates terms that have coefficients less than 10.0E-10.

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    tol : constant double_float := 10.0E-10;
    res : Poly := Null_Poly;

    procedure Scan_Term ( t : in Term; cont : out boolean ) is
    begin
      if AbsVal(t.cf) > tol
       then Add(res,t);
      end if;
      cont := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    Clear(p);
    p := res;
  end Purify;

  procedure Purify ( p : in out DoblDobl_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Eliminates terms that have coefficients less than 10.0E-10.
  
    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    tol : constant double_float := 10.0E-10;
    res : Poly := Null_Poly;

    procedure Scan_Term ( t : in Term; cont : out boolean ) is
    begin
      if AbsVal(t.cf) > tol
       then Add(res,t);
      end if;
      cont := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    Clear(p);
    p := res;
  end Purify;

  procedure Purify ( p : in out QuadDobl_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Eliminates terms that have coefficients less than 10.0E-10.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    tol : constant double_float := 10.0E-10;
    res : Poly := Null_Poly;

    procedure Scan_Term ( t : in Term; cont : out boolean ) is
    begin
      if AbsVal(t.cf) > tol
       then Add(res,t);
      end if;
      cont := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    Clear(p);
    p := res;
  end Purify;

-- TARGET ROUTINES :

  function Schubert_Pattern
             ( n : natural32; b1,b2 : Bracket )
             return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(1..integer32(n),b1'range);
    d : constant natural32 := natural32(b1'last);

    use Standard_Complex_Polynomials;

  begin
    for i in res'range(1) loop
      for j in b1'range loop
        if ((natural32(i) < b2(j)) or (natural32(i) > n+1-b1(b1'last+1-j)))
         then res(i,j) := Null_Poly;
         else res(i,j) := Matrix_Indeterminates.Monomial
                            (n,d,natural32(i),natural32(j));
        end if;
      end loop;
    end loop;
    return res;
  end Schubert_Pattern;

  function Schubert_Pattern
             ( n : natural32; b1,b2 : Bracket )
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    res : DoblDobl_Complex_Poly_Matrices.Matrix(1..integer32(n),b1'range);
    d : constant natural32 := natural32(b1'last);

    use DoblDobl_Complex_Polynomials;

  begin
    for i in res'range(1) loop
      for j in b1'range loop
        if ((natural32(i) < b2(j)) or (natural32(i) > n+1-b1(b1'last+1-j)))
         then res(i,j) := Null_Poly;
         else res(i,j) := Matrix_Indeterminates.Monomial
                            (n,d,natural32(i),natural32(j));
        end if;
      end loop;
    end loop;
    return res;
  end Schubert_Pattern;

  function Schubert_Pattern
             ( n : natural32; b1,b2 : Bracket )
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    res : QuadDobl_Complex_Poly_Matrices.Matrix(1..integer32(n),b1'range);
    d : constant natural32 := natural32(b1'last);

    use QuadDobl_Complex_Polynomials;

  begin
    for i in res'range(1) loop
      for j in b1'range loop
        if ((natural32(i) < b2(j)) or (natural32(i) > n+1-b1(b1'last+1-j)))
         then res(i,j) := Null_Poly;
         else res(i,j) := Matrix_Indeterminates.Monomial
                            (n,d,natural32(i),natural32(j));
        end if;
      end loop;
    end loop;
    return res;
  end Schubert_Pattern;

  function Localization_Pattern
             ( n : natural32; top,bottom : Bracket )
             return Standard_Complex_Poly_Matrices.Matrix is

    p : constant integer32 := top'length;
    res : Standard_Complex_Poly_Matrices.Matrix(1..integer32(n),1..p);

    use Standard_Complex_Polynomials;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if natural32(i) < top(j) or natural32(i) > bottom(j)
         then res(i,j) := Null_Poly;
         else res(i,j) := Matrix_Indeterminates.Monomial
                            (n,natural32(p),natural32(i),natural32(j));
        end if;
      end loop;
    end loop;
    return res;
  end Localization_Pattern;

  function Localization_Pattern
             ( n : natural32; top,bottom : Bracket )
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    p : constant integer32 := top'length;
    res : DoblDobl_Complex_Poly_Matrices.Matrix(1..integer32(n),1..p);

    use DoblDobl_Complex_Polynomials;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if natural32(i) < top(j) or natural32(i) > bottom(j)
         then res(i,j) := Null_Poly;
         else res(i,j) := Matrix_Indeterminates.Monomial
                            (n,natural32(p),natural32(i),natural32(j));
        end if;
      end loop;
    end loop;
    return res;
  end Localization_Pattern;

  function Localization_Pattern
             ( n : natural32; top,bottom : Bracket )
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    p : constant integer32 := top'length;
    res : QuadDobl_Complex_Poly_Matrices.Matrix(1..integer32(n),1..p);

    use QuadDobl_Complex_Polynomials;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if natural32(i) < top(j) or natural32(i) > bottom(j)
         then res(i,j) := Null_Poly;
         else res(i,j) := Matrix_Indeterminates.Monomial
                            (n,natural32(p),natural32(i),natural32(j));
        end if;
      end loop;
    end loop;
    return res;
  end Localization_Pattern;

-- SYMBOLIC REPRESENTATION OF THE EQUATIONS :

  function Number_of_Choices ( n,k : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns n!/(k!*(n-k)!), or the number of choices of k out of n.

    res : natural32 := 1;

  begin
    for i in k+1..n loop    -- compute n!/k!
      res := res*i;
    end loop;
    for i in 2..n-k loop    -- divide by (n-k)!
      res := res/i;
    end loop;
    return res;
  end Number_of_Choices;

  function Number_of_Maximal_Minors ( n,m : natural32 ) return natural32 is
  begin
    return Number_of_Choices(n,m);
  end Number_of_Maximal_Minors;

  function Number_of_Minors ( n,m,s : natural32 ) return natural32 is
  begin
    return Number_of_Choices(n,s)*Number_of_Choices(m,s);
  end Number_of_Minors;

  function Maximal_Minors ( n,m : natural32 ) return Bracket_Monomial is

    res : Bracket_Monomial;
    first : boolean := true;
    accu : Bracket(1..integer32(m));

    procedure Enumerate_Brackets ( k,start : integer32 ) is

    -- DESCRIPTION :
    --   Enumerate all m-brackets with entries between 1 and n,
    --   beginning at start, and currently filling in the kth entry.

    begin
      if k > integer32(m) then
        if first
         then res := Create(accu); first := false;
         else Multiply(res,accu);
        end if;
      else
        for i in start..integer32(n) loop
          accu(k) := natural32(i);
          Enumerate_Brackets(k+1,i+1);
        end loop;
      end if;
    end Enumerate_Brackets;

  begin
    Enumerate_Brackets(1,1);
    return res;
  end Maximal_Minors;

  function Minors ( n,m,s : natural32 ) return Bracket_Polynomial is

    res : Bracket_Polynomial;
    accu_row,accu_col : Bracket(1..integer32(s));

    procedure Enumerate_Rows ( row,start_row : integer32 ) is

    -- DESCRIPTION :
    --   Enumerates all s-brackets with entries between 1 and n,
    --   beginning at start_row, and currently filling in at row.

      bt : Bracket_Term;

      procedure Enumerate_Columns ( col,start_col : integer32 ) is

      -- DESCRIPTION :
      --   Enumerates all s-brackets with entries between 1 and m,
      --   beginning at start_col, and currently filling in at col.

      begin
        if col > integer32(s) then
          Append(bt.monom,accu_col);
        else
          for i in start_col..integer32(m) loop
            accu_col(col) := natural32(i);
            Enumerate_Columns(col+1,i+1);
          end loop;
        end if;  
      end Enumerate_Columns;

    begin
      if row > integer32(s) then
        bt.coeff := Standard_Complex_Numbers.Create(1.0);
        bt.monom := Create(accu_row);
        Enumerate_Columns(1,1);
        Frontal_Construct(res,bt); Clear(bt);
      else
        for i in start_row..integer32(n) loop
          accu_row(row) := natural32(i);
          Enumerate_Rows(row+1,i+1);
        end loop;
      end if;
    end Enumerate_Rows;

  begin
    Enumerate_Rows(1,1);
    return res;
  end Minors;

  function Minor_Equations
             ( m,d : natural32; bm : Bracket_Monomial )
             return Bracket_System is

    res : Bracket_System(0..integer32(Number_of_Brackets(bm)));
    gen : constant Bracket_Polynomial := Laplace_Expansion(m,d);
    cnt : integer32 := 0;

    procedure Substitute_Bracket ( b : in Bracket; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      res(cnt) := Substitute(gen,b);
      continue := true;
    end Substitute_Bracket;
    procedure Substitute_Brackets is
      new Enumerate_Brackets(Substitute_Bracket);

  begin
    res(0) := gen;
    Substitute_Brackets(bm);
    return res;
  end Minor_Equations;

  function Expanded_Minor
             ( m : Standard_Complex_Poly_Matrices.Matrix; b : Bracket )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    acc : Poly;

  begin
    if b'length = 1 then
      Copy(m(integer32(b(1)),1),res);
    elsif b'length = 2 then
      if (m(integer32(b(1)),1) /= Null_Poly)
         and (m(integer32(b(2)),2) /= Null_Poly)
       then res := m(integer32(b(1)),1)*m(integer32(b(2)),2);
      end if;
      if (m(integer32(b(2)),1) /= Null_Poly)
          and (m(integer32(b(1)),2) /= Null_Poly) then
        acc := m(integer32(b(2)),1)*m(integer32(b(1)),2);
        Sub(res,acc);
        Clear(acc);
      end if;
    else
      res := General_Expanded_Minor(m,b);
    end if;
    Purify(res);
    return res;
  end Expanded_Minor;

  function Expanded_Minor
             ( m : DoblDobl_Complex_Poly_Matrices.Matrix; b : Bracket )
             return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    acc : Poly;

  begin
    if b'length = 1 then
      Copy(m(integer32(b(1)),1),res);
    elsif b'length = 2 then
      if (m(integer32(b(1)),1) /= Null_Poly)
         and (m(integer32(b(2)),2) /= Null_Poly)
       then res := m(integer32(b(1)),1)*m(integer32(b(2)),2);
      end if;
      if (m(integer32(b(2)),1) /= Null_Poly)
          and (m(integer32(b(1)),2) /= Null_Poly) then
        acc := m(integer32(b(2)),1)*m(integer32(b(1)),2);
        Sub(res,acc);
        Clear(acc);
      end if;
    else
      res := General_Expanded_Minor(m,b);
    end if;
    Purify(res);
    return res;
  end Expanded_Minor;

  function Expanded_Minor
             ( m : QuadDobl_Complex_Poly_Matrices.Matrix; b : Bracket )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    acc : Poly;

  begin
    if b'length = 1 then
      Copy(m(integer32(b(1)),1),res);
    elsif b'length = 2 then
      if (m(integer32(b(1)),1) /= Null_Poly)
         and (m(integer32(b(2)),2) /= Null_Poly)
       then res := m(integer32(b(1)),1)*m(integer32(b(2)),2);
      end if;
      if (m(integer32(b(2)),1) /= Null_Poly)
          and (m(integer32(b(1)),2) /= Null_Poly) then
        acc := m(integer32(b(2)),1)*m(integer32(b(1)),2);
        Sub(res,acc);
        Clear(acc);
      end if;
    else
      res := General_Expanded_Minor(m,b);
    end if;
    Purify(res);
    return res;
  end Expanded_Minor;

  function Extend_Zero_Lifting
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Extend_Term ( t : in Term; continue : out boolean ) is

      et : Term;
    
    begin
      et.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last+1);
      et.dg(t.dg'range) := t.dg.all;
      et.dg(et.dg'last) := 0;
      et.cf := t.cf;
      Add(res,et);
      Clear(et.dg);
      continue := true;
    end Extend_Term;
    procedure Extend_Terms is new Visiting_Iterator(Extend_Term);

  begin
    Extend_Terms(p);
    return res;
  end Extend_Zero_Lifting;

end Symbolic_Minor_Equations;
