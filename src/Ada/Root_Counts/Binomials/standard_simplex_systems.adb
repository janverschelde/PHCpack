with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Complex_Numbers;         use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Laurentials;
with Standard_Complex_Polynomials;
with Standard_Binomial_Systems;        use Standard_Binomial_Systems;

package body Standard_Simplex_Systems is

-- UTILITIES :

  function Row_Minimum ( A : Standard_Integer_Matrices.Matrix;
                         i : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the minimum value at row i of A.

  -- REQUIRED : i is in A'range(1).

    min : integer32 := A(i,A'first(2));

  begin
    for j in A'first(2)+1..A'last(2) loop
      if A(i,j) < min
       then min := A(i,j);
      end if;
    end loop;
    return min;
  end Row_Minimum;

  function Is_In ( A : Standard_Integer_Matrices.Matrix; k : integer32;
                   v : Standard_Integer_Vectors.Vector ) return integer32 is

  -- DESCRIPTION :
  --   Returns 0 if d does not belong among the first k columns of A,
  --   otherwise the corresponding column index j of A is returned for
  --   which d equals that jth column of A.

    res : integer32;

  begin
    for j in 1..k loop
      res := j;
      for i in A'range(1) loop
        if A(i,j) /= v(i)
         then res := 0; exit;
        end if;
      end loop;
      if res > 0
       then return res;
      end if;
    end loop;
    return 0;
  end Is_In;

  function Least ( p : Standard_Complex_Polynomials.Poly )
                 return Standard_Complex_Polynomials.Term is

  -- DESCRIPTION :
  --   Returns the last term of the polynomial p,
  --   which is least in the graded lexicographic ordering.

    use Standard_Complex_Polynomials;

    res : Term;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      res := t;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Least;

  function Least ( p : Standard_Complex_Laurentials.Poly )
                 return Standard_Complex_Laurentials.Term is

  -- DESCRIPTION :
  --   Returns the last term of the polynomial p,
  --   which is least in the graded lexicographic ordering.

    use Standard_Complex_Laurentials;

    res : Term;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      res := t;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Least;

  procedure Parse ( p : in Standard_Complex_Polynomials.Poly;
                    row : in integer32;
                    shift : in Standard_Natural_Vectors.Vector;
		    cnt : in out integer32;
                    A : in out Standard_Integer_Matrices.Matrix;
                    C : in out Standard_Complex_Matrices.Matrix;
                    fail : out boolean ) is

  -- DESCRIPTION :
  --   Parses the polynomial p at the given row, updating the
  --   exponent matrix and the coefficients in C.
  --   Failure is returned if cnt exceeds the column range of A.

  -- ON ENTRY :
  --   p        a polynomial in nv variables;
  --   row      current row of the system, p = p(row);
  --   shift    used to substract exponent vectors with;
  --   cnt      current number of different monomials;
  --   A        current matrix of exponents, filled with cnt columns;
  --   C        current coefficient matrix, filled with cnt columns.

  -- ON RETURN :
  --   cnt      updated counter of columns;
  --   A        updated matrix of exponents;
  --   C        updated coefficient matrix;
  --   fail     true if too many monomials.

    use Standard_Complex_Polynomials;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      wrk : Standard_Integer_Vectors.Vector(shift'range);
      col : integer32;

    begin
      if Standard_Natural_Vectors.Equal(t.dg.all,shift) then
        continue := false;           -- reached least monomial
      else
        for i in wrk'range loop
          wrk(i) := integer32(t.dg(i) - shift(i));
        end loop;
        col := Is_In(A,cnt,wrk);
        if col > 0 then              -- found a stored monomial
          C(row,col) := t.cf;
          continue := true;
        elsif cnt < A'last(2) then
          cnt := cnt + 1;            -- found new monomial
          for i in wrk'range loop    -- copy the shifted exponent vector
            A(i,cnt) := wrk(i);
          end loop;
          C(row,cnt) := t.cf;        -- update the coefficient matrix
          continue := true;
        else 
          fail := true;              -- too many monomials failure
          continue := false;
        end if;
      end if;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin 
    fail := false;
    Scan_Terms(p);
  end Parse;

  procedure Parse ( p : in Standard_Complex_Laurentials.Poly;
                    row : in integer32;
                    shift : in Standard_Integer_Vectors.Vector;
		    cnt : in out integer32;
                    A : in out Standard_Integer_Matrices.Matrix;
                    C : in out Standard_Complex_Matrices.Matrix;
                    fail : out boolean ) is

  -- DESCRIPTION :
  --   Parses the polynomial p at the given row, updating the
  --   exponent matrix and the coefficients in C.
  --   Failure is returned if cnt exceeds the column range of A.

  -- ON ENTRY :
  --   p        a polynomial in nv variables;
  --   row      current row of the system, p = p(row);
  --   shift    used to substract exponent vectors with;
  --   cnt      current number of different monomials;
  --   A        current matrix of exponents, filled with cnt columns;
  --   C        current coefficient matrix, filled with cnt columns.

  -- ON RETURN :
  --   cnt      updated counter of columns;
  --   A        updated matrix of exponents;
  --   C        updated coefficient matrix;
  --   fail     true if too many monomials.

    use Standard_Complex_Laurentials;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      wrk : Standard_Integer_Vectors.Vector(shift'range);
      col : integer32;

    begin
      if Standard_Integer_Vectors.Equal(t.dg.all,shift) then
        continue := false;           -- reached least monomial
      else
        for i in wrk'range loop
          wrk(i) := t.dg(i) - shift(i);
        end loop;
        col := Is_In(A,cnt,wrk);
        if col > 0 then              -- found a stored monomial
          C(row,col) := t.cf;
          continue := true;
        elsif cnt < A'last(2) then
          cnt := cnt + 1;            -- found new monomial
          for i in wrk'range loop    -- copy the shifted exponent vector
            A(i,cnt) := wrk(i);
          end loop;
          C(row,cnt) := t.cf;        -- update the coefficient matrix
          continue := true;
        else 
          fail := true;              -- too many monomials failure
          continue := false;
        end if;
      end if;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin 
    fail := false;
    Scan_Terms(p);
  end Parse;

-- FORMAT of a SIMPLEX SYSTEM : p(x) = 0 => C*x^A = b

  procedure Parse ( p : in Poly_Sys; nv : in integer32;
                    A : out Standard_Integer_Matrices.Matrix;
                    C : out Standard_Complex_Matrices.Matrix;
                    b : out Standard_Complex_Vectors.Vector;
                    fail : out boolean ) is

    use Standard_Complex_Polynomials;

    cnt : integer32 := 0;
    tb : Term;

  begin
    fail := false;
    for i in p'range loop
      tb := Least(p(i));
      b(i) := -tb.cf;
      for j in C'range(2) loop
        C(i,j) := Create(0.0);
      end loop;
      Parse(p(i),i,tb.dg.all,cnt,A,C,fail);
      exit when fail;
    end loop;
    if fail then
      return;
    elsif cnt < A'last(2) then           -- very few monomials
      for j in cnt+1..A'last(2) loop     -- padding with zeroes
        for i in A'range(1) loop
          A(i,j) := 0;
        end loop;
      end loop;
    end if;
  end Parse;

  procedure Parse ( p : in Laur_Sys; nv : in integer32;
                    A : out Standard_Integer_Matrices.Matrix;
                    C : out Standard_Complex_Matrices.Matrix;
                    b : out Standard_Complex_Vectors.Vector;
                    fail : out boolean ) is

    use Standard_Complex_Laurentials;

    cnt : integer32 := 0;
    tb : Term;

  begin
    fail := false;
    for i in p'range loop
      tb := Least(p(i));
      b(i) := -tb.cf;
      for j in C'range(2) loop
        C(i,j) := Create(0.0);
      end loop;
      Parse(p(i),i,tb.dg.all,cnt,A,C,fail);
      exit when fail;
    end loop;
    if fail then
      return;
    elsif cnt < A'last(2) then           -- very few monomials
      for j in cnt+1..A'last(2) loop     -- padding with zeroes
        for i in A'range(1) loop
          A(i,j) := 0;
        end loop;
      end loop;
    end if;
  end Parse;

  function Create ( A : Standard_Integer_Matrices.Matrix;
                    C : Standard_Complex_Matrices.Matrix;
                    b : Standard_Complex_Vectors.Vector ) return Poly_Sys is

    use Standard_Complex_Polynomials;

    res : Poly_Sys(C'range(1));
    ta,tb : Term;
    min : integer32;
    AA : Standard_Integer_Matrices.Matrix(A'range(1),A'range(2));

  begin
    ta.dg := new Standard_Natural_Vectors.Vector(A'range(1));
    tb.dg := new Standard_Natural_Vectors.Vector(A'range(1));
    for i in A'range(1) loop
      min := Row_Minimum(A,i);
      if min >= 0
       then tb.dg(i) := 0;
            for j in AA'range(2) loop        -- copy row i into AA
              AA(i,j) := A(i,j);
            end loop;
       else tb.dg(i) := natural32(-min);
            for j in AA'range(2) loop        -- shift row i into AA
              AA(i,j) := A(i,j) - min;       -- so all elements are >= 0
            end loop;
      end if;
    end loop;
    for i in res'range loop
      tb.cf := -b(i);
      res(i) := Create(tb);
      for j in AA'range(2) loop
        for k in AA'range(1) loop
          ta.dg(k) := natural32(AA(k,j));
        end loop;
        ta.cf := C(i,j);
        Add(res(i),ta);
      end loop;
    end loop;
    Clear(ta); Clear(tb);
    return res;
  end Create;

  function Create ( A : Standard_Integer_Matrices.Matrix;
                    C : Standard_Complex_Matrices.Matrix;
                    b : Standard_Complex_Vectors.Vector ) return Laur_Sys is

    use Standard_Complex_Laurentials;

    res : Laur_Sys(C'range(1));
    ta,tb : Term;

  begin
    ta.dg := new Standard_Integer_Vectors.Vector(A'range(1));
    tb.dg := new Standard_Integer_Vectors.Vector(A'range(1));
    for i in res'range loop
      tb.cf := -b(i);
      res(i) := Create(tb);
      for j in A'range(2) loop
        for k in A'range(1) loop
          ta.dg(k) := A(k,j);
        end loop;
        ta.cf := C(i,j);
        Add(res(i),ta);
      end loop;
    end loop;
    Clear(ta); Clear(tb);
    return res;
  end Create;

-- EVALUATION of a SIMPLEX SYSTEM :

  function Eval ( A : Standard_Integer_Matrices.Matrix;
                  C : Standard_Complex_Matrices.Matrix;
                  b,x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector is

    y : constant Standard_Complex_Vectors.Vector(A'range(2)) := Eval(A,x);
    r : Standard_Complex_Vectors.Vector(C'range(1));

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

  begin
    r := C*y;
    return r - b;
  end Eval;

end Standard_Simplex_Systems;
