with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Random_Numbers;           use Standard_Random_Numbers;

package body Standard_Random_Polynomials is

  function Random_Coefficient ( c : natural32 ) return Complex_Number is
 
    res : Complex_Number;
    ranflt : double_float;

  begin
    case c is
      when 1 => res := Create(1.0);
      when 2 => ranflt := Random;
                res := Create(ranflt);
      when others => res := Random1;
    end case;
    return res;
  end Random_Coefficient;

  function Random_Monomial
             ( n,d : natural32 ) return Standard_Complex_Polynomials.Term is

    use Standard_Complex_Polynomials;

    res : Term;
    deg,pos : integer32;

  begin
    res.cf := Create(1.0);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    for i in 1..d loop
      deg := Random(0,1);
      pos := Random(1,integer32(n));
      res.dg(pos) := res.dg(pos) + natural32(deg);
    end loop;
    return res;
  end Random_Monomial;

  function Random_Monomial
             ( n : natural32; d,e : integer32 )
             return Standard_Complex_Laurentials.Term is

    use Standard_Complex_Laurentials;

    res : Term;
    deg,pos : integer32;

  begin
    res.cf := Create(1.0);
    res.dg := new Standard_Integer_Vectors.Vector'(1..integer32(n) => 0);
    if d < 0 then
      for i in 1..(-d) loop
        deg := Random(-1,0);
        pos := Random(1,integer32(n));
        res.dg(pos) := res.dg(pos) + deg;
      end loop;
      if e < 0 then
        for i in 1..(-e) loop
          deg := Random(-1,0);
          pos := Random(1,integer32(n));
          res.dg(pos) := res.dg(pos) + deg;
        end loop;
      else
        for i in 1..e loop
          deg := Random(0,1);
          pos := Random(1,integer32(n));
          res.dg(pos) := res.dg(pos) + deg;
        end loop;
      end if;
    else -- as d >= 0 also e >= 0
      for i in 1..(d+e) loop
        deg := Random(0,1);
        pos := Random(1,integer32(n));
        res.dg(pos) := res.dg(pos) + deg;
      end loop;
    end if;
    return res;
  end Random_Monomial;

  function Random_Term
             ( n,d,c : natural32 ) return Standard_Complex_Polynomials.Term is

    use Standard_Complex_Polynomials;

    res : Term := Random_Monomial(n,d);

  begin
    res.cf := Random_Coefficient(c);
    return res;
  end Random_Term;

  function Random_Term
             ( n : natural32; d,e : integer32; c : natural32 )
             return Standard_Complex_Laurentials.Term is

    use Standard_Complex_Laurentials;

    res : Term := Random_Monomial(n,d,e);

  begin
    res.cf := Random_Coefficient(c);
    return res;
  end Random_Term;

  function Random_Dense_Poly
             ( n,d,c : natural32 ) return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    t : Term;

    procedure Generate_Monomials
                ( accu : in out Term; k,sum : in natural32 ) is

    -- DESCRIPTION :
    --   Accumulating procedure to generate all monomials up to degree d.

    -- ON ENTRY :
    --   accu     accumulator contains the current exponent vector;
    --   k        current component;
    --   sum      sum of current entries in the accumulator.

    -- ON RETURN :
    --   accu     accumulator determined up to k component, k included.

    begin
      if k > n then
        accu.cf := Random_Coefficient(c);
        Add(res,accu);
      else
        for i in 0..d loop
          if sum + i <= d then
            accu.dg(integer32(k)) := i;
            Generate_Monomials(accu,k+1,sum+i);
            accu.dg(integer32(k)) := 0;
          end if;
        end loop;
      end if;
    end Generate_Monomials;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    Generate_Monomials(t,1,0);
    return res;
  end Random_Dense_Poly;

  function Random_Dense_Poly
             ( n : natural32; d,e : integer32; c : natural32 ) 
             return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;
    t : Term;

    procedure Generate_Monomials
                ( accu : in out Term;
                  k : in natural32; sum : in integer32 ) is

    -- DESCRIPTION :
    --   Accumulating procedure to generate all monomials 
    --   with degrees between d and e.

    -- ON ENTRY :
    --   accu     accumulator contains the current exponent vector;
    --   k        current component;
    --   sum      sum of current entries in the accumulator.

    -- ON RETURN :
    --   accu     accumulator determined up to k component, k included.

    begin
      if k > n then
        accu.cf := Random_Coefficient(c);
        Add(res,accu);
      else
        for i in d..e loop
          if sum + i <= e then
            accu.dg(integer32(k)) := i;
            Generate_Monomials(accu,k+1,sum+i);
            accu.dg(integer32(k)) := 0;
          end if;
        end loop;
      end if;
    end Generate_Monomials;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..integer32(n) => 0);
    Generate_Monomials(t,1,0);
    return res;
  end Random_Dense_Poly;

  function Random_Sparse_Poly
             ( n,d,m,c : natural32 )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;

  begin
    for i in 1..m loop
      declare
        rt : Term := Random_Term(n,d,c);
      begin
        Add(res,rt);
        Clear(rt);
      end;
    end loop;
    return res;
  end Random_Sparse_Poly;

  function Random_Sparse_Poly
             ( n : natural32; d,e : integer32; m,c : natural32 )
             return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;

  begin
    for i in 1..m loop
      declare
        rt : Term := Random_Term(n,d,e,c);
      begin
        Add(res,rt);
        Clear(rt);
      end;
    end loop;
    return res;
  end Random_Sparse_Poly;

end Standard_Random_Polynomials;
