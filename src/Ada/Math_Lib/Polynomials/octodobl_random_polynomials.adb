with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with Standard_Random_Numbers;
with OctoDobl_Random_Numbers;           use OctoDobl_Random_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;

package body OctoDobl_Random_Polynomials is

  function Random_Coefficient ( c : natural32 ) return Complex_Number is
 
    res : Complex_Number;
    ranflt : octo_double;

  begin
    case c is
      when 1 => res := Create(integer(1));
      when 2 => ranflt := Random;
                res := Create(ranflt);
      when others => res := Random1;
    end case;
    return res;
  end Random_Coefficient;

  function Random_Monomial ( n,d : natural32 ) return Term is

    res : Term;
    deg,pos : integer32;

  begin
    res.cf := Create(integer(1));
    res.dg := new Vector'(1..integer32(n) => 0);
    for i in 1..d loop
      deg := Standard_Random_Numbers.Random(0,1);
      pos := Standard_Random_Numbers.Random(1,integer32(n));
      res.dg(pos) := res.dg(pos) + natural32(deg);
    end loop;
    return res;
  end Random_Monomial;

  function Random_Term ( n,d,c : natural32 ) return Term is

    res : Term := Random_Monomial(n,d);

  begin
    res.cf := Random_Coefficient(c);
    return res;
  end Random_Term;

  function Random_Dense_Poly ( n,d,c : natural32 ) return Poly is

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
    t.dg := new Vector'(1..integer32(n) => 0);
    Generate_Monomials(t,1,0);
    return res;
  end Random_Dense_Poly;

  function Random_Sparse_Poly ( n,d,m,c : natural32 ) return Poly is

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

end OctoDobl_Random_Polynomials;
