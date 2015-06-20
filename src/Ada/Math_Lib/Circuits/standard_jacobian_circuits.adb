with unchecked_deallocation;
with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with Standard_Complex_VecVecs;
with Lists_of_Integer_Vectors;
with Lexicographical_Supports;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Polynomial_Flatteners;
with Coefficient_Supported_Polynomials;
with Standard_Jacobian_Evaluations;      use Standard_Jacobian_Evaluations;

package body Standard_Jacobian_Circuits is

-- DATA STRUCTURE :

  type Circuit_Rep ( m,nq,nv : integer32 ) is record
    -- m  : distinct number of monomials 
    -- nq : the number of equations 
    -- nv : the number of variables
    b : Standard_Natural_VecVecs.VecVec(1..m);   -- positions are bits
    f : Standard_Natural_VecVecs.Link_to_VecVec; -- common factors
    c : Standard_Complex_VecVecs.VecVec(1..nq);  -- coefficients
    k : Standard_Natural_VecVecs.VecVec(1..nq);  -- indices
  end record;

-- CONSTRUCTORS :

  function Create ( p : Poly_Sys ) return Circuit is

    res : Circuit;
    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    c : Standard_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    s,e : Lists_of_Integer_Vectors.List;
    nof : boolean;

  begin
    s := Standard_Polynomial_Flatteners.Distinct_Supports(p);
    e := Lexicographical_Supports.Sort(s);
    declare
      m : constant integer32
        := integer32(Lists_of_Integer_Vectors.Length_Of(e));
      v : Standard_Integer_VecVecs.VecVec(1..m)
        := Lists_of_Integer_Vectors.Shallow_Create(e);
      w : Standard_Natural_VecVecs.VecVec(1..m)
        := Standard_Jacobian_Evaluations.Integer_to_Natural(v);
      b,f : Standard_Natural_VecVecs.VecVec(1..m);
      rep : Circuit_Rep(m,nq,nv);
    begin
      Coefficient_Supported_Polynomials.Split_Common_Factors(w,f,b,nof);
      rep.b := b;
      if nof then
        Standard_Natural_VecVecs.Clear(f);
        rep.f := null;
      else
        rep.f := new Standard_Natural_VecVecs.VecVec'(f);
      end if;
      Standard_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
      rep.k := k;
      rep.c := c;
      res := new Circuit_Rep'(rep);
    end;
    return res;
  end Create;

-- SELECTORS :

  function Number_of_Polynomials ( c : Circuit ) return natural32 is
  begin
    if c = null
     then return 0;
     else return natural32(c.nq);
    end if;
  end Number_of_Polynomials;

  function Number_of_Variables ( c : Circuit ) return natural32 is
  begin
    if c = null
     then return 0;
     else return natural32(c.nv);
    end if;
  end Number_of_Variables;

  function Number_of_Monomials ( c : Circuit ) return natural32 is
  begin
    if c = null
     then return 0;
     else return natural32(c.m);
    end if;
  end Number_of_Monomials;

  function Number_of_Terms ( c : Circuit; i : integer32 ) return natural32 is
  begin
    if c = null then
      return 0;
    elsif ((i < 1) or (i > c.nq)) then
      return 0;
    else
      return natural32(c.k(i)'last);
    end if;
  end Number_of_Terms;

  function Coefficients ( c : Circuit; i : integer32 )
                        return Standard_Complex_Vectors.Link_to_Vector is
  begin
    if c = null then
      return null;
    elsif ((i < 1) or (i > c.nq)) then
      return null;
    else
      return c.c(i);
    end if;
  end Coefficients;

  function Coefficient ( c : Circuit; i,j : integer32 )
                       return Complex_Number is

    res : Complex_Number;
    cff : constant Standard_Complex_Vectors.Link_to_Vector
        := Coefficients(c,i);

    use Standard_Complex_Vectors;

  begin
    if cff = null
     then res := Create(0.0);
     else res := cff(j); 
    end if;
    return res;
  end Coefficient;

  function Product ( c : Circuit; i,j : integer32 )
                   return Standard_Natural_Vectors.Link_to_Vector is

    ind : integer32;

  begin
    if c = null then
      return null;
    elsif ((i < 1) or (i > c.nq)) then
      return null;
    else
      ind := integer32(c.k(i)(j));
      return c.b(ind);
    end if;
  end Product;

  function Factor ( c : Circuit; i,j : integer32 )
                  return Standard_Natural_Vectors.Link_to_Vector is

    ind : integer32;

    use Standard_Natural_VecVecs;

  begin
    if c = null then
      return null;
    elsif ((i < 1) or (i > c.nq)) then
      return null;
    elsif c.f = null then
      return null;
    else
      ind := integer32(c.k(i)(j));
      return c.f(ind);
    end if;
  end Factor;

-- DESTRUCTOR :

  procedure Clear ( c : in out Circuit ) is

    procedure free is new unchecked_deallocation(Circuit_Rep,Circuit);

  begin
    if c /= null then
      free(c);
    end if;
  end Clear;

end Standard_Jacobian_Circuits;
