with unchecked_deallocation;
with Coefficient_Supported_Polynomials; use Coefficient_Supported_Polynomials;
with QuadDobl_Gradient_Evaluations;     use QuadDobl_Gradient_Evaluations;

package body QuadDobl_Gradient_Circuits is

-- DATA STRUCTURE :

  type Circuit_Rep ( m : integer32 ) is record   -- m is number of terms
    n : natural32;                               -- number of variables
    c : QuadDobl_Complex_Vectors.Vector(1..m);   -- coefficients
    b : Standard_Natural_VecVecs.VecVec(1..m);   -- positions are bits
    f : Standard_Natural_VecVecs.Link_to_VecVec; -- common factors
  end record;

  -- NOTE :
  --   The factors are common to all partial derivatives and separating
  --   those factors from the products of variables absorbds the costs
  --   of higher powers in the polynomials.  The absence of the factors
  --   occurs frequently enough that the factors are stored via a pointer.
  --   The absence of common factors is then simply stored with null,
  --   instead of a vector of zero vectors.

-- CONSTRUCTORS :  

  function Create ( n : natural32;
                    c : QuadDobl_Complex_Vectors.Vector;
                    b : Standard_Natural_VecVecs.VecVec )
                  return Circuit is

    res : Circuit;
    rep : Circuit_Rep(c'last);

  begin
    rep.n := n;
    rep.c := c;
    rep.b := Standard_Natural_VecVecs.Create_Copy(b);
    rep.f := null;
    res := new Circuit_Rep'(rep);
    return res;
  end Create;

  function Create ( n : natural32;
                    c : QuadDobl_Complex_Vectors.Vector;
                    b,f : Standard_Natural_VecVecs.VecVec )
                  return Circuit is

    res : Circuit;
    rep : Circuit_Rep(c'last);

  begin
    rep.n := n;
    rep.c := c;
    rep.b := Standard_Natural_VecVecs.Create_Copy(b);
    declare
      cf : constant Standard_Natural_VecVecs.VecVec(f'range)
         := Standard_Natural_VecVecs.Create_Copy(f);
    begin
      rep.f := new Standard_Natural_VecVecs.VecVec'(cf);
    end;
    res := new Circuit_Rep'(rep);
    return res;
  end Create;

  function Create ( p : Poly ) return Circuit is

    res : Circuit;
    m : constant integer32 
      := integer32(QuadDobl_Complex_Polynomials.Number_of_Terms(p));
    rep : Circuit_Rep(m);
    n : constant natural32
      := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p);
    c : QuadDobl_Complex_Vectors.Vector(1..m);
    e,f,b : Standard_Natural_VecVecs.VecVec(1..m);
    nof : boolean;

  begin
    Coefficients_and_Supports(p,c,e);
    Split_Common_Factors(e,f,b,nof);
    rep.n := n;
    rep.c := c;
    rep.b := b;
    if nof then
      Standard_Natural_VecVecs.Clear(f);
      rep.f := null;
    else
      rep.f := new Standard_Natural_VecVecs.VecVec'(f);
    end if;
    res := new Circuit_Rep'(rep);
    return res;
  end Create;

-- SELECTORS :

  function Number_of_Terms ( c : Circuit ) return natural32 is
  begin
    if c = null
     then return 0;
     else return natural32(c.m);
    end if;
  end Number_of_Terms;

  function Number_of_Variables ( c : Circuit ) return natural32 is
  begin
    if c = null
     then return 0;
     else return natural32(c.n);
    end if;
  end Number_of_Variables;

  function Coefficients
             ( c : Circuit ) return QuadDobl_Complex_Vectors.Vector is
  begin
    return c.c;
  end Coefficients;

  function Coefficient
             ( c : Circuit; k : integer32 ) return Complex_Number is
  begin
    return c.c(k);
  end Coefficient;

  function Positions
             ( c : Circuit ) return Standard_Natural_VecVecs.VecVec is
  begin
    return c.b;
  end Positions;

  function Positions
             ( c : Circuit; k : integer32 )
             return Standard_Natural_Vectors.Link_to_Vector is
  begin
    return c.b(k);
  end Positions;

  function Factors
             ( c : Circuit )
             return Standard_Natural_VecVecs.Link_to_VecVec is
  begin
    return c.f;
  end Factors;

  function Factors
             ( c : Circuit; k : integer32 )
             return Standard_Natural_Vectors.Link_to_Vector is
  begin
    return c.f(k);
  end Factors;

-- EVALUATION AND DIFFERENTIATION :

  function WorkSpace ( c : Circuit )
                     return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..c.m);
    dim : integer32 := integer32(c.n);

  begin
    for k in res'range loop
      res(k) := new QuadDobl_Complex_Vectors.Vector(0..dim);
    end loop;
    return res;
  end WorkSpace;

  procedure EvalDiff ( c : in Circuit;
                       x : in QuadDobl_Complex_Vectors.Vector;
                       wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                       ydx : out QuadDobl_Complex_Vectors.Vector ) is
    
    use Standard_Natural_VecVecs;

  begin
    if c.f = null
     then Gradient_of_Polynomial(c.b,c.c,x,wrk,ydx);
     else Gradient_of_Polynomial(c.f.all,c.b,c.c,x,wrk,ydx);
    end if;
  end EvalDiff;

  function EvalDiff ( c : Circuit; x : QuadDobl_Complex_Vectors.Vector )
                    return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(0..x'last);
    wrk : QuadDobl_Complex_VecVecs.VecVec(1..c.m) := WorkSpace(c);

  begin
    EvalDiff(c,x,wrk,res);
    QuadDobl_Complex_VecVecs.Clear(wrk);
    return res;
  end EvalDiff;

-- DESTRUCTORS :

  procedure Clear ( c : in out Circuit_Rep ) is
  begin
    Standard_Natural_VecVecs.Clear(c.b);
    Standard_Natural_VecVecs.Deep_Clear(c.f);
  end Clear;

  procedure Clear ( c : in out Circuit ) is

    procedure free is new unchecked_deallocation(Circuit_Rep,Circuit);

  begin
    if c /= null then
      Clear(c.all);
      free(c);
    end if;
  end Clear;

end QuadDobl_Gradient_Circuits;
