with unchecked_deallocation;
with Coefficient_Supported_Polynomials; use Coefficient_Supported_Polynomials;
with Standard_Gradient_Evaluations;     use Standard_Gradient_Evaluations;

package body Standard_Gradient_Circuits is

-- DATA STRUCTURE :

  type Circuit_Rep ( m : integer32 ) is record   -- m is number of terms
    n : natural32;                               -- number of variables
    c : Standard_Complex_Vectors.Vector(1..m);   -- coefficients
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
                    c : Standard_Complex_Vectors.Vector;
                    b : Standard_Natural_VecVecs.VecVec )
                  return Circuit is

    res : Circuit;
    rep : Circuit_Rep(c'last);

  begin
    rep.n := n;
    rep.c := c;
    rep.b := b;
    rep.f := null;
    res := new Circuit_Rep'(rep);
    return res;
  end Create;

  function Create ( n : natural32;
                    c : Standard_Complex_Vectors.Vector;
                    b,f : Standard_Natural_VecVecs.VecVec )
                  return Circuit is

    res : Circuit;
    rep : Circuit_Rep(c'last);

  begin
    rep.n := n;
    rep.c := c;
    rep.b := b;
    rep.f := new Standard_Natural_VecVecs.VecVec'(f);
    res := new Circuit_Rep'(rep);
    return res;
  end Create;

  function Create ( p : Poly ) return Circuit is

    res : Circuit;
    m : constant integer32 
      := integer32(Standard_Complex_Polynomials.Number_of_Terms(p));
    n : constant natural32
      := Standard_Complex_Polynomials.Number_of_Unknowns(p);
    c : Standard_Complex_Vectors.Vector(1..m);
    e,f,b : Standard_Natural_VecVecs.VecVec(1..m);
    nof : boolean;

  begin
    Coefficients_and_Supports(p,c,e);
    Split_Common_Factors(e,f,b,nof);
    if nof then
      Standard_Natural_VecVecs.Clear(f);
      res := Create(n,c,b);
    else
      res := Create(n,c,b,f);
    end if;
    return res;
  end Create;

-- SELECTORS :

  function Number_of_Terms ( crc : Circuit ) return natural32 is
  begin
    if crc = null
     then return 0;
     else return natural32(crc.m);
    end if;
  end Number_of_Terms;

  function Number_of_Variables ( crc : Circuit ) return natural32 is
  begin
    if crc = null
     then return 0;
     else return natural32(crc.n);
    end if;
  end Number_of_Variables;

  function Coefficients
             ( crc : Circuit ) return Standard_Complex_Vectors.Vector is
  begin
    return crc.c;
  end Coefficients;

  function Coefficient
             ( crc : Circuit; k : integer32 ) return Complex_Number is
  begin
    return crc.c(k);
  end Coefficient;

  function Positions
             ( crc : Circuit ) return Standard_Natural_VecVecs.VecVec is
  begin
    return crc.b;
  end Positions;

  function Positions
             ( crc : Circuit; k : integer32 )
             return Standard_Natural_Vectors.Link_to_Vector is
  begin
    return crc.b(k);
  end Positions;

  function Factors
             ( crc : Circuit )
             return Standard_Natural_VecVecs.Link_to_VecVec is
  begin
    return crc.f;
  end Factors;

  function Factors
             ( crc : Circuit; k : integer32 )
             return Standard_Natural_Vectors.Link_to_Vector is
  begin
    return crc.f(k);
  end Factors;

-- EVALUATION AND DIFFERENTIATION :

  function WorkSpace ( crc : Circuit )
                     return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(1..crc.m);
    dim : integer32 := integer32(crc.n);

  begin
    for k in res'range loop
      res(k) := new Standard_Complex_Vectors.Vector(0..dim);
    end loop;
    return res;
  end WorkSpace;

  procedure EvalDiff ( crc : in Circuit;
                       x : in Standard_Complex_Vectors.Vector;
                       wrk : in out Standard_Complex_VecVecs.VecVec;
                       ydx : out Standard_Complex_Vectors.Vector ) is
    
    use Standard_Natural_VecVecs;

  begin
    if crc.f = null
     then Gradient_of_Polynomial(crc.b,crc.c,x,wrk,ydx);
     else Gradient_of_Polynomial(crc.f.all,crc.b,crc.c,x,wrk,ydx);
    end if;
  end EvalDiff;

  function EvalDiff ( c : Circuit; x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(0..x'last);
    wrk : Standard_Complex_VecVecs.VecVec(1..c.m) := WorkSpace(c);

  begin
    EvalDiff(c,x,wrk,res);
    Standard_Complex_VecVecs.Clear(wrk);
    return res;
  end EvalDiff;

-- DESTRUCTORS :

  procedure Clear ( crc : in out Circuit_Rep ) is
  begin
    Standard_Natural_VecVecs.Clear(crc.b);
    Standard_Natural_VecVecs.Deep_Clear(crc.f);
  end Clear;

  procedure Clear ( crc : in out Circuit ) is

    procedure free is new unchecked_deallocation(Circuit_Rep,Circuit);

  begin
    if crc /= null then
      Clear(crc.all);
      free(crc);
    end if;
  end Clear;

end Standard_Gradient_Circuits;
