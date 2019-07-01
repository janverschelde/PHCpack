with unchecked_deallocation;
with QuadDobl_Rational_Approximations;

package body QuadDobl_Pade_Approximants is

-- CONSTRUCTORS :

  function Create ( num,den : QuadDobl_Complex_Vectors.Vector )
                  return Pade_Rep is

    res : Pade_Rep(num'last,den'last);

  begin
    res.num := num;
    res.den := den;
    return res;
  end Create;

  function Create ( num,den : QuadDobl_Complex_Vectors.Vector )
                  return Pade is

    rep : constant Pade_Rep := Create(num,den);
    res : constant Pade := new Pade_Rep'(rep);

  begin
    return res;
  end Create;

  function Coefficients ( srv : QuadDobl_Complex_Series_Vectors.Vector;
                          idx : integer32 )
                        return QuadDobl_Complex_Vectors.Vector is

    dim : constant integer32 := srv(idx).deg;
    res : QuadDobl_Complex_Vectors.Vector(0..dim);

  begin
    for i in res'range loop
      res(i) := srv(idx).cff(i);
    end loop;
    return res;
  end Coefficients;

  function Create ( numdeg,dendeg : integer32;
                    srv : QuadDobl_Complex_Series_Vectors.Vector;
                    verbose : boolean := false )
                  return Pade_Vector is

    res : Pade_Vector(srv'range);

  begin
    for i in srv'range loop
      declare
        cff : constant QuadDobl_Complex_Vectors.Vector
            := Coefficients(srv,i);
        num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
        den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
        info : integer32;
      begin
        QuadDobl_Rational_Approximations.Pade
          (numdeg,dendeg,cff,num,den,info,verbose);
        res(i) := Create(num,den);
      end;
    end loop;
    return res;
  end Create;

  function Allocate ( numdeg,dendeg : integer32 ) return Pade is

    zero : constant Complex_Number := Create(integer32(0));
    num : constant QuadDobl_Complex_Vectors.Vector(0..numdeg)
        := (0..numdeg => zero);
    den : constant QuadDobl_Complex_Vectors.Vector(0..dendeg)
        := (0..dendeg => zero);
    res : constant Pade := Create(num,den);

  begin
    return res;
  end Allocate;

  function Allocate ( dim,numdeg,dendeg : integer32 ) return Pade_Vector is

    res : Pade_Vector(1..dim);

  begin
    for i in 1..dim loop
      res(i) := Allocate(numdeg,dendeg);
    end loop;
    return res;
  end Allocate;

  procedure Create ( pv : in out Pade_Vector;
                     srv : in QuadDobl_Complex_Series_Vectors.Vector;
                     verbose : in boolean := false ) is

    numdeg : constant integer32 := pv(pv'first).numdeg;
    dendeg : constant integer32 := pv(pv'first).dendeg;
    info : integer32;

  begin
    for i in srv'range loop
      declare
        cff : constant QuadDobl_Complex_Vectors.Vector
            := Coefficients(srv,i);
      begin
        QuadDobl_Rational_Approximations.Pade
          (numdeg,dendeg,cff,pv(i).num,pv(i).den,info,verbose);
      end;
    end loop;
  end Create;

-- SELECTORS :

  function Numerator_Degree ( p : Pade ) return integer32 is
  begin
    if p = null
     then return -1;
     else return p.numdeg;
    end if;
  end Numerator_Degree;

  function Denominator_Degree ( p : Pade ) return integer32 is
  begin
    if p = null
     then return -1;
     else return p.dendeg;
    end if;
  end Denominator_Degree;

  function Numerator_Coefficients
             ( p : Pade ) return QuadDobl_Complex_Vectors.Vector is
  begin
    return p.num;
  end Numerator_Coefficients;

  function Denominator_Coefficients
             ( p : Pade ) return QuadDobl_Complex_Vectors.Vector is
  begin
    return p.den;
  end Denominator_Coefficients;

-- EVALUATORS :

  function Eval ( p : Pade; x : quad_double ) return Complex_Number is

    cx : constant Complex_Number := Create(x);

  begin
    return Eval(p,cx);
  end Eval;

  function Eval ( p : Pade_Vector; x : quad_double )
                return QuadDobl_Complex_Vectors.Vector is

    cx : constant Complex_Number := Create(x);

  begin
    return Eval(p,cx);
  end Eval;

  function Eval ( p : Pade; x : Complex_Number ) return Complex_Number is

    num : constant QuadDobl_Complex_Vectors.Vector
        := Numerator_Coefficients(p);
    den : constant QuadDobl_Complex_Vectors.Vector
        := Denominator_Coefficients(p);
    res : constant Complex_Number
        := QuadDobl_Rational_Approximations.Evaluate(num,den,x);

  begin
    return res;
  end Eval;

  function Eval ( p : Pade_Vector; x : Complex_Number )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      res(k) := Eval(p(k),x);
    end loop;
    return res;
  end Eval;

-- DESTRUCTORS :

  procedure Clear ( p : in out Pade ) is

    procedure free is new unchecked_deallocation(Pade_Rep,Pade);

  begin
    free(p);
  end Clear;

  procedure Clear ( p : in out Pade_Vector ) is
  begin
    for k in p'range loop
      Clear(p(k));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Pade_Vector ) is

    procedure free is
      new unchecked_deallocation(Pade_Vector,Link_to_Pade_Vector);

  begin
    if p /= null then
      for k in p'range loop
        Clear(p(k));
      end loop;
      free(p);
    end if;
  end Clear;

  procedure Clear ( p : in out Pade_VecVec ) is
  begin
    for k in p'range loop
      Clear(p(k));
    end loop;
  end Clear;

end QuadDobl_Pade_Approximants;
