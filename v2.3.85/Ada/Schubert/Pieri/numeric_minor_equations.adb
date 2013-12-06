with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Brackets;                           use Brackets;
with Plane_Representations;              use Plane_Representations;
with Evaluated_Minors;                   use Evaluated_Minors;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;

package body Numeric_Minor_Equations is

  tol : constant double_float := 10.0**(-10);

-- EXPANDING ACCORDING A BRACKET MONOMIAL :

  function Expanded_Minors
               ( cffmat : Standard_Floating_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bm : Bracket_Monomial ) return Poly is

    res : Poly := Null_Poly;
    first : boolean := true;

    procedure Visit_Bracket ( b : in Bracket; continue : out boolean ) is

      factor : double_float;
      minor : Poly;

    begin
      if first then
        declare
          bb : constant Bracket(b'first..b'last-1) := b(b'first+1..b'last);
        begin
          factor := Determinant(cffmat,bb);
        end;
        first := false;
      else
        minor := Expanded_Minor(polmat,b);
        if (minor /= Null_Poly) and (abs(factor) > tol) then
          Mul(minor,Create(factor));
          Add(res,minor);
        end if;
        Clear(minor);
      end if;
      continue := true;
    end Visit_Bracket;
    procedure Visit_Brackets is new Enumerate_Brackets(Visit_Bracket);

  begin
    Visit_Brackets(bm);
    return res;
  end Expanded_Minors;

  function Expanded_Minors
               ( cffmat : Standard_Complex_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bm : Bracket_Monomial ) return Poly is

    res : Poly := Null_Poly;
    first : boolean := true;
    factor : Complex_Number;

    procedure Visit_Bracket ( b : in Bracket; continue : out boolean ) is

      minor : Poly;

    begin
      if first then
        declare
          bb : constant Bracket(b'first..b'last-1) := b(b'first+1..b'last);
        begin
          factor := Determinant(cffmat,bb);
        end;
        first := false;
      else
        minor := Expanded_Minor(polmat,b);
        if (minor /= Null_Poly) and (AbsVal(factor) > tol) then
          Mul(minor,factor);
          Add(res,minor);
        end if;
        Clear(minor);
      end if;
      continue := true;
    end Visit_Bracket;
    procedure Visit_Brackets is new Enumerate_Brackets(Visit_Bracket);

  begin
    Visit_Brackets(bm);
    return res;
  end Expanded_Minors;

  function Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bm : Bracket_Monomial ) return Poly is

    res : Poly := Null_Poly;
    first : boolean := true;
    factor : Poly;

    procedure Visit_Bracket ( b : in Bracket; continue : out boolean ) is

      minor : Poly;

    begin
      if first then
        declare
          bb : constant Bracket(b'first..b'last-1) := b(b'first+1..b'last);
        begin
          factor := Expanded_Minor(cntmat,bb);
        end;
        first := false;
      else
        minor := Expanded_Minor(polmat,b);
        if (minor /= Null_Poly) and (factor /= Null_Poly) then
          Mul(minor,factor);
          Add(res,minor);
        end if;
        Clear(factor); Clear(minor);
      end if;
      continue := true;
    end Visit_Bracket;
    procedure Visit_Brackets is new Enumerate_Brackets(Visit_Bracket);

  begin
    Visit_Brackets(bm);
    return res;
  end Expanded_Minors;

  function Lifted_Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bm : Bracket_Monomial ) return Poly is

    res : Poly := Null_Poly;
    first : boolean := true;
    factor : Poly;

    procedure Visit_Bracket ( b : in Bracket; continue : out boolean ) is

      minor,extmin : Poly;

    begin
      if first then
        declare
          bb : constant Bracket(b'first..b'last-1) := b(b'first+1..b'last);
        begin
          factor := Expanded_Minor(cntmat,bb);
        end;
        first := false;
      else
        minor := Expanded_Minor(polmat,b);
        if (minor /= Null_Poly) and (factor /= Null_Poly) then
          extmin := Extend_Zero_Lifting(minor); 
          Mul(extmin,factor);
          Add(res,extmin);
        end if;
        Clear(factor); Clear(minor); Clear(extmin);
      end if;
      continue := true;
    end Visit_Bracket;
    procedure Visit_Brackets is new Enumerate_Brackets(Visit_Bracket);

  begin
    Visit_Brackets(bm);
    return res;
  end Lifted_Expanded_Minors;

-- EXPANDING ACCORDING A BRACKET POLYNOMIAL :

  function Expanded_Minors
               ( cffmat : Standard_Floating_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bp : Bracket_Polynomial ) return Poly is

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Bracket_Term; continue : out boolean ) is

      minor : Poly := Expanded_Minors(cffmat,polmat,t.monom);

    begin
      if REAL_PART(t.coeff) > 0.0
       then Add(res,minor);
       else Sub(res,minor);
      end if;
      Clear(minor);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Enumerate_Terms(Visit_Term);

  begin
    Visit_Terms(bp);
    return res;
  end Expanded_Minors;

  function Expanded_Minors
               ( cffmat : Standard_Complex_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bp : Bracket_Polynomial ) return Poly is

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Bracket_Term; continue : out boolean ) is

      minor : Poly := Expanded_Minors(cffmat,polmat,t.monom);

    begin
      if REAL_PART(t.coeff) > 0.0
       then Add(res,minor);
       else Sub(res,minor);
      end if;
      Clear(minor);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Enumerate_Terms(Visit_Term);

  begin
    Visit_Terms(bp);
    return res;
  end Expanded_Minors;

  function Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bp : Bracket_Polynomial ) return Poly is

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Bracket_Term; continue : out boolean ) is

      minor : Poly := Expanded_Minors(cntmat,polmat,t.monom);

    begin
      if REAL_PART(t.coeff) > 0.0
       then Add(res,minor);
       else Sub(res,minor);
      end if;
      Clear(minor);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Enumerate_Terms(Visit_Term);

  begin
    Visit_Terms(bp);
    return res;
  end Expanded_Minors;

  function Lifted_Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bp : Bracket_Polynomial ) return Poly is

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Bracket_Term; continue : out boolean ) is

      minor : Poly := Lifted_Expanded_Minors(cntmat,polmat,t.monom);

    begin
      if REAL_PART(t.coeff) > 0.0
       then Add(res,minor);
       else Sub(res,minor);
      end if;
      Clear(minor);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Enumerate_Terms(Visit_Term);

  begin
    Visit_Terms(bp);
    return res;
  end Lifted_Expanded_Minors;

-- EXPANDING TO CONSTRUCT POLYNOMIAL SYSTEMS :

  function Expanded_Minors
               ( cffmat : Standard_Floating_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bs : Bracket_System ) return Poly_Sys is

    res : Poly_Sys(bs'first+1..bs'last);

  begin
    for i in res'range loop
      res(i) := Expanded_Minors(cffmat,polmat,bs(i));
    end loop;
    return res;
  end Expanded_Minors;

  function Expanded_Minors
               ( cffmat : Standard_Complex_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bs : Bracket_System ) return Poly_Sys is

    res : Poly_Sys(bs'first+1..bs'last);

  begin
    for i in res'range loop
      res(i) := Expanded_Minors(cffmat,polmat,bs(i));
    end loop;
    return res;
  end Expanded_Minors;

  function Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bs : Bracket_System ) return Poly_Sys is

    res : Poly_Sys(bs'first+1..bs'last);

  begin
    for i in res'range loop
      res(i) := Expanded_Minors(cntmat,polmat,bs(i));
    end loop;
    return res;
  end Expanded_Minors;

  function Lifted_Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bs : Bracket_System ) return Poly_Sys is

    res : Poly_Sys(bs'first+1..bs'last);

  begin
    for i in res'range loop
      res(i) := Lifted_Expanded_Minors(cntmat,polmat,bs(i));
    end loop;
    return res;
  end Lifted_Expanded_Minors;

  function Evaluate ( p : Poly; x : Standard_Complex_Matrices.Matrix )
                    return Complex_Number is

    xv : constant Standard_Complex_Vectors.Vector := Vector_Rep(x);

  begin
    return Eval(p,xv);
  end Evaluate;

  function Evaluate ( p : Poly_Sys; x : Standard_Complex_Matrices.Matrix )
                    return Standard_Complex_Vectors.Vector is  

    xv : constant Standard_Complex_Vectors.Vector := Vector_Rep(x);

  begin
    return Eval(p,xv);
  end Evaluate;

  procedure Embed ( t : in out Term ) is

    dg : constant Degrees
       := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last+1);

  begin
    dg(t.dg'range) := t.dg.all;
    dg(dg'last) := 0;
    Clear(t.dg);
    t.dg := dg;
  end Embed;

  function Embed ( t : Term ) return Term is

    res : Term;

  begin
    res.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last+1);
    res.dg(t.dg'range) := t.dg.all;
    res.dg(res.dg'last) := 0;
    res.cf := t.cf;
    return res;
  end Embed;

  procedure Embed ( p : in out Poly ) is

    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      et : Term := Embed(t);

    begin
      Add(res,et);
      Clear(et);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(p);
    p := res;
  end Embed;

  procedure Embed ( p : in out Poly_Sys ) is
  begin
    for i in p'range loop
      Embed(p(i));
    end loop;
  end Embed;

  procedure Embed ( m : in out Standard_Complex_Poly_Matrices.Matrix ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j) /= Null_Poly
         then Embed(m(i,j));
        end if;
      end loop;
    end loop;
  end Embed;

  function Linear_Homotopy ( target,start : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Embed_Target_Term ( t : in Term; continue : out boolean ) is

      et : Term := Embed(t);

    begin
      et.dg(et.dg'last) := 1;
      Add(res,et);
      Clear(et);
      continue := true;
    end Embed_Target_Term;
    procedure Embed_Target_Terms is new Visiting_Iterator(Embed_Target_Term);

    procedure Embed_Start_Term ( t : in Term; continue : out boolean ) is

      et : Term := Embed(t);

    begin
      Add(res,et);
      et.dg(et.dg'last) := 1;
      Sub(res,et);
      Clear(et);
      continue := true;
    end Embed_Start_Term;
    procedure Embed_Start_Terms is new Visiting_Iterator(Embed_Start_Term);

  begin
    Embed_Target_Terms(target);
    Embed_Start_Terms(start);
    return res;
  end Linear_Homotopy;

  function Linear_Interpolation
              ( target,start : Poly; k : integer32 ) return Poly is

    res : Poly := Null_Poly;

    procedure Embed_Target_Term ( t : in Term; continue : out boolean ) is

      et : Term;

    begin
      Copy(t,et);
      et.dg(k) := et.dg(k) + 1;         -- multiply with t
      Add(res,et);
      Clear(et);
      continue := true;
    end Embed_Target_Term;
    procedure Embed_Target_Terms is new Visiting_Iterator(Embed_Target_Term);

    procedure Embed_Start_Term ( t : in Term; continue : out boolean ) is

      et : Term;

    begin
      Copy(t,et);
      Add(res,et);                       -- res := res + et
      et.dg(k) := et.dg(k) + 1;          -- multiply with t
      Sub(res,et);                       -- res := res + et - t*et
      Clear(et);
      continue := true;
    end Embed_Start_Term;
    procedure Embed_Start_Terms is new Visiting_Iterator(Embed_Start_Term);

  begin
    Embed_Target_Terms(target);
    Embed_Start_Terms(start);
    return res;
  end Linear_Interpolation;

  procedure Divide_Common_Factor ( p : in out Poly; k : in integer32 ) is

    first : boolean := true;
    min : natural32;

    procedure Min_Power ( t : in Term; continue : out boolean ) is
    begin
      if first then
        first := false;
        min := t.dg(k);              -- initialize minimal power
      else
        if t.dg(k) < min
         then min := t.dg(k);
        end if;
      end if;
      continue := true;
    end Min_Power;
    procedure Find_Min_Power is new Visiting_Iterator(Min_Power);

    procedure Divide ( t : in out Term; continue : out boolean ) is
    begin
      t.dg(k) := t.dg(k) - min;
      continue := true;
    end Divide;
    procedure Divide_Min_Power is new Changing_Iterator(Divide);

  begin
    Find_Min_Power(p);
    if min > 0
     then Divide_Min_Power(p);
    end if;
  end Divide_Common_Factor;

end Numeric_Minor_Equations;
