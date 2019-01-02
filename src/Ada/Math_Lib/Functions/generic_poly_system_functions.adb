with unchecked_deallocation;

package body Generic_Poly_System_Functions is 

  use Polynomials;

-- CREATORS :

  function Create ( p : Poly_Sys ) return Eval_Poly_Sys is
    
    res : Eval_Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Create(p(k));
    end loop;
    return res;
  end Create;

  function Create ( p : Poly_Sys ) return Eval_Coeff_Poly_Sys is
    
    res : Eval_Coeff_Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Create(p(k));
    end loop;
    return res;
  end Create;

  function Coeff ( p : Poly_Sys ) return VecVec is

    res : VecVec(p'range);

  begin
    for i in p'range loop
      res(i) := new Vector'(Coeff(p(i)));
    end loop;
    return res;
  end Coeff;

-- EVALUATORS :

  function Eval ( p : Poly_Sys; x : number; i : integer32 ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for j in p'range loop
      res(j) := Eval(p(j),x,i);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Poly_Sys; x : Vector ) return Vector is

    res : Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Eval(p(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Eval_Poly_Sys; x : Vector ) return Vector is

    res : Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Eval(p(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Eval_Coeff_Poly_Sys; c : VecVec; x : Vector )
                return Vector is

    res : Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Eval(p(i),c(i).all,x);
    end loop;
    return res;
  end Eval;

-- DESTRUCTORS :

  procedure Clear ( p : in out Eval_Poly_Sys ) is
  begin
    for k in p'range loop
      Clear(p(k));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Eval_Poly_Sys ) is

    procedure free is 
      new unchecked_deallocation(Eval_Poly_Sys,Link_to_Eval_Poly_Sys);

  begin
    if p /= null
     then Clear(p.all);
          free(p);
    end if;
  end Clear;

  procedure Clear ( p : in out Eval_Coeff_Poly_Sys ) is
  begin
    for k in p'range loop
      Clear(p(k));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Eval_Coeff_Poly_Sys ) is

    procedure free is
      new unchecked_deallocation(Eval_Coeff_Poly_Sys,
                                 Link_to_Eval_Coeff_Poly_Sys);

  begin
    if p /= null
     then Clear(p.all);
          free(p);
    end if;
  end Clear;

end Generic_Poly_System_Functions;
