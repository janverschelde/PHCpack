with unchecked_deallocation;

package body Generic_Polynomial_Vectors is

-- EVALUATION and DIFFERENTIATION :

  function Eval ( p : Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector is

    res : Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomials.Eval(p(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Link_to_Polynomial_Vector;
                  x : Vectors.Vector ) return Vectors.Vector is

    res : Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomials.Eval(p(i),x);
    end loop;
    return res;
  end Eval;

  procedure Diff ( p : in Polynomial_Vector; x : in Vectors.Vector;
                   m : out Matrices.Matrix ) is

    yd,wrk : Vectors.Vector(x'range);

  begin
    for i in p'range loop
      Polynomials.Diff(p(i).mons,x,yd,wrk);
      for j in yd'range loop
        Ring.Copy(yd(j),m(i,j));
      end loop;
    end loop;
  end Diff;

  procedure Diff ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                   m : out Matrices.Matrix ) is
  begin
    if p /= null
     then Diff(p.all,x,m);
    end if;
  end Diff;

  procedure Speel ( p : in Polynomial_Vector; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    for i in p'range loop
      Polynomials.Speel(p(i).mons,x,y(i),yd,wrk);
      for j in yd'range loop
        Ring.Copy(yd(j),m(i,j));
      end loop;
    end loop;
  end Speel;

  procedure Speel ( p : in Polynomial_Vector; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix ) is

    yd,wrk : Vectors.Vector(x'range);

  begin
    Speel(p,x,y,m,yd,wrk);
  end Speel;

  procedure Speel ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    if p /= null
     then Speel(p.all,x,y,m,yd,wrk);
    end if;
  end Speel;

  procedure Speel ( p : in Link_to_Polynomial_Vector; x : in Vectors.Vector;
                    y : out Vectors.Vector; m : out Matrices.Matrix ) is
  begin
    if p /= null
     then Speel(p.all,x,y,m);
    end if;
  end Speel;

-- DESTRUCTORS :

  procedure Clear ( p : in out Polynomial_Vector ) is
  begin
    for i in p'range loop
      Polynomials.Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Polynomial_Vector ) is

    procedure free is
      new unchecked_deallocation(Polynomial_Vector,Link_to_Polynomial_Vector);

  begin
    if p /= null then
      Clear(p.all);
      free(p);
    end if;
  end Clear;

end Generic_Polynomial_Vectors;
