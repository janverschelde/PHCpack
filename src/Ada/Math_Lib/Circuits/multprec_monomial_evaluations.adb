with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Monomial_Evaluations;

package body Multprec_Monomial_Evaluations is

  function Eval ( e : Standard_Natural_Vectors.Vector;
                  x : Multprec_Complex_Vectors.Vector )
                return Complex_Number is

    res : Complex_Number := Create(integer(1));

  begin
    for i in e'range loop
       if e(i) > 0 then
         for j in 1..e(i) loop
           Mul(res,x(i));
         end loop;
       end if;
    end loop;
    return res;
  end Eval;

  function Eval ( e : Standard_Natural_VecVecs.VecVec;
                  x : Multprec_Complex_Vectors.Vector )
                return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(e'range);

  begin
    for i in e'range loop
      res(i) := Eval(e(i).all,x);
    end loop;
    return res;
  end Eval;

  function Power_Table
                ( n,m : integer32;
                  d : Standard_Natural_Vectors.Vector;
                  x : Multprec_Complex_Vectors.Vector )
                return Multprec_Complex_Matrices.Matrix is

    res : Multprec_Complex_Matrices.Matrix(1..n,1..m);

  begin
    for i in 1..n loop
      Copy(x(i),res(i,1));
      for j in 2..integer32(d(i)) loop
        res(i,j) := res(i,j-1)*x(i);
      end loop;
    end loop;
    return res;
  end Power_Table;

  function Eval ( e : Standard_Natural_Vectors.Vector;
                  p : Multprec_Complex_Matrices.Matrix )
                return Complex_Number is

    res : Complex_Number := Create(integer(1));

  begin
    for i in e'range loop
      if e(i) > 0
       then Mul(res,p(i,integer32(e(i))));
      end if;
    end loop;
    return res;
  end Eval;
 
  function Eval ( d : Standard_Natural_Vectors.Vector;
                  e : Standard_Natural_VecVecs.VecVec;
                  x : Multprec_Complex_Vectors.Vector )
                return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(e'range);
    maxdeg : constant integer32
           := integer32(Standard_Monomial_Evaluations.Max(d));
    n : constant integer32 := d'last;

  begin
    if maxdeg = 0 then
      declare
        powtab : Multprec_Complex_Matrices.Matrix(1..n,1..1);
      begin
        for i in 1..n loop
          Copy(x(i),powtab(i,1));
        end loop;
        for i in e'range loop
          res(i) := Eval(e(i).all,powtab);
        end loop;
        Multprec_Complex_Matrices.Clear(powtab);
      end;
    else
      declare
        powtab : Multprec_Complex_Matrices.Matrix(1..n,1..maxdeg);
      begin
        powtab := Power_Table(n,maxdeg,d,x);
        for i in e'range loop
          res(i) := Eval(e(i).all,powtab);
        end loop;
        Multprec_Complex_Matrices.Clear(powtab);
      end;
    end if;
    return res;
  end Eval;
 
  function Eval_with_Power_Table 
                ( e : Standard_Natural_VecVecs.VecVec;
                  x : Multprec_Complex_Vectors.Vector )
                return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(e'range);
    d : constant Standard_Natural_Vectors.Vector(x'range)
      := Standard_Monomial_Evaluations.Largest_Degrees(x'last,e);

  begin
    res := Eval(d,e,x);
    return res;
  end Eval_with_Power_Table;

end Multprec_Monomial_Evaluations;
