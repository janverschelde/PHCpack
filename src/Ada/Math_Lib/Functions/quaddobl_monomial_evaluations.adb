with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Monomial_Evaluations;

package body QuadDobl_Monomial_Evaluations is

  function Eval ( e : Standard_Natural_Vectors.Vector;
                  x : QuadDobl_Complex_Vectors.Vector )
                return Complex_Number is

    res : Complex_Number := Create(integer(1));

  begin
    for i in e'range loop
       if e(i) > 0 then
         for j in 1..e(i) loop
           res := res*x(i);
         end loop;
       end if;
    end loop;
    return res;
  end Eval;

  function Eval ( e : Standard_Natural_VecVecs.VecVec;
                  x : QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(e'range);

  begin
    for i in e'range loop
      res(i) := Eval(e(i).all,x);
    end loop;
    return res;
  end Eval;

  function Power_Table
                ( n,m : integer32;
                  d : Standard_Natural_Vectors.Vector;
                  x : QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(1..n,1..m);

  begin
    for i in 1..n loop
      res(i,1) := x(i);
      for j in 2..integer32(d(i)) loop
        res(i,j) := res(i,j-1)*x(i);
      end loop;
    end loop;
    return res;
  end Power_Table;

  function Eval ( e : Standard_Natural_Vectors.Vector;
                  p : QuadDobl_Complex_Matrices.Matrix )
                return Complex_Number is

    res : Complex_Number := Create(integer(1));

  begin
    for i in e'range loop
      if e(i) > 0
       then res := res*p(i,integer32(e(i)));
      end if;
    end loop;
    return res;
  end Eval;
 
  function Eval ( d : Standard_Natural_Vectors.Vector;
                  e : Standard_Natural_VecVecs.VecVec;
                  x : QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(e'range);
    maxdeg : constant integer32
           := integer32(Standard_Monomial_Evaluations.Max(d));
    n : constant integer32 := d'last;

  begin
    if maxdeg = 0 then
      declare
        powtab : QuadDobl_Complex_Matrices.Matrix(1..n,1..1);
      begin
        for i in 1..n loop
          powtab(i,1) := x(i);
        end loop;
        for i in e'range loop
          res(i) := Eval(e(i).all,powtab);
        end loop;
      end;
    else
      declare
        powtab : QuadDobl_Complex_Matrices.Matrix(1..n,1..maxdeg);
      begin
        powtab := Power_Table(n,maxdeg,d,x);
        for i in e'range loop
          res(i) := Eval(e(i).all,powtab);
        end loop;
      end;
    end if;
    return res;
  end Eval;
 
  function Eval_with_Power_Table 
                ( e : Standard_Natural_VecVecs.VecVec;
                  x : QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(e'range);
    d : constant Standard_Natural_Vectors.Vector(x'range)
      := Standard_Monomial_Evaluations.Largest_Degrees(x'last,e);

  begin
    res := Eval(d,e,x);
    return res;
  end Eval_with_Power_Table;

end QuadDobl_Monomial_Evaluations;
