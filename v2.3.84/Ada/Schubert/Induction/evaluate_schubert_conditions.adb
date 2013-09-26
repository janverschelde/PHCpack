with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;

package body Evaluate_Schubert_Conditions is

-- DESCRIPTION :
--   This package provides functions to evaluate the Schubert
--   intersection conditions.

  function Eval ( X,P,B : Standard_Complex_Matrices.Matrix;
                  h,m : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector is

    n : constant natural := X'last(1);
    np1 : constant natural := n+1;
    res : Standard_Complex_Vectors.Vector(1..np1);
    Btm : Standard_Complex_Vectors.Vector(1..n);
    ind : integer;

  begin
    for i in 1..n loop                   -- compute B times m
      Btm(i) := Create(0.0);
      for j in B'range(2) loop
        Btm(i) := Btm(i) + B(i,j)*m(j);
      end loop;
    end loop;
    for i in 1..n loop                   -- compute [X|P] times B*m
      res(i) := Create(0.0); 
      ind := Btm'first-1;
      for j in X'range(2) loop
        ind := ind + 1;
        res(i) := res(i) + X(i,j)*Btm(ind);
      end loop;
      for j in P'range(2) loop
        ind := ind + 1;
        res(i) := res(i) + P(i,j)*Btm(ind);
      end loop;
    end loop;
    res(np1) := Create(-1.0);
    for i in h'range loop                -- compute -1+<h,m>
      res(np1) := res(np1) + h(i)*m(i);
    end loop;
    return res;
  end Eval;

end Evaluate_Schubert_Conditions;
