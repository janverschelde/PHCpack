with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Gradient_Evaluations;

package body Standard_Jacobian_Evaluations is

  function Integer_to_Natural
              ( v : Standard_Integer_VecVecs.VecVec )
              return Standard_Natural_VecVecs.VecVec is

    res : Standard_Natural_VecVecs.VecVec(v'range);

  begin
    for i in v'range loop
      res(i) := new Standard_Natural_Vectors.Vector(v(i)'range);
      for j in v(i)'range loop
        res(i)(j) := natural32(v(i)(j));
      end loop;
    end loop;
    return res;
  end Integer_to_Natural;

  procedure Standard_Jacobian_Evaluation
              ( f,b : in Standard_Natural_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix ) is

    y : Standard_Complex_VecVecs.VecVec(f'range);
    ind : integer32;
    use Standard_Complex_Numbers;
    cff : Complex_Number;

  begin
    y := Standard_Gradient_Evaluations.Gradient_Monomials(f,b,x);
    for i in z'range loop
      z(i) := Create(0.0);
      for j in A'range(2) loop
        A(i,j) := Create(0.0);
      end loop;
      for j in c(i)'range loop
        ind := integer32(k(i)(j));
        cff := c(i)(j);
        z(i) := z(i) + cff*y(ind)(0);
        for k in A'range(2) loop
          A(i,k) := A(i,k) + cff*y(ind)(k);
        end loop;
      end loop;
    end loop;
  end Standard_Jacobian_Evaluation;

  procedure Standard_Jacobian_Evaluation
              ( v : in Standard_Natural_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix ) is

    f,b : Standard_Natural_VecVecs.VecVec(v'range);

  begin
    Standard_Gradient_Evaluations.Split_Common_Factors(v,f,b);
    Standard_Jacobian_Evaluation(f,b,c,k,x,z,A);
  end Standard_Jacobian_Evaluation;

  procedure Standard_Jacobian_Evaluation
              ( v : in Standard_Integer_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix ) is

    e : Standard_Natural_VecVecs.VecVec(v'range) := Integer_to_Natural(v);

  begin
    Standard_Jacobian_Evaluation(e,c,k,x,z,A);
    Standard_Natural_VecVecs.Clear(e);
  end Standard_Jacobian_Evaluation;

end Standard_Jacobian_Evaluations;
