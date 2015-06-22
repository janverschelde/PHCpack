with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Multprec_Gradient_Evaluations;

package body Multprec_Jacobian_Evaluations is

  procedure EvalDiff
              ( b : in Standard_Natural_VecVecs.VecVec;
                c : in Multprec_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Multprec_Complex_Vectors.Vector;
                z : out Multprec_Complex_Vectors.Vector;
                y : out Multprec_Complex_VecVecs.VecVec;
                A : out Multprec_Complex_Matrices.Matrix ) is

    ind : integer32;
    use Multprec_Complex_Numbers;
    cff : Complex_Number;
    zero : Floating_Number := create(0.0);

  begin
    Multprec_Gradient_Evaluations.Gradient_Monomials(b,x,y);
    for i in z'range loop
      z(i) := Create(zero);
      for j in A'range(2) loop
        A(i,j) := Create(zero);
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
  end EvalDiff;

  procedure EvalDiff
              ( f,b : in Standard_Natural_VecVecs.VecVec;
                c : in Multprec_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Multprec_Complex_Vectors.Vector;
                z : out Multprec_Complex_Vectors.Vector;
                y : out Multprec_Complex_VecVecs.VecVec;
                A : out Multprec_Complex_Matrices.Matrix ) is

    ind : integer32;
    use Multprec_Complex_Numbers;
    cff : Complex_Number;
    zero : Floating_Number := create(0.0);

  begin
    Multprec_Gradient_Evaluations.Gradient_Monomials(f,b,x,y);
    for i in z'range loop
      z(i) := Create(zero);
      for j in A'range(2) loop
        A(i,j) := Create(zero);
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
  end EvalDiff;

  procedure EvalDiff
              ( b : in Standard_Natural_VecVecs.VecVec;
                c : in Multprec_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Multprec_Complex_Vectors.Vector;
                z : out Multprec_Complex_Vectors.Vector;
                y : out Multprec_Complex_VecVecs.VecVec;
                A : in Multprec_Complex_VecVecs.VecVec ) is

    ind : integer32;
    use Multprec_Complex_Numbers;
    cff : Complex_Number;
    yind,Arow : Multprec_Complex_Vectors.Link_to_Vector;
    zero : Floating_Number := create(0.0);

  begin
    Multprec_Gradient_Evaluations.Gradient_Monomials(b,x,y);
    for i in z'range loop
      z(i) := Create(zero);
      Arow := A(i);
      for j in Arow'range loop
        Arow(j) := Create(zero);
      end loop;
      for j in c(i)'range loop
        ind := integer32(k(i)(j));
        cff := c(i)(j);
        yind := y(ind);
        Arow := A(i);
        z(i) := z(i) + cff*yind(0);
        for k in Arow'range loop
          Arow(k) := Arow(k) + cff*yind(k);
        end loop;
      end loop;
    end loop;
  end EvalDiff;

  procedure EvalDiff
              ( f,b : in Standard_Natural_VecVecs.VecVec;
                c : in Multprec_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Multprec_Complex_Vectors.Vector;
                z : out Multprec_Complex_Vectors.Vector;
                y : out Multprec_Complex_VecVecs.VecVec;
                A : in Multprec_Complex_VecVecs.VecVec ) is

    ind : integer32;
    use Multprec_Complex_Numbers;
    cff : Complex_Number;
    yind,Arow : Multprec_Complex_Vectors.Link_to_Vector;
    zero : Floating_Number := create(0.0);

  begin
    Multprec_Gradient_Evaluations.Gradient_Monomials(f,b,x,y);
    for i in z'range loop
      z(i) := Create(zero);
      Arow := A(i);
      for j in Arow'range loop
        Arow(j) := Create(zero);
      end loop;
      for j in c(i)'range loop
        ind := integer32(k(i)(j));
        cff := c(i)(j);
        yind := y(ind);
        Arow := A(i);
        z(i) := z(i) + cff*yind(0);
        for k in Arow'range loop
          Arow(k) := Arow(k) + cff*yind(k);
        end loop;
      end loop;
    end loop;
  end EvalDiff;

end Multprec_Jacobian_Evaluations;
