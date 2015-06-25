with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with Coefficient_Supported_Polynomials;  use Coefficient_Supported_Polynomials;
with Standard_Gradient_Evaluations;
with Standard_Jacobian_Evaluations;
with DoblDobl_Gradient_Evaluations;

package body DoblDobl_Jacobian_Evaluations is

  procedure EvalDiff
              ( b : in Standard_Natural_VecVecs.VecVec;
                c : in DoblDobl_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in DoblDobl_Complex_Vectors.Vector;
                z : out DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_VecVecs.VecVec;
                A : out DoblDobl_Complex_Matrices.Matrix ) is

    ind : integer32;
    use DoblDobl_Complex_Numbers;
    cff : Complex_Number;
    zero : constant double_double := create(0.0);

  begin
    DoblDobl_Gradient_Evaluations.Gradient_Monomials(b,x,y);
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
                c : in DoblDobl_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in DoblDobl_Complex_Vectors.Vector;
                z : out DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_VecVecs.VecVec;
                A : out DoblDobl_Complex_Matrices.Matrix ) is

    ind : integer32;
    use DoblDobl_Complex_Numbers;
    cff : Complex_Number;
    zero : constant double_double := create(0.0);

  begin
    DoblDobl_Gradient_Evaluations.Gradient_Monomials(f,b,x,y);
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
                c : in DoblDobl_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in DoblDobl_Complex_Vectors.Vector;
                z : out DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_VecVecs.VecVec;
                A : in DoblDobl_Complex_VecVecs.VecVec ) is

    ind : integer32;
    use DoblDobl_Complex_Numbers;
    cff : Complex_Number;
    yind,Arow : DoblDobl_Complex_Vectors.Link_to_Vector;
    zero : constant double_double := create(0.0);

  begin
    DoblDobl_Gradient_Evaluations.Gradient_Monomials(b,x,y);
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
                c : in DoblDobl_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in DoblDobl_Complex_Vectors.Vector;
                z : out DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_VecVecs.VecVec;
                A : in DoblDobl_Complex_VecVecs.VecVec ) is

    ind : integer32;
    use DoblDobl_Complex_Numbers;
    cff : Complex_Number;
    yind,Arow : DoblDobl_Complex_Vectors.Link_to_Vector;
    zero : constant double_double := create(0.0);

  begin
    DoblDobl_Gradient_Evaluations.Gradient_Monomials(f,b,x,y);
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

  procedure DoblDobl_Jacobian_Evaluation
              ( v : in Standard_Natural_VecVecs.VecVec;
                c : in DoblDobl_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in DoblDobl_Complex_Vectors.Vector;
                z : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix ) is

    f,b : Standard_Natural_VecVecs.VecVec(v'range);
    y : DoblDobl_Complex_VecVecs.VecVec(v'range);
    ind : integer32;
    use DoblDobl_Complex_Numbers;
    cff : Complex_Number;

  begin
    Split_Common_Factors(v,f,b);
    y := DoblDobl_Gradient_Evaluations.Gradient_Monomials(f,b,x);
    for i in z'range loop
      z(i) := Create(integer(0));
      for j in A'range(2) loop
        A(i,j) := Create(integer(0));
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
  end DoblDobl_Jacobian_Evaluation;

  procedure DoblDobl_Jacobian_Evaluation
              ( v : in Standard_Integer_VecVecs.VecVec;
                c : in DoblDobl_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in DoblDobl_Complex_Vectors.Vector;
                z : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix ) is

    e : Standard_Natural_VecVecs.VecVec(v'range)
      := Standard_Jacobian_Evaluations.Integer_to_Natural(v);

  begin
    DoblDobl_Jacobian_Evaluation(e,c,k,x,z,A);
    Standard_Natural_VecVecs.Clear(e);
  end DoblDobl_Jacobian_Evaluation;

end DoblDobl_Jacobian_Evaluations;
