with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Multprec_Natural_Numbers;          use Multprec_Natural_Numbers;
with Multprec_Integer_Numbers;          use Multprec_Integer_Numbers;
with Standard_Floating_Linear_Solvers;  use Standard_Floating_Linear_Solvers;
with Double_Double_Linear_Solvers;      use Double_Double_Linear_Solvers;
with Quad_Double_Linear_Solvers;        use Quad_Double_Linear_Solvers;
with Multprec_Floating_Linear_Solvers;  use Multprec_Floating_Linear_Solvers;
with VarbPrec_Matrix_Conversions;

package body Varbprec_Floating_Linear_Solvers is

  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Standard_Floating_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out double_float; loss : out integer32 ) is

    dim : constant integer32 := mtx'last(1);

  begin
    lufco(mtx,dim,piv,rco);
    if rco > 0.0
     then loss := integer32(log10(rco));
     else loss := -2**30;
    end if;
  end Estimated_Loss_of_Decimal_Places;

  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Double_Double_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out double_double; loss : out integer32 ) is

    dim : constant integer32 := mtx'last(1);

  begin
    lufco(mtx,dim,piv,rco);
    if rco > 0.0
     then loss := integer32(log10(to_double(rco)));
     else loss := -2**30;
    end if;
  end Estimated_Loss_of_Decimal_Places;

  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Quad_Double_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out quad_double; loss : out integer32 ) is

    dim : constant integer32 := mtx'last(1);

  begin
    lufco(mtx,dim,piv,rco);
    if rco > 0.0
     then loss := integer32(log10(to_double(rco)));
     else loss := -2**30;
    end if;
  end Estimated_Loss_of_Decimal_Places;

  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Multprec_Floating_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out Floating_Number; loss : out integer32 ) is

    dim : constant integer32 := mtx'last(1);

  begin
    lufco(mtx,dim,piv,rco);
    if Equal(rco,0.0) then
      loss := -2**30;
    else
      declare
        dp : constant natural32 := Decimal_Places(Unsigned(Fraction(rco)));
        exprco : Integer_Number := Exponent(rco);
      begin
        loss := Create(exprco) + integer32(dp-1);
      end;
    end if;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : Standard_Floating_Matrices.Matrix )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last(1);
    wrk : Standard_Floating_Matrices.Matrix(mtx'range(1),mtx'range(2)) := mtx;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : double_float;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : Double_Double_Matrices.Matrix )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last(1);
    wrk : Double_Double_Matrices.Matrix(mtx'range(1),mtx'range(2)) := mtx;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : double_double;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : Quad_Double_Matrices.Matrix )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last(1);
    wrk : Quad_Double_Matrices.Matrix(mtx'range(1),mtx'range(2)) := mtx;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : quad_double;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : Multprec_Floating_Matrices.Matrix )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last(1);
    wrk : Multprec_Floating_Matrices.Matrix(mtx'range(1),mtx'range(2));
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : Floating_Number;

  begin
    Multprec_Floating_Matrices.Copy(mtx,wrk);
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    Clear(rco);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Standard_Floating_Matrices.Matrix;
                rhs : in out Standard_Floating_Vectors.Vector;
                want : in integer32; fail : out boolean;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out double_float; loss : out integer32 ) is

    precision : constant integer32 := 16;

  begin
    Estimated_Loss_of_Decimal_Places(mtx,piv,rco,loss);
    fail := (precision + loss < want);
    if not fail
     then lusolve(mtx,mtx'last(1),piv,rhs);
    end if;
  end Solve_to_Wanted_Decimal_Places;

  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Double_Double_Matrices.Matrix;
                rhs : in out Double_Double_Vectors.Vector;
                want : in integer32; fail : out boolean;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out double_double; loss : out integer32 ) is

    precision : constant integer32 := 32;

  begin
    Estimated_Loss_of_Decimal_Places(mtx,piv,rco,loss);
    fail := (precision + loss < want);
    if not fail
     then lusolve(mtx,mtx'last(1),piv,rhs);
    end if;
  end Solve_to_Wanted_Decimal_Places;

  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Quad_Double_Matrices.Matrix;
                rhs : in out Quad_Double_Vectors.Vector;
                want : in integer32; fail : out boolean;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out quad_double; loss : out integer32 ) is

    precision : constant integer32 := 64;

  begin
    Estimated_Loss_of_Decimal_Places(mtx,piv,rco,loss);
    fail := (precision + loss < want);
    if not fail
     then lusolve(mtx,mtx'last(1),piv,rhs);
    end if;
  end Solve_to_Wanted_Decimal_Places;

  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Multprec_Floating_Matrices.Matrix;
                rhs : in out Multprec_Floating_Vectors.Vector;
                want : in integer32;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out Floating_Number; loss : out integer32 ) is

    num1 : Floating_Number := mtx(1,1);
    mat_precision : integer32 := integer32(Decimal_Places_Fraction(num1));
    precision : integer32 := mat_precision;
    size : natural32;

  begin
    Estimated_Loss_of_Decimal_Places(mtx,piv,rco,loss);
    precision := precision + loss;
    if precision < want then
      precision := mat_precision + (want - precision);
      size := Multprec_Floating_Numbers.Decimal_to_Size(natural32(precision));
      VarbPrec_Matrix_Conversions.Set_Size(mtx,size);
    end if;
    lusolve(mtx,mtx'last(1),piv,rhs);
  end Solve_to_Wanted_Decimal_Places;

end Varbprec_Floating_Linear_Solvers;
