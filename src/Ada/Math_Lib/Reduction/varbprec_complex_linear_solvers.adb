with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Multprec_Natural_Numbers;          use Multprec_Natural_Numbers;
with Multprec_Integer_Numbers;          use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Standard_Complex_Linear_Solvers;   use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;   use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;   use QuadDobl_Complex_Linear_Solvers;
with Multprec_Complex_Linear_Solvers;   use Multprec_Complex_Linear_Solvers;
with Varbprec_VecVec_Conversions;
with Varbprec_Matrix_Conversions;

package body Varbprec_Complex_Linear_Solvers is

  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Standard_Complex_Matrices.Matrix;
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
              ( mtx : in out DoblDobl_Complex_Matrices.Matrix;
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
              ( mtx : in out QuadDobl_Complex_Matrices.Matrix;
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
              ( mtx : in out Multprec_Complex_Matrices.Matrix;
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

  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Standard_Complex_VecVecs.VecVec;
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
              ( mtx : in out DoblDobl_Complex_VecVecs.VecVec;
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
              ( mtx : in out QuadDobl_Complex_VecVecs.VecVec;
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
              ( mtx : in out Multprec_Complex_VecVecs.VecVec;
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
              ( mtx : Standard_Complex_Matrices.Matrix )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last(1);
    wrk : Standard_Complex_Matrices.Matrix(mtx'range(1),mtx'range(2)) := mtx;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : double_float;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : DoblDobl_Complex_Matrices.Matrix )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last(1);
    wrk : DoblDobl_Complex_Matrices.Matrix(mtx'range(1),mtx'range(2)) := mtx;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : double_double;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : QuadDobl_Complex_Matrices.Matrix )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last(1);
    wrk : QuadDobl_Complex_Matrices.Matrix(mtx'range(1),mtx'range(2)) := mtx;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : quad_double;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : Multprec_Complex_Matrices.Matrix )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last(1);
    wrk : Multprec_Complex_Matrices.Matrix(mtx'range(1),mtx'range(2));
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : Floating_Number;

  begin
    Multprec_Complex_Matrices.Copy(mtx,wrk);
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    Clear(rco);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : Standard_Complex_VecVecs.VecVec )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last;
    wrk : Standard_Complex_VecVecs.VecVec(mtx'range)
        := Standard_Complex_VecVecs.Create_Copy(mtx);
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : double_float;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    Standard_Complex_VecVecs.Clear(wrk);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : DoblDobl_Complex_VecVecs.VecVec )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last;
    wrk : DoblDobl_Complex_VecVecs.VecVec(mtx'range)
        := DoblDobl_Complex_VecVecs.Create_Copy(mtx);
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : double_double;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    DoblDobl_Complex_VecVecs.Clear(wrk);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : QuadDobl_Complex_VecVecs.VecVec )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last;
    wrk : QuadDobl_Complex_VecVecs.VecVec(mtx'range)
        := QuadDobl_Complex_VecVecs.Create_Copy(mtx);
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : quad_double;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    QuadDobl_Complex_VecVecs.Clear(wrk);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  function Estimated_Loss_of_Decimal_Places
              ( mtx : Multprec_Complex_VecVecs.VecVec )
              return integer32 is

    res : integer32;
    dim : constant integer32 := mtx'last;
    wrk : Multprec_Complex_VecVecs.VecVec(mtx'range)
        := Multprec_Complex_VecVecs.Create_Copy(mtx);
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : Floating_Number;

  begin
    Estimated_Loss_of_Decimal_Places(wrk,piv,rco,res);
    Multprec_Complex_VecVecs.Clear(wrk);
    Clear(rco);
    return res;
  end Estimated_Loss_of_Decimal_Places;

  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Standard_Complex_Matrices.Matrix;
                rhs : in out Standard_Complex_Vectors.Vector;
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
              ( mtx : in out DoblDobl_Complex_Matrices.Matrix;
                rhs : in out DoblDobl_Complex_Vectors.Vector;
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
              ( mtx : in out QuadDobl_Complex_Matrices.Matrix;
                rhs : in out QuadDobl_Complex_Vectors.Vector;
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
              ( mtx : in out Multprec_Complex_Matrices.Matrix;
                rhs : in out Multprec_Complex_Vectors.Vector;
                want : in integer32;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out Floating_Number; loss : out integer32 ) is

    num1 : Multprec_Complex_Numbers.Complex_Number := mtx(1,1);
    rpt1 : Floating_Number := Multprec_Complex_Numbers.REAL_PART(num1);
    mat_precision : integer32 := integer32(Decimal_Places_Fraction(rpt1));
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
    Multprec_Floating_Numbers.Clear(rpt1);
  end Solve_to_Wanted_Decimal_Places;

  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Standard_Complex_VecVecs.VecVec;
                rhs : in out Standard_Complex_Vectors.Vector;
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
              ( mtx : in out DoblDobl_Complex_VecVecs.VecVec;
                rhs : in out DoblDobl_Complex_Vectors.Vector;
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
              ( mtx : in out QuadDobl_Complex_VecVecs.VecVec;
                rhs : in out QuadDobl_Complex_Vectors.Vector;
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
              ( mtx : in out Multprec_Complex_VecVecs.VecVec;
                rhs : in out Multprec_Complex_Vectors.Vector;
                want : in integer32;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out Floating_Number; loss : out integer32 ) is

    num1 : Multprec_Complex_Numbers.Complex_Number := mtx(1)(1);
    rpt1 : Floating_Number := Multprec_Complex_Numbers.REAL_PART(num1);
    mat_precision : integer32 := integer32(Decimal_Places_Fraction(rpt1));
    precision : integer32 := mat_precision;
    size : natural32;

  begin
    Estimated_Loss_of_Decimal_Places(mtx,piv,rco,loss);
    precision := precision + loss;
    if precision < want then
      precision := mat_precision + (want - precision);
      size := Multprec_Floating_Numbers.Decimal_to_Size(natural32(precision));
      Varbprec_VecVec_Conversions.Set_Size(mtx,size);
    end if;
    lusolve(mtx,mtx'last(1),piv,rhs);
    Multprec_Floating_Numbers.Clear(rpt1);
  end Solve_to_Wanted_Decimal_Places;

end Varbprec_Complex_Linear_Solvers;
