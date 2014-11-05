with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Complex_Linear_Solvers;   use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;   use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;   use QuadDobl_Complex_Linear_Solvers;

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

end Varbprec_Complex_Linear_Solvers;
