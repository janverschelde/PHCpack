with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Floating_Linear_Solvers;  use Standard_Floating_Linear_Solvers;
with Double_Double_Linear_Solvers;      use Double_Double_Linear_Solvers;
with Quad_Double_Linear_Solvers;        use Quad_Double_Linear_Solvers;

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

end Varbprec_Floating_Linear_Solvers;
