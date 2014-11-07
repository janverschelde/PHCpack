with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Multprec_Floating_Numbers;
with Standard_Random_Matrices;
with DoblDobl_Random_Matrices;
with QuadDobl_Random_Matrices;
with Multprec_Random_Matrices;
with VarbPrec_Matrix_Conversions;       use VarbPrec_Matrix_Conversions;

package body Random_Conditioned_Matrices is

  function Singular_Value_Matrix
             ( n : integer32; c : double_float )
             return Standard_Floating_Matrices.Matrix is

   res : Standard_Floating_Matrices.Matrix(1..n,1..n);
   deg : constant integer32 := integer32(log10(abs(c)));

  begin
    for i in 1..n loop
      for j in 1..n loop
        res(i,j) := 0.0;
      end loop;
      if i = 1 then
        res(i,i) := 1.0;
      elsif i = n then
        res(i,i) := c;
      else
        res(i,i) := 10.0**(integer((i-1)*deg/(n-1)));
      end if;
    end loop;
    return res;
  end Singular_Value_Matrix;

  function Singular_Value_Matrix
             ( n : integer32; c : double_float )
             return Standard_Complex_Matrices.Matrix is

   res : Standard_Complex_Matrices.Matrix(1..n,1..n);
   deg : constant integer32 := integer32(log10(abs(c)));
   rdi : double_float;

  begin
    for i in 1..n loop
      for j in 1..n loop
        res(i,j) := Standard_Complex_Numbers.Create(0.0);
      end loop;
      if i = 1 then
        res(i,i) := Standard_Complex_Numbers.Create(1.0);
      elsif i = n then
        res(i,i) := Standard_Complex_Numbers.Create(c);
      else
        rdi := 10.0**(integer((i-1)*deg/(n-1)));
        res(i,i) := Standard_Complex_Numbers.Create(rdi);
      end if;
    end loop;
    return res;
  end Singular_Value_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Standard_Floating_Matrices.Matrix is

    use Standard_Random_Matrices;

    res : Standard_Floating_Matrices.Matrix(1..n,1..n);
    svm : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
        := Singular_Value_Matrix(n,c);
    rq1 : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 
    rq2 : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 

    use Standard_Floating_Matrices;

  begin
    res := rq1*svm*rq2;
    return res;
  end Random_Conditioned_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Standard_Complex_Matrices.Matrix is

    use Standard_Random_Matrices;

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);
    svm : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := Singular_Value_Matrix(n,c);
    rq1 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 
    rq2 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 

    use Standard_Complex_Matrices;

  begin
    res := rq1*svm*rq2;
    return res;
  end Random_Conditioned_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Double_Double_Matrices.Matrix is

    use DoblDobl_Random_Matrices;

    res : Double_Double_Matrices.Matrix(1..n,1..n);
    st_svm : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
           := Singular_Value_Matrix(n,c);
    dd_svm : constant Double_Double_Matrices.Matrix(1..n,1..n)
           := d2dd(st_svm);
    rq1 : constant Double_Double_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 
    rq2 : constant Double_Double_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 

    use Double_Double_Matrices;

  begin
    res := rq1*dd_svm*rq2;
    return res;
  end Random_Conditioned_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return DoblDobl_Complex_Matrices.Matrix is

    use DoblDobl_Random_Matrices;

    res : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    st_svm : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
           := Singular_Value_Matrix(n,c);
    dd_svm : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
           := d2dd(st_svm);
    rq1 : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 
    rq2 : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 

    use DoblDobl_Complex_Matrices;

  begin
    res := rq1*dd_svm*rq2;
    return res;
  end Random_Conditioned_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Quad_Double_Matrices.Matrix is

    use QuadDobl_Random_Matrices;

    res : Quad_Double_Matrices.Matrix(1..n,1..n);
    st_svm : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
           := Singular_Value_Matrix(n,c);
    dd_svm : constant Quad_Double_Matrices.Matrix(1..n,1..n)
           := d2qd(st_svm);
    rq1 : constant Quad_Double_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 
    rq2 : constant Quad_Double_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 

    use Quad_Double_Matrices;

  begin
    res := rq1*dd_svm*rq2;
    return res;
  end Random_Conditioned_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return QuadDobl_Complex_Matrices.Matrix is

    use QuadDobl_Random_Matrices;

    res : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    st_svm : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
           := Singular_Value_Matrix(n,c);
    dd_svm : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
           := d2qd(st_svm);
    rq1 : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 
    rq2 : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 

    use QuadDobl_Complex_Matrices;

  begin
    res := rq1*dd_svm*rq2;
    return res;
  end Random_Conditioned_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Multprec_Floating_Matrices.Matrix is

    res : Multprec_Floating_Matrices.Matrix(1..n,1..n);
    magc : constant natural32 := natural32(log10(abs(c)));
    deci : natural32 := 2*magc;
    size : natural32;
    st_svm : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
           := Singular_Value_Matrix(n,c);
    mp_svm : Multprec_Floating_Matrices.Matrix(1..n,1..n)
           := d2mp(st_svm);
    rq1,rq2 : Multprec_Floating_Matrices.Matrix(1..n,1..n);

    use Multprec_Floating_Matrices,Multprec_Random_Matrices;

  begin
    if deci < 16
     then deci := 16; -- at least double precision !
    end if;
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Set_Size(mp_svm,size);
    Set_Size(rq1,size);
    rq1 := Random_Orthogonal_Matrix(natural32(n),natural32(n),size);
    rq2 := Random_Orthogonal_Matrix(natural32(n),natural32(n),size);
    res := rq1*mp_svm*rq2;
    Multprec_Floating_Matrices.Clear(mp_svm);
    Multprec_Floating_Matrices.Clear(rq1);
    Multprec_Floating_Matrices.Clear(rq2);
    return res;
  end Random_Conditioned_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Multprec_Complex_Matrices.Matrix is

    res : Multprec_Complex_Matrices.Matrix(1..n,1..n);
    magc : constant natural32 := natural32(log10(abs(c)));
    deci : natural32 := 2*magc;
    size : natural32;
    st_svm : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
           := Singular_Value_Matrix(n,c);
    mp_svm : Multprec_Complex_Matrices.Matrix(1..n,1..n)
           := d2mp(st_svm);
    rq1,rq2 : Multprec_Complex_Matrices.Matrix(1..n,1..n);

    use Multprec_Complex_Matrices,Multprec_Random_Matrices;

  begin
    if deci < 16
     then deci := 16; -- at least double precision !
    end if;
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Set_Size(mp_svm,size);
    Set_Size(rq1,size);
    rq1 := Random_Orthogonal_Matrix(natural32(n),natural32(n),size);
    rq2 := Random_Orthogonal_Matrix(natural32(n),natural32(n),size);
    res := rq1*mp_svm*rq2;
    Multprec_Complex_Matrices.Clear(mp_svm);
    Multprec_Complex_Matrices.Clear(rq1);
    Multprec_Complex_Matrices.Clear(rq2);
    return res;
  end Random_Conditioned_Matrix;

end Random_Conditioned_Matrices;
