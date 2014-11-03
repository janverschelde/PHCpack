with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Random_Matrices;
with DoblDobl_Random_Matrices;
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

end Random_Conditioned_Matrices;
