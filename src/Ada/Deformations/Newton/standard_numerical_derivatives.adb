with Standard_Floating_Numbers;         use Standard_Floating_Numbers;

package body Standard_Numerical_Derivatives is

  function Diff1 ( f : Complex_Multivariate_Function; x : Vector;
                   i : integer32; h : double_float ) return Complex_Number is

    res : Complex_Number;
    x_pls,x_min : Vector(x'range) := x;

  begin
    x_pls(i) := x(i) + Create(h);
    x_min(i) := x(i) - Create(h);
    res := (f(x_pls) - f(x_min))/Create(2.0*h);
    return res;
  end Diff1;

  function Diff2 ( f : Complex_Multivariate_Function; x : Vector;
                   i : integer32; h : double_float ) return Complex_Number is

    res : Complex_Number;

  begin
    res := (Diff1(f,x,i,h/2.0) - Diff1(f,x,i,h)/Create(4.0))/Create(0.75);
    return res;
  end Diff2;

  function Diff3 ( f : Complex_Multivariate_Function; x : Vector;
                   i : integer32; h : double_float ) return Complex_Number is

    res : Complex_Number;

  begin
    res := (Diff2(f,x,i,h/2.0) - Diff2(f,x,i,h)/Create(16.0))
           /Create(15.0/16.0);
    return res;
  end Diff3;

  function Diff1 ( f : Complex_Multivariate_System; n : integer32;
                   x : Vector; h : double_float ) return Matrix is

    res : Matrix(1..n,x'range);
    x_pls,x_min : Vector(x'range);
    y_pls,y_min : Vector(1..n);

  begin
    for j in x'range loop
      x_pls := x; x_pls(j) := x(j) + Create(h);
      x_min := x; x_min(j) := x(j) - Create(h);
      y_pls := f(x_pls);
      y_min := f(x_min);
      for i in 1..n loop
        res(i,j) := (y_pls(i) - y_min(i))/Create(2.0*h);
      end loop;
    end loop;
    return res;
  end Diff1;

  function Diff2 ( f : Complex_Multivariate_System; n : integer32;
                   x : Vector; h : double_float ) return Matrix is

    res : Matrix(1..n,x'range);
    dh1 : Matrix(1..n,x'range) := Diff1(f,n,x,h);
    dh2 : Matrix(1..n,x'range) := Diff1(f,n,x,h/2.0);

  begin
    for i in 1..n loop
      for j in x'range loop
        res(i,j) := (dh2(i,j) - dh1(i,j)/Create(4.0))/Create(0.75);
      end loop;
    end loop;
    return res;
  end Diff2;

  function Diff3 ( f : Complex_Multivariate_System; n : integer32;
                   x : Vector; h : double_float ) return Matrix is

    res : Matrix(1..n,x'range);
    dh1 : Matrix(1..n,x'range) := Diff2(f,n,x,h);
    dh2 : Matrix(1..n,x'range) := Diff2(f,n,x,h/2.0);

  begin
    for i in 1..n loop
      for j in x'range loop
        res(i,j) := (dh2(i,j) - dh1(i,j)/Create(16.0))/Create(15.0/16.0);
      end loop;
    end loop;
    return res;
  end Diff3;

end Standard_Numerical_Derivatives;
