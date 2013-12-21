with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices_io;    use Standard_Integer64_Matrices_io;
with Standard_Integer_Linear_Solvers;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Standard_Integer_Matrix_Inverse;
with Standard_Smith_Normal_Form;
with Standard_Integer_Kernel;
with Evaluated_Minors; -- to test whether determinant is plus or minus one

package body Standard_Exponent_Transformations is

  procedure UV_Coordinate_Transformation
              ( U,V : in Standard_Integer_Matrices.Matrix;
                M : out Standard_Integer_Matrices.Matrix ) is

    Uinv : Standard_Integer_Matrices.Matrix(U'range(1),U'range(2));
    Vinv : Standard_Integer_Matrices.Matrix(V'range(1),V'range(2));
    E : Standard_Integer_Matrices.Matrix(V'range(1),V'range(2));
    d : integer32;
    use Standard_Integer_Matrices;

  begin
    Standard_Integer_Matrix_Inverse.Inverse(U,d,Uinv);
    for i in Uinv'range(1) loop
      for j in Uinv'range(2) loop
        E(i,j) := Uinv(i,j);
      end loop;
      for j in Uinv'last(2)+1..E'last(2) loop
        E(i,j) := 0;
      end loop;
    end loop;
    for i in Uinv'last(1)+1..E'last(1) loop
      for j in Uinv'range(2) loop
        E(i,j) := 0;
      end loop;
      for j in Uinv'last(2)+1..E'last(2) loop
        if i = j
         then E(i,j) := 1;
         else E(i,j) := 0;
        end if;
      end loop;
    end loop;
    Standard_Integer_Matrix_Inverse.Inverse(V,d,Vinv);
    M := E*Vinv;
  end UV_Coordinate_Transformation;

  procedure UV_Coordinate_Transformation
              ( U,V : in Standard_Integer64_Matrices.Matrix;
                M : out Standard_Integer64_Matrices.Matrix ) is

    Uinv : Standard_Integer64_Matrices.Matrix(U'range(1),U'range(2));
    Vinv : Standard_Integer64_Matrices.Matrix(V'range(1),V'range(2));
    E : Standard_Integer64_Matrices.Matrix(V'range(1),V'range(2));
    use Standard_Integer64_Matrices;
    d : integer64;

  begin
    Standard_Integer_Matrix_Inverse.Inverse(U,d,Uinv);
    for i in Uinv'range(1) loop
      for j in Uinv'range(2) loop
        E(i,j) := Uinv(i,j);
      end loop;
      for j in Uinv'last(2)+1..E'last(2) loop
        E(i,j) := 0;
      end loop;
    end loop;
    for i in Uinv'last(1)+1..E'last(1) loop
      for j in Uinv'range(2) loop
        E(i,j) := 0;
      end loop;
      for j in Uinv'last(2)+1..E'last(2) loop
        if i = j
         then E(i,j) := 1;
         else E(i,j) := 0;
        end if;
      end loop;
    end loop;
    Standard_Integer_Matrix_Inverse.Inverse(V,d,Vinv);
    M := E*Vinv;
  end UV_Coordinate_Transformation;

  function Diagonal_Product 
             ( A : Standard_Integer_Matrices.Matrix ) return integer32 is

    res : integer32 := 1;

  begin
    for i in A'range(1) loop
      exit when (i > A'last(2));
      res := res*A(i,i);
    end loop;
    return res;
  end Diagonal_Product;

  function Diagonal_Product 
             ( A : Standard_Integer64_Matrices.Matrix ) return integer64 is

    res : integer64 := 1;

  begin
    for i in A'range(1) loop
      exit when (i > A'last(2));
      res := res*A(i,i);
    end loop;
    return res;
  end Diagonal_Product;

  procedure Unimodular_Coordinate_Transformation
              ( C : in Standard_Integer_Matrices.Matrix;
                U,T,V,M : out Standard_Integer_Matrices.Matrix;
                p : out integer32; uv,fail : out boolean ) is

    d : integer32;

  begin
    T := Standard_Integer_Matrices.Transpose(C);
    U := Standard_Smith_Normal_Form.Identity(natural32(T'last(1)));
    V := Standard_Smith_Normal_Form.Identity(natural32(T'last(2)));
    Standard_Smith_Normal_Form.Diagonalize(U,T,V);
    p := Diagonal_Product(T);
    if Standard_Integer_Matrix_Inverse.Is_Identity(U) then
      Standard_Integer_Matrix_Inverse.Inverse(V,d,M);
      fail := false; uv := false;
    elsif p = 1 then
      UV_Coordinate_Transformation(U,V,M);
      fail := false; uv := true;
    else
      fail := true; uv := false;
    end if;
  end Unimodular_Coordinate_Transformation;

  procedure Unimodular_Coordinate_Transformation
               ( C : in Standard_Integer64_Matrices.Matrix;
                 U,T,V,M : out Standard_Integer64_Matrices.Matrix;
                 p : out integer64; uv,fail : out boolean ) is

    d : integer64;

  begin
    T := Standard_Integer64_Matrices.Transpose(C);
    U := Standard_Smith_Normal_Form.Identity(natural32(T'last(1)));
    V := Standard_Smith_Normal_Form.Identity(natural32(T'last(2)));
    Standard_Smith_Normal_Form.Diagonalize(U,T,V);
    p := Diagonal_Product(T);
    if Standard_Integer_Matrix_Inverse.Is_Identity(U) then
      Standard_Integer_Matrix_Inverse.Inverse(V,d,M);
      fail := false; uv := false;
    elsif p = 1 then
      UV_Coordinate_Transformation(U,V,M);
      fail := false; uv := true;
    else
      fail := true; uv := false;
    end if;
  end Unimodular_Coordinate_Transformation;

  procedure Unimodular_Coordinate_Transformation
              ( C : in Standard_Integer_Matrices.Matrix;
                fail : out boolean;
                M : out Standard_Integer_Matrices.Matrix ) is

    T : Standard_Integer_Matrices.Matrix(C'range(2),C'range(1))
      := Standard_Integer_Matrices.Transpose(C); 
    U : Standard_Integer_Matrices.Matrix(T'range(1),T'range(1));
     -- := Standard_Smith_Normal_Form.Identity(T'last(1));
    V : Standard_Integer_Matrices.Matrix(T'range(2),T'range(2));
     -- := Standard_Smith_Normal_Form.Identity(T'last(2));
    p : integer32;
    uv : boolean;

  begin
    Unimodular_Coordinate_Transformation(C,U,T,V,M,p,uv,fail);
  end Unimodular_Coordinate_Transformation;

  procedure Unimodular_Coordinate_Transformation
               ( file : in file_type;
                 C : in Standard_Integer_Matrices.Matrix;
                 fail : out boolean;
                 M : out Standard_Integer_Matrices.Matrix ) is

    T : Standard_Integer_Matrices.Matrix(C'range(2),C'range(1))
      := Standard_Integer_Matrices.Transpose(C); 
    U : Standard_Integer_Matrices.Matrix(T'range(1),T'range(1));
     -- := Standard_Smith_Normal_Form.Identity(T'last(1));
    V : Standard_Integer_Matrices.Matrix(T'range(2),T'range(2));
     -- := Standard_Smith_Normal_Form.Identity(T'last(2));
    p : integer32;
    uv : boolean;

  begin
    Unimodular_Coordinate_Transformation(C,U,T,V,M,p,uv,fail);
    put_line(file,"The Smith normal form of the tropisms :"); put(file,T);
    put(file,"The product of the diagonal of the Smith form : ");
    put(file,p,1); put_line(file,".");
    if fail then
      put_line(file,"No integer parameterization possible.");
    else
      if not uv then
        put_line(file,"U is the identity matrix => M = V^-1.");
      else
        put_line(file,"U is not the identity matrix.");
        put(file,"Smith form diagonal product = 1 => M = E V^-1,");
        put_line(file," E contains U^-1.");
      end if;
      put_line(file,"The unimodular transformation M :"); put(M);
    end if;
  end Unimodular_Coordinate_Transformation;

  procedure Unimodular_Coordinate_Transformation
               ( C : in Standard_Integer64_Matrices.Matrix;
                 fail : out boolean;
                 M : out Standard_Integer64_Matrices.Matrix ) is

    T : Standard_Integer64_Matrices.Matrix(C'range(2),C'range(1))
      := Standard_Integer64_Matrices.Transpose(C); 
    U : Standard_Integer64_Matrices.Matrix(T'range(1),T'range(1));
     -- := Standard_Smith_Normal_Form.Identity(T'last(1));
    V : Standard_Integer64_Matrices.Matrix(T'range(2),T'range(2));
     -- := Standard_Smith_Normal_Form.Identity(T'last(2));
    p : integer64;
    uv : boolean;

  begin
    Unimodular_Coordinate_Transformation(C,U,T,V,M,p,uv,fail);
  end Unimodular_Coordinate_Transformation;

  procedure Unimodular_Coordinate_Transformation
               ( file : in file_type;
                 C : in Standard_Integer64_Matrices.Matrix;
                 fail : out boolean;
                 M : out Standard_Integer64_Matrices.Matrix ) is

    T : Standard_Integer64_Matrices.Matrix(C'range(2),C'range(1))
      := Standard_Integer64_Matrices.Transpose(C); 
    U : Standard_Integer64_Matrices.Matrix(T'range(1),T'range(1));
     -- := Standard_Smith_Normal_Form.Identity(T'last(1));
    V : Standard_Integer64_Matrices.Matrix(T'range(2),T'range(2));
     -- := Standard_Smith_Normal_Form.Identity(T'last(2));
    p : integer64;
    uv : boolean;

  begin
    Unimodular_Coordinate_Transformation(C,U,T,V,M,p,uv,fail);
    put_line(file,"The Smith normal form of the tropisms :"); put(file,T);
    put(file,"The product of the diagonal of the Smith form : ");
    Standard_Integer_Numbers_io.put(file,p,1); put_line(file,".");
    if fail then
      put_line(file,"No integer parameterization possible.");
    else
      if not uv then
        put_line(file,"U is the identity matrix => M = V^-1.");
      else
        put_line(file,"U is not the identity matrix.");
        put(file,"Smith form diagonal product = 1 => M = E V^-1,");
        put_line(file," E contains U^-1.");
      end if;
      put_line(file,"The unimodular transformation M :"); put(file,M);
    end if;
  end Unimodular_Coordinate_Transformation;

  procedure Pivots_of_Rational_Coordinate_Transformation
               ( C : in Standard_Integer_Matrices.Matrix;
                 H : out Standard_Integer_Matrices.Matrix;
                 p : out Standard_Integer_Vectors.Vector ) is

    r : integer32;

  begin
    H := Standard_Integer_Matrices.Transpose(C);
    Standard_Integer_Linear_Solvers.Upper_Triangulate(H);
    Standard_Integer_Kernel.Pivots_in_Upper(H,r,p);
  end Pivots_of_Rational_Coordinate_Transformation;

  procedure Pivots_of_Rational_Coordinate_Transformation
               ( file : in file_type;
                 C : in Standard_Integer_Matrices.Matrix;
                 H : out Standard_Integer_Matrices.Matrix;
                 p : out Standard_Integer_Vectors.Vector ) is

    r : integer32;

  begin
    H := Standard_Integer_Matrices.Transpose(C);
    put_line(file,"The matrix "); put(file,H);
    Standard_Integer_Linear_Solvers.Upper_Triangulate(H);
    put_line("has the Hermite Normal Form equal to"); put(file,H);
    Standard_Integer_Kernel.Pivots_in_Upper(H,r,p);
    put(file,"Pivots of Hermite form : ");
    for i in p'range loop
      put(file," "); put(file,H(i,p(i)),1);
    end loop;
    new_line(file);
  end Pivots_of_Rational_Coordinate_Transformation;

  function Is_Zero_Row
              ( A : Standard_Integer_Matrices.Matrix; i : integer32 )
              return boolean is
  begin
    for j in A'range(2) loop
      if A(i,j) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero_Row;

  function Is_Zero_Row
              ( A : Standard_Integer64_Matrices.Matrix; i : integer32 )
              return boolean is
  begin
    for j in A'range(2) loop
      if A(i,j) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero_Row;

  procedure Test_Unimodular_Coordinate_Transformation
               ( A,M : in Standard_Integer_Matrices.Matrix;
                 dim : in integer32; nzr : out natural32;
                 fail : out boolean ) is

    use Standard_Integer_Matrices;
    p : constant Matrix(M'range(1),A'range(2)) := M*A;

  begin
    fail := false;
    nzr := 0;
    for i in 1..dim loop
      if not Is_Zero_Row(p,i)
       then fail := true;
       else nzr := nzr + 1;
      end if;
    end loop;
    for i in dim+1..p'last(1) loop
      if Is_Zero_Row(p,i)
       then nzr := nzr + 1;
      end if;
    end loop;
  end Test_Unimodular_Coordinate_Transformation;

  procedure Test_Unimodular_Coordinate_Transformation
               ( file : in file_type;
                 A,M : in Standard_Integer_Matrices.Matrix;
                 dim : in integer32; nzr : out natural32;
                 fail : out boolean ) is

    use Standard_Integer_Matrices;
    p : constant Matrix(M'range(1),A'range(2)) := M*A;

  begin
    put(file,"Unimodular transformation M multiplied by exponent matrix A, ");
    put_line(file,"M*A :"); put(file,p);
    fail := false;
    nzr := 0;
    for i in 1..dim loop
      if not Is_Zero_Row(p,i) then
        put(file,"Row "); put(file,i,1);
        put_line(file," of M*A is nonzero, bug!");
        fail := true;
      else
        nzr := nzr + 1;
      end if;
    end loop;
    for i in dim+1..p'last(1) loop
      if Is_Zero_Row(p,i)
       then nzr := nzr + 1;
      end if;
    end loop;
    put(file,"The number of nonzero rows in M*A : ");
    put(file,nzr,1); new_line(file);
  end Test_Unimodular_Coordinate_Transformation;

  procedure Test_Unimodular_Coordinate_Transformation
               ( A,M : in Standard_Integer64_Matrices.Matrix;
                 dim : in integer32; nzr : out natural32;
                 fail : out boolean ) is

    use Standard_Integer64_Matrices;
    p : constant Matrix(M'range(1),A'range(2)) := M*A;

  begin
    fail := false;
    nzr := 0;
    for i in 1..dim loop
      if not Is_Zero_Row(p,i)
       then fail := true;
       else nzr := nzr + 1;
      end if;
    end loop;
    for i in dim+1..p'last(1) loop
      if Is_Zero_Row(p,i)
       then nzr := nzr + 1;
      end if;
    end loop;
  end Test_Unimodular_Coordinate_Transformation;

  procedure Test_Unimodular_Coordinate_Transformation
               ( file : in file_type;
                 A,M : in Standard_Integer64_Matrices.Matrix;
                 dim : in integer32; nzr : out natural32;
                 fail : out boolean ) is

    use Standard_Integer64_Matrices;
    p : constant Matrix(M'range(1),A'range(2)) := M*A;

  begin
    put(file,"Unimodular transformation M multiplied by exponent matrix A, ");
    put_line(file,"M*A :"); put(file,p);
    fail := false;
    nzr := 0;
    for i in 1..dim loop
      if not Is_Zero_Row(p,i) then
        put(file,"Row "); put(file,i,1);
        put_line(file," of M*A is nonzero, bug!");
        fail := true;
      else
        nzr := nzr + 1;
      end if;
    end loop;
    for i in dim+1..p'last(1) loop
      if Is_Zero_Row(p,i)
       then nzr := nzr + 1;
      end if;
    end loop;
    put(file,"The number of nonzero rows in M*A : ");
    put(file,nzr,1); new_line(file);
  end Test_Unimodular_Coordinate_Transformation;

  function Rational_Coordinate_Transformation
               ( V : in Standard_Integer_Matrices.Matrix;
                 pivots : in Standard_Integer_Vectors.Vector )
               return Standard_Integer_Matrices.Matrix is

    res : Standard_Integer_Matrices.Matrix(V'range(1),V'range(1))
        := Standard_Smith_Normal_Form.Identity(natural32(V'last(1)));
    tmp : integer32;

  begin
    for j in V'range(2) loop     -- copy j-th tropism into res
      for i in V'range(1) loop   -- placing it on row j
        res(j,i) := V(i,j);
      end loop;
    end loop;
    for j in pivots'range loop   -- take pivots into account
      if pivots(j) /= j then     -- swap columns of added basis vectors
        for i in V'last(2)+1..V'last(1) loop
          tmp := res(i,j);
          res(i,j) := res(i,pivots(j));
          res(i,pivots(j)) := tmp;
        end loop;
      end if;
    end loop;
    return res;
  end Rational_Coordinate_Transformation;

  procedure Test_Rational_Coordinate_Transformation
              ( file : in file_type;
                M : in Standard_Integer_Matrices.Matrix;
                w : in Standard_Integer_Vectors.Vector;
                tol : in double_float; fail : out boolean ) is

    fM : Standard_Floating_Matrices.Matrix(M'range(1),M'range(2));
    det : double_float;

  begin
    for i in w'range loop
      for j in M'range(2) loop
        fM(i,j) := double_float(M(i,j))/double_float(w(i));
      end loop;
    end loop;
    for i in w'last(1)+1..M'last(1) loop
      for j in M'range(2) loop
        fM(i,j) := double_float(M(i,j));
      end loop;
    end loop;
    put_line(file,"The rational coordinate transformation : ");
    put(file,fM,3);
    det := Evaluated_Minors.Determinant(fM);
    put("the determinant : "); put(det);
    fail := (abs(abs(det) - 1.0) > tol);
    if not fail
     then put_line("  okay");
     else put_line("  BUG?!");
    end if;
  end Test_Rational_Coordinate_Transformation;

end Standard_Exponent_Transformations;
