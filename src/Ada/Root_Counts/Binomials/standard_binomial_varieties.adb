with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Norms_Equals;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices_io;    use Standard_Integer64_Matrices_io;
with Standard_Random_Vectors;
with Standard_Integer_Linear_Solvers;
with Standard_Integer64_Linear_Solvers;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Integer_Kernel;
with Standard_Integer64_Kernel;
with Standard_Exponent_Transformations;
with Standard_Binomial_Systems;
with Standard_Binomial_Solvers;
with Standard_Binomial_Varieties_io;
with Standard_Integer32_Triangulations; use Standard_Integer32_Triangulations;
with Standard_Dynamic32_Triangulations; use Standard_Dynamic32_Triangulations;

package body Standard_Binomial_Varieties is

-- STAGE I

  procedure Cone_of_Tropisms 
               ( A : in Standard_Integer_Matrices.Matrix;
                 rank : out integer32;
                 V : out Standard_Integer_Matrices.Link_to_Matrix ) is

    n : constant integer32 := A'last(1);
    k : constant integer32 := A'last(2);
    B : constant Standard_Integer_Matrices.Matrix(1..k,1..n) 
      := Standard_Integer_Matrices.Transpose(A);

  begin
    Standard_Integer_Kernel.Kernel(B,rank,V);
  end Cone_of_Tropisms;

  procedure Cone_of_Tropisms 
               ( A : in Standard_Integer64_Matrices.Matrix;
                 rank : out integer32;
                 V : out Standard_Integer64_Matrices.Link_to_Matrix ) is

    n : constant integer32 := A'last(1);
    k : constant integer32 := A'last(2);
    B : constant Standard_Integer64_Matrices.Matrix(1..k,1..n) 
      := Standard_Integer64_Matrices.Transpose(A);

  begin
    Standard_Integer64_Kernel.Kernel(B,rank,V);
  end Cone_of_Tropisms;

  procedure Check_Inner_Products
               ( A : in Standard_Integer_Matrices.Matrix;
                 V : in Standard_Integer_Matrices.Matrix;
                 bug : out boolean ) is

    s : integer32;

  begin
    bug := false;
    for j in V'range(2) loop    -- all inner products of j-th tropism
      for k in A'range(2) loop  -- with all columns of A
        s := 0;
        for i in A'range(1) loop
          s := s + A(i,k)*V(i,j);
        end loop;
        bug := (s /= 0);
        exit when bug;
      end loop;
    end loop;
  end Check_Inner_Products;

  procedure Check_Inner_Products
               ( file : in file_type;
                 A : in Standard_Integer_Matrices.Matrix;
                 V : in Standard_Integer_Matrices.Matrix;
                 bug : out boolean ) is

    s : integer32;

  begin
    bug := false;
    for j in V'range(2) loop    -- all inner products of j-th tropism
      for k in A'range(2) loop  -- with all columns of A
        s := 0;
        for i in A'range(1) loop
          s := s + A(i,k)*V(i,j);
        end loop;
        put(file," "); put(file,s,1);
        bug := (s /= 0);
        exit when bug;
      end loop;
      new_line(file);
    end loop;
  end Check_Inner_Products;

  procedure Check_Inner_Products
               ( A : in Standard_Integer64_Matrices.Matrix;
                 V : in Standard_Integer64_Matrices.Matrix;
                 bug : out boolean ) is

    s : integer64;

  begin
    bug := false;
    for j in V'range(2) loop    -- all inner products of j-th tropism
      for k in A'range(2) loop  -- with all columns of A
        s := 0;
        for i in A'range(1) loop
          s := s + A(i,k)*V(i,j);
        end loop;
        bug := (s /= 0);
        exit when bug;
      end loop;
    end loop;
  end Check_Inner_Products;

  procedure Check_Inner_Products
               ( file : in file_type;
                 A : in Standard_Integer64_Matrices.Matrix;
                 V : in Standard_Integer64_Matrices.Matrix;
                 bug : out boolean ) is

    s : integer64;

  begin
    bug := false;
    for j in V'range(2) loop    -- all inner products of j-th tropism
      for k in A'range(2) loop  -- with all columns of A
        s := 0;
        for i in A'range(1) loop
          s := s + A(i,k)*V(i,j);
        end loop;
        put(file," "); put(file,s,1);
        bug := (s /= 0);
        exit when bug;
      end loop;
      new_line(file);
    end loop;
  end Check_Inner_Products;

  procedure Check_Rank
               ( V : in Standard_Integer_Matrices.Matrix;
                 rank : in integer32; bug : out boolean ) is

    n : constant integer32 := V'last(1);
    d : constant integer32 := V'last(2);
    r : integer32;

  begin
    if rank /= n-d then
      bug := true;
    else
      r := Standard_Integer_Linear_Solvers.Rank(V);
      bug := not (r = d);
    end if;
  end Check_Rank;

  procedure Check_Rank
               ( file : in file_type;
                 V : in Standard_Integer_Matrices.Matrix;
                 rank : in integer32; bug : out boolean ) is

    n : constant integer32 := V'last(1);
    d : constant integer32 := V'last(2);
    r : integer32;

  begin
    if rank /= n-d then
      bug := true;
    else
      r := Standard_Integer_Linear_Solvers.Rank(V);
      put(file,"Rank of tropisms cone : ");
      put(file,r,1); new_line(file);
      bug := not (r = d);
    end if;
  end Check_Rank;

  procedure Check_Rank
               ( V : in Standard_Integer64_Matrices.Matrix;
                 rank : in integer32; bug : out boolean ) is

    n : constant integer32 := V'last(1);
    d : constant integer32 := V'last(2);
    r : integer32;

  begin
    if rank /= n-d then
      bug := true;
    else
      r := Standard_Integer64_Linear_Solvers.Rank(V);
      bug := not (r = d);
    end if;
  end Check_Rank;

  procedure Check_Rank
               ( file : in file_type;
                 V : in Standard_Integer64_Matrices.Matrix;
                 rank : in integer32; bug : out boolean ) is

    n : constant integer32 := V'last(1);
    d : constant integer32 := V'last(2);
    r : integer32;

  begin
    if rank /= n-d then
      bug := true;
    else
      r := Standard_Integer64_Linear_Solvers.Rank(V);
      put(file,"Rank of tropisms cone : ");
      put(file,r,1); new_line(file);
      bug := not (r = d);
    end if;
  end Check_Rank;

  procedure Check_Cone 
               ( A : in Standard_Integer_Matrices.Matrix;
                 V : in Standard_Integer_Matrices.Matrix;
                 rank : in integer32; bug : out boolean ) is
  begin
    Check_Inner_Products(A,V,bug);
    if not bug
     then Check_Rank(V,rank,bug);
    end if;
  end Check_Cone;

  procedure Check_Cone 
               ( file : in file_type;
                 A : in Standard_Integer_Matrices.Matrix;
                 V : in Standard_Integer_Matrices.Matrix;
                 rank : in integer32; bug : out boolean ) is
  begin
    Check_Inner_Products(file,A,V,bug);
    if not bug
     then Check_Rank(file,V,rank,bug);
    end if;
  end Check_Cone;

  procedure Check_Cone 
               ( A : in Standard_Integer64_Matrices.Matrix;
                 V : in Standard_Integer64_Matrices.Matrix;
                 rank : in integer32; bug : out boolean ) is
  begin
    Check_Inner_Products(A,V,bug);
    if not bug
     then Check_Rank(V,rank,bug);
    end if;
  end Check_Cone;

  procedure Check_Cone 
               ( file : in file_type;
                 A : in Standard_Integer64_Matrices.Matrix;
                 V : in Standard_Integer64_Matrices.Matrix;
                 rank : in integer32; bug : out boolean ) is
  begin
    Check_Inner_Products(file,A,V,bug);
    if not bug
     then Check_Rank(file,V,rank,bug);
    end if;
  end Check_Cone;

  procedure Expected_Dimension
               ( A : in Standard_Integer_Matrices.Matrix;
                 rnk : in integer32; dim : out integer32 ) is

    n : constant integer32 := A'last(1);
    d : integer32;

  begin
    d := n - rnk;
    if d > 0
     then dim := d;
     else dim := 0;
    end if;
  end Expected_Dimension;

  procedure Expected_Dimension
               ( file : in file_type;
                 A,V : in Standard_Integer_Matrices.Matrix;
                 rnk : in integer32; dim : out integer32 ) is

    n : constant integer32 := A'last(1);
    k : constant integer32 := A'last(2);
    d,m : integer32;

  begin
    d := n - rnk;
    if rnk = k then
      put_line(file,"The rank equals the codimension.");
      put(file,"=> For nonzero coefficients, we have a ");
      put(file,d,1); put_line(file,"-dimensional solution set.");
    else
      put_line(file,"The rank does not equal the codimension.");
    end if;
    if d > 0 then
      put_line(file,"The cone of tropisms is spanned by"); put(file,V); 
      m := V'last(2);
      put(file,"The expected dimension is "); put(file,d,1);
      put(file,", found "); put(file,m,1); put_line(file," tropisms.");
      dim := d;
    else
      dim := 0;
    end if;
  end Expected_Dimension;

  procedure Expected_Dimension
               ( A : in Standard_Integer64_Matrices.Matrix;
                 rnk : in integer32; dim : out integer32 ) is

    n : constant integer32 := A'last(1);
    d : integer32;

  begin
    d := n - rnk;
    if d > 0
     then dim := d;
     else dim := 0;
    end if;
  end Expected_Dimension;

  procedure Expected_Dimension
               ( file : in file_type;
                 A,V : in Standard_Integer64_Matrices.Matrix;
                 rnk : in integer32; dim : out integer32 ) is

    n : constant integer32 := A'last(1);
    k : constant integer32 := A'last(2);
    d,m : integer32;

  begin
    d := n - rnk;
    if rnk = k then
      put_line(file,"The rank equals the codimension.");
      put(file,"=> For nonzero coefficients, we have a ");
      put(file,d,1);
      put_line(file,"-dimensional solution set.");
    else
      put_line(file,"The rank does not equal the codimension.");
    end if;
    if d > 0 then
      put_line(file,"The cone of tropisms is spanned by"); put(file,V); 
      m := V'last(2);
      put(file,"The expected dimension is "); put(file,d,1);
      put(file,", found "); put(file,m,1); put_line(file," tropisms.");
      dim := d;
    else
      dim := 0;
    end if;
  end Expected_Dimension;

-- STAGE II, see the package standard_exponent_transformations.

-- STAGE III

  procedure Upper_Transformed_Exponents
               ( A,M : in Standard_Integer_Matrices.Matrix;
                 dim : in integer32;
                 U : out Standard_Integer_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector ) is

    use Standard_Integer_Matrices;
    p : constant Matrix(M'range(1),A'range(2)) := M*A;
    W : Matrix(1..M'last(1)-dim,A'range(2));

  begin
    for i in W'range(1) loop
      for j in W'range(2) loop
        W(i,j) := p(i+dim,j);
      end loop;
    end loop;
    Standard_Integer_Linear_Solvers.Upper_Triangulate(W);
    Standard_Integer_Kernel.Pivots_in_Upper(W,rank,pivots);
    U := W;
  end Upper_Transformed_Exponents;

  procedure Upper_Transformed_Exponents
               ( file : in file_type;
                 A,M : in Standard_Integer_Matrices.Matrix;
                 dim : in integer32;
                 U : out Standard_Integer_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector ) is

    use Standard_Integer_Matrices;
    p : constant Matrix(M'range(1),A'range(2)) := M*A;
    W : Matrix(1..M'last(1)-dim,A'range(2));

  begin
    for i in W'range(1) loop
      for j in W'range(2) loop
        W(i,j) := p(i+dim,j);
      end loop;
    end loop;
    Standard_Integer_Linear_Solvers.Upper_Triangulate(W);
    put_line(file,"The upper triangulated M*A : "); put(file,W);
    Standard_Integer_Kernel.Pivots_in_Upper(W,rank,pivots);
    put(file,"The rank of the matrix M*A : ");
    put(file,rank,1); new_line(file);
    put(file,"The pivots : "); put(file,pivots); new_line(file);
    U := W;
  end Upper_Transformed_Exponents;

  procedure Upper_Transformed_Exponents
               ( A,M : in Standard_Integer64_Matrices.Matrix;
                 dim : in integer32;
                 U : out Standard_Integer64_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector ) is

    use Standard_Integer64_Matrices;
    p : constant Matrix(M'range(1),A'range(2)) := M*A;
    W : Matrix(1..M'last(1)-dim,A'range(2));

  begin
    for i in W'range(1) loop
      for j in W'range(2) loop
        W(i,j) := p(i+dim,j);
      end loop;
    end loop;
    Standard_Integer64_Linear_Solvers.Upper_Triangulate(W);
    Standard_Integer64_Kernel.Pivots_in_Upper(W,rank,pivots);
    U := W;
  end Upper_Transformed_Exponents;

  procedure Upper_Transformed_Exponents
               ( file : in file_type;
                 A,M : in Standard_Integer64_Matrices.Matrix;
                 dim : in integer32;
                 U : out Standard_Integer64_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector ) is

    use Standard_Integer64_Matrices;
    p : constant Matrix(M'range(1),A'range(2)) := M*A;
    W : Matrix(1..M'last(1)-dim,A'range(2));

  begin
    for i in W'range(1) loop
      for j in W'range(2) loop
        W(i,j) := p(i+dim,j);
      end loop;
    end loop;
    Standard_Integer64_Linear_Solvers.Upper_Triangulate(W);
    put_line(file,"The upper triangulated M*A : "); put(file,W);
    Standard_Integer64_Kernel.Pivots_in_Upper(W,rank,pivots);
    put(file,"The rank of the matrix M*A : ");
    put(file,rank,1); new_line(file);
    put(file,"The pivots : "); put(file,pivots); new_line(file); 
    U := W;
  end Upper_Transformed_Exponents;

  procedure Nonpivot_Selection
               ( A : in Standard_Integer_Matrices.Matrix;
                 p : in Standard_Integer_Vectors.Vector;
                 B : out Standard_Integer_Matrices.Matrix ) is

    kp,kB : integer32;

  begin
    for j in A'range(2) loop
      kp := p'first;
      kB := B'first(1)-1;
      for i in A'range(1) loop
        if kp <= p'last and then i = p(kp)
         then kp := kp + 1;
         else kB := kB + 1; B(kB,j) := A(i,j);
        end if;
      end loop;
    end loop;
  end Nonpivot_Selection;

  procedure Nonpivot_Selection
               ( file : in file_type;
                 A : in Standard_Integer_Matrices.Matrix;
                 p : in Standard_Integer_Vectors.Vector;
                 B : out Standard_Integer_Matrices.Matrix ) is
  begin
    put_line(file,"the exponent matrix A : "); put(file,A);
    put(file,"the pivots : "); put(file,p); new_line(file);
    Nonpivot_Selection(A,p,B);
    put_line(file,"after nonpivot selection : "); put(file,B);
  end Nonpivot_Selection;

  function Product_of_Pivots
              ( U : Standard_Integer_Matrices.Matrix;
                pivots : Standard_Integer_Vectors.Vector ) 
              return integer32 is

    res : integer32 := 1;

  begin
    for i in pivots'range loop
      exit when i > U'last(1);
      res := res*U(i,pivots(i));
    end loop;
    return res;
  end Product_of_Pivots;

  function Product_of_Pivots
              ( U : Standard_Integer64_Matrices.Matrix;
                pivots : Standard_Integer_Vectors.Vector ) 
              return Standard_Integer_Numbers.integer64 is

    res : integer64 := 1;

  begin
    for i in pivots'range loop
      exit when i > U'last(1);
      res := res*U(i,pivots(i));
    end loop;
    return res;
  end Product_of_Pivots;

  procedure Solve_Leading_Coefficient_System
               ( n,d : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 pivots : Standard_Integer_Vectors.Vector;
                 c : in Standard_Complex_Vectors.Vector;
                 sols : out Solution_List ) is

    dim : constant integer32 := n-d;
    cp : Standard_Complex_Vectors.Vector(1..dim);
    Ap : Standard_Integer_Matrices.Matrix(1..dim,1..dim);
    M,U : Standard_Integer_Matrices.Matrix(1..dim,1..dim);
    r : integer32;

  begin
    for i in 1..dim loop
      cp(i) := c(pivots(i));
      for j in 1..dim loop
        Ap(i,j) := A(i,pivots(j));
      end loop;
    end loop;
    Standard_Binomial_Solvers.Solve(Ap,cp,r,M,U,sols);
  end Solve_Leading_Coefficient_System;

  procedure Solve_Leading_Coefficient_System
               ( file : in file_type; n,d : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 pivots : Standard_Integer_Vectors.Vector;
                 c : in Standard_Complex_Vectors.Vector;
                 sols : out Solution_List ) is

    dim : constant integer32 := n-d;
    cp : Standard_Complex_Vectors.Vector(1..dim);
    Ap : Standard_Integer_Matrices.Matrix(1..dim,1..dim);
    M,U : Standard_Integer_Matrices.Matrix(1..dim,1..dim);
    r : integer32;

  begin
    for i in 1..dim loop
      cp(i) := c(pivots(i));
      for j in 1..dim loop
        Ap(i,j) := A(i,pivots(j));
      end loop;
    end loop;
    put_line(file,"Ap : "); put(file,Ap);
    Standard_Binomial_Solvers.Solve(Ap,cp,r,M,U,sols);
    put(file,"Found "); put(file,Length_Of(sols),1);
    put_line(file," leading coefficients.");
    put_line(file,"The solutions :"); put(file,sols);
    put_line(file,"The residuals : "); 
    Standard_Binomial_Solvers.Write_Residuals(file,Ap,cp,sols);
  end Solve_Leading_Coefficient_System;

  procedure Solve_Projected_Coefficient_System 
               ( n : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector;
                 sols : out Solution_List ) is

  -- NOTE : Unlike the Solve_Leading_Coefficient_System,
  --   this solver is set up to work only for square systems.

    M,U : Standard_Integer_Matrices.Matrix(1..n,1..n);
    r : integer32;

  begin
    Standard_Binomial_Solvers.Solve(A,c,r,M,U,sols);
  end Solve_Projected_Coefficient_System;

  procedure Solve_Projected_Coefficient_System 
               ( file : in file_type; n : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector;
                 sols : out Solution_List ) is

  -- NOTE : Unlike the Solve_Leading_Coefficient_System,
  --   this solver is set up to work only for square systems.

    M,U : Standard_Integer_Matrices.Matrix(1..n,1..n);
    r : integer32;

  begin
    Standard_Binomial_Solvers.Solve(A,c,r,M,U,sols);
    put(file,"Found "); put(file,Length_Of(sols),1);
    put_line(file," leading coefficients.");
    put_line(file,"The solutions :"); put(file,sols);
    put_line(file,"The residuals : "); 
    Standard_Binomial_Solvers.Write_Residuals(file,A,c,sols);
  end Solve_Projected_Coefficient_System;

-- STAGE IV

  function Evaluate_Tropisms
               ( T : Standard_Integer_Matrices.Matrix;
                 z : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(T'range(1));
 
  begin
    for i in res'range loop
      res(i) := Create(1.0);
    end loop;
    for i in res'range loop
      for j in z'range loop
        if T(i,j) > 0 then
          for k in 1..T(i,j) loop
            res(i) := res(i)*z(j);
          end loop;
        elsif T(i,j) < 0 then
          for k in 1..(-T(i,j)) loop
            res(i) := res(i)/z(j);
          end loop;
        end if;
      end loop;
    end loop;
    return res;
  end Evaluate_Tropisms;

  function Transform_Coefficients
               ( d : integer32; M : Standard_Integer_Matrices.Matrix;
                 c : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(M'range(1));
 
  begin
    for j in M'range(2) loop
      res(j) := Create(1.0);
      for i in d+1..M'last(1) loop
        if M(i,j) > 0 then
          for k in 1..M(i,j) loop
            res(j) := res(j)*c(i-d);
          end loop;
        elsif M(i,j) < 0 then
          for k in 1..(-M(i,j)) loop
            res(j) := res(j)/c(i-d);
          end loop;
        end if;
      end loop;
    end loop;
    return res;
  end Transform_Coefficients;

  function Evaluate_Algebraic_Set
               ( M : Standard_Integer_Matrices.Matrix;
                 c,z : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(M'range(1));
    d : constant integer32 := z'last;
    T : Standard_Integer_Matrices.Matrix(M'range(1),z'range);
    x,y : Standard_Complex_Vectors.Vector(M'range(1));

  begin
    for i in 1..d loop
      for j in T'range(1) loop
        T(j,i) := M(i,j);
      end loop;
    end loop;
    x := Evaluate_Tropisms(T,z);
    y := Transform_Coefficients(d,M,c);
    for i in res'range loop
      res(i) := x(i)*y(i);
    end loop;
    return res;
  end Evaluate_Algebraic_Set;

  function Evaluate_Binomial_System
               ( A : Standard_Integer_Matrices.Matrix;
                 b,x : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(b'range);

  begin
    for j in res'range loop
      res(j) := Create(1.0);
      for i in x'range loop
        if A(i,j) > 0 then
          for k in 1..A(i,j) loop
            res(j) := res(j)*x(i);
          end loop;
        elsif A(i,j) < 0 then
          for k in 1..(-A(i,j)) loop
            res(j) := res(j)/x(i);
          end loop;
        end if;
      end loop;     
      res(j) := res(j) - b(j);
    end loop;
    return res;
  end Evaluate_Binomial_System;

  procedure Residual_Test
               ( A,M : in Standard_Integer_Matrices.Matrix;
                 b,c,z : in Standard_Complex_Vectors.Vector;
                 y : out Standard_Complex_Vectors.Vector;
                 r : out double_float ) is

    x : constant Standard_Complex_Vectors.Vector
      := Evaluate_Algebraic_Set(M,c,z);

  begin
    y := Evaluate_Binomial_System(A,b,x);
    r := Standard_Complex_Norms_Equals.Max_Norm(y);
  end Residual_Test;

  procedure Random_Point_Filter
               ( A : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 d : in integer32; M : in Standard_Integer_Matrices.Matrix;
                 c : in Solution_List; tol : in double_float;
                 filtered : out Solution_List ) is

    tmp : Solution_List := c;
    z : constant Standard_Complex_Vectors.Vector(1..d)
      := Standard_Random_Vectors.Random_Vector(1,d);
    y : Standard_Complex_Vectors.Vector(b'range);
    r : double_float;
    cnt : natural32 := 0;
    true_c,true_c_last : Solution_List;
 
  begin
    Standard_Binomial_Varieties_io.Fill_Symbol_Table(natural32(M'last(1)));
    for i in 1..Length_Of(c) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        cf : constant Standard_Complex_Vectors.Vector := ls.v;
      begin
        Residual_Test(A,M,b,cf,z,y,r);
        if r < tol then
          cnt := cnt + 1;
          Append(true_c,true_c_last,ls.all);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    filtered := true_c;
  end Random_Point_Filter;

  procedure Random_Point_Filter
               ( file : in file_type;
                 A : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 d : in integer32; M : in Standard_Integer_Matrices.Matrix;
                 c : in Solution_List; tol : in double_float;
                 filtered : out Solution_List ) is

    tmp : Solution_List := c;
    z : constant Standard_Complex_Vectors.Vector(1..d)
      := Standard_Random_Vectors.Random_Vector(1,d);
    y : Standard_Complex_Vectors.Vector(b'range);
    r : double_float;
    cnt : natural32 := 0;
    true_c,true_c_last : Solution_List;
 
  begin
    Standard_Binomial_Varieties_io.Fill_Symbol_Table(natural32(M'last(1)));
    for i in 1..Length_Of(c) loop
      put(file,"testing solution "); put(file,i,1); put_line(file," ...");
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        cf : constant Standard_Complex_Vectors.Vector := ls.v;
      begin
        Residual_Test(A,M,b,cf,z,y,r);
        put(file,"The max norm : "); put(file,r); new_line(file);
        if r < tol then
          cnt := cnt + 1;
          Standard_Binomial_Varieties_io.Write_Header
            (file,natural32(M'last(1)),natural32(d));
          Standard_Binomial_Varieties_io.Write_Solution
            (file,natural32(d),M,cf);
          Append(true_c,true_c_last,ls.all);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"Found "); put(file,cnt,1); put_line(file," true solutions.");
    filtered := true_c;
  end Random_Point_Filter;

-- THE SOLVERS :

  procedure Solve_with_Integer_Coordinate_Transformation
               ( d : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 tol : in double_float;
                 M : in Standard_Integer_Matrices.Matrix;
                 c : out Solution_List ) is

  -- DESCRIPTION :
  --   Auxiliary routine to Solve in case an integer unimodular
  --   coordinate transformation is defined by M.

    U : Standard_Integer_Matrices.Matrix(1..M'last(1)-d,A'range(2));
    fail : boolean;
    nzr : natural32;
    rank : integer32;
    pivots : Standard_Integer_Vectors.Vector(b'range);
    cff : Solution_List;
    use Standard_Integer_Matrices,Standard_Binomial_Varieties;
    use Standard_Exponent_Transformations;

  begin
    Test_Unimodular_Coordinate_Transformation(A,M,d,nzr,fail);
    pivots := (pivots'range => 0);
    Upper_Transformed_Exponents(A,M,d,U,rank,pivots);
    declare
      p : constant Matrix := M*A;
      q : Matrix(1..p'last(1)-d,p'range(2));
      n : constant integer32 := A'last(1);
    begin
      for i in d+1..p'last(1) loop
        for j in p'range(2) loop
          q(i-d,j) := p(i,j);
        end loop;
      end loop;
      Solve_Leading_Coefficient_System(n,d,q,pivots,b,cff);
    end;
    Random_Point_Filter(A,b,d,M,cff,tol,c);
  end Solve_with_Integer_Coordinate_Transformation;

  procedure Solve_with_Integer_Coordinate_Transformation
               ( file : in file_type; d : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 tol : in double_float;
                 M : in Standard_Integer_Matrices.Matrix;
                 c : out Solution_List ) is

  -- DESCRIPTION :
  --   Auxiliary routine to Solve in case an integer unimodular
  --   coordinate transformation is defined by M.

    U : Standard_Integer_Matrices.Matrix(1..M'last(1)-d,A'range(2));
    fail : boolean;
    nzr : natural32;
    rank : integer32;
    pivots : Standard_Integer_Vectors.Vector(b'range);
    pp : integer32;
    cff : Solution_List;
    use Standard_Integer_Matrices,Standard_Exponent_Transformations;

  begin
    Test_Unimodular_Coordinate_Transformation(file,A,M,d,nzr,fail);
    pivots := (pivots'range => 0);
    Upper_Transformed_Exponents(file,A,M,d,U,rank,pivots);
    pp := Product_of_Pivots(U,pivots);
    put(file,"The product of pivots : ");
    put(file,pp,1); new_line(file);
    declare
      p : constant Matrix := M*A;
      q : Matrix(1..p'last(1)-d,p'range(2));
      n : constant integer32 := A'last(1);
    begin
      for i in d+1..p'last(1) loop
        for j in p'range(2) loop
          q(i-d,j) := p(i,j);
        end loop;
      end loop;
      Solve_Leading_Coefficient_System(file,n,d,q,pivots,b,cff);
    end;
    Random_Point_Filter(file,A,b,d,M,cff,tol,c);
  end Solve_with_Integer_Coordinate_Transformation;

  procedure Solve
               ( d : in integer32;
                 A,V : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 tol : in double_float;
                 M : out Standard_Integer_Matrices.Matrix;
                 w : out Standard_Integer_Vectors.Vector;
                 c : out Solution_List ) is

    fail : boolean;
    use Standard_Integer_Matrices,Standard_Exponent_Transformations;

  begin
    Unimodular_Coordinate_Transformation(V,fail,M);
    if fail then
      declare
        H : Standard_Integer_Matrices.Matrix(V'range(2),V'range(1));
        Hpiv : Standard_Integer_Vectors.Vector(H'range(1));
        HB : Matrix(A'first(1)..A'last(1)-Hpiv'last,A'range(2));
      begin
        Pivots_of_Rational_Coordinate_Transformation(V,H,Hpiv);
        for i in Hpiv'range loop
          w(i) := H(i,Hpiv(i));
        end loop;
        Nonpivot_Selection(A,Hpiv,HB);
        Solve_Projected_Coefficient_System(HB'last(1),HB,b,c);
        M := Rational_Coordinate_Transformation(V,Hpiv);
      end;
    else
      Solve_with_Integer_Coordinate_Transformation(d,A,b,tol,M,c);
      w := (w'range => 0);
    end if;
  end Solve;

  procedure Solve
               ( file : in file_type; d : in integer32;
                 A,V : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 tol : in double_float;
                 M : out Standard_Integer_Matrices.Matrix;
                 w : out Standard_Integer_Vectors.Vector;
                 c : out Solution_List ) is

    fail : boolean;
    use Standard_Integer_Matrices,Standard_Exponent_Transformations;

  begin
    Unimodular_Coordinate_Transformation(file,V,fail,M);
    if fail then
      declare
        H : Standard_Integer_Matrices.Matrix(V'range(2),V'range(1));
        Hpiv : Standard_Integer_Vectors.Vector(H'range(1));
        HB : Matrix(A'first(1)..A'last(1)-Hpiv'last,A'range(2));
      begin
        Pivots_of_Rational_Coordinate_Transformation(file,V,H,Hpiv);
        for i in Hpiv'range loop
          w(i) := H(i,Hpiv(i));
        end loop;
        Nonpivot_Selection(file,A,Hpiv,HB);
        Solve_Projected_Coefficient_System(file,HB'last(1),HB,b,c);
        M := Rational_Coordinate_Transformation(V,Hpiv);
        Test_Rational_Coordinate_Transformation(file,M,w,1.0E-8,fail);
      end;
    else
      Solve_with_Integer_Coordinate_Transformation(file,d,A,b,tol,M,c);
      w := (w'range => 0);
    end if;
  end Solve;

  procedure Black_Box_Solver
              ( p : in Laur_Sys;
                fail : out boolean; d : out integer32;
                M : out Standard_Integer_Matrices.Link_to_Matrix;
                c : out Solution_List ) is

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    r : integer32;
    V : Standard_Integer_Matrices.Link_to_Matrix;
 
  begin
    Standard_Binomial_Systems.Parse(p,nq,A,b,fail);
    if not fail then
      Cone_of_Tropisms(A,r,V);
      Expected_Dimension(A,r,d);
      if d > 0 then
        declare
          U : Standard_Integer_Matrices.Matrix(V'range(1),V'range(1));
          w : Standard_Integer_Vectors.Vector(V'range(2));
        begin
          Solve(d,A,V.all,b,1.0E-8,U,w,c);
          M := new Standard_Integer_Matrices.Matrix'(U);
        end;
      end if;
    end if;
  end Black_Box_Solver;

  procedure Black_Box_Solver
              ( file : in file_type; p : in Laur_Sys;
                fail : out boolean; d : out integer32;
                M : out Standard_Integer_Matrices.Link_to_Matrix;
                c : out Solution_List ) is

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    r : integer32;
    V : Standard_Integer_Matrices.Link_to_Matrix;
 
  begin
    Standard_Binomial_Systems.Parse(p,nq,A,b,fail);
    if fail then
      null; -- we do not want the message for phc -b
      -- put_line(file,"The system is not a binomial system!");
    else
      Cone_of_Tropisms(A,r,V);
      Check_Cone(file,A,V.all,r,fail);
      if fail
       then put_line("Bug in computation of cone of tropisms.");
       else put_line("Computation of cone of tropisms is okay.");
      end if;
      Expected_Dimension(file,A,V.all,r,d);
      if d > 0 then
        declare
          U : Standard_Integer_Matrices.Matrix(V'range(1),V'range(1));
          w : Standard_Integer_Vectors.Vector(V'range(2));
        begin
          Solve(file,d,A,V.all,b,1.0E-8,U,w,c);
          M := new Standard_Integer_Matrices.Matrix'(U);
        end;
      end if;
    end if;
  end Black_Box_Solver;

-- COMPUTING THE DEGREE OF A SOLUTION :

  function Support ( T : Standard_Integer_Matrices.Matrix )
                   return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    v : Standard_Integer_Vectors.Vector(T'range(2));

  begin
    for i in T'range(1) loop
      for j in T'range(2) loop
        v(j) := T(i,j);
      end loop;
      Lists_of_Integer_Vectors.Append(res,res_last,v);
    end loop;
    v := (v'range => 0);
    if not Lists_of_Integer_Vectors.Is_In(res,v)
     then Lists_of_Integer_Vectors.Append(res,res_last,v);
    end if;
    return res;
  end Support;

  function Volume ( A : Lists_of_Integer_Vectors.List ) return natural32 is

    res : natural32 := 0;
    lifted,lifted_last : Lists_of_Integer_Vectors.List;
    t : Triangulation;

  begin
    Dynamic_Lifting(A,false,false,0,lifted,lifted_last,t);
    res := Volume(t);
    Lists_of_Integer_Vectors.Clear(lifted);
    Clear(t);
    return res;
  end Volume;

  function Degree ( T : Standard_Integer_Matrices.Matrix ) return natural32 is

    A : Lists_of_Integer_Vectors.List := Support(T);
    v : constant natural32 := Volume(A);

  begin
    Lists_of_Integer_Vectors.Clear(A);
    return v;
  end Degree;

end Standard_Binomial_Varieties;
