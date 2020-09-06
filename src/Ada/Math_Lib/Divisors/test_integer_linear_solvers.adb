with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;
with Standard_Integer_Linear_Solvers;
with Standard_Integer64_Vectors;
with Standard_Integer64_Vectors_io;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;
with Standard_Integer64_Linear_Solvers;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Multprec_Random_Matrices;           use Multprec_Random_Matrices;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Multprec_Integer_Vectors;
with Multprec_Integer_Vectors_io;
with Multprec_Integer_Matrices_io;
with Multprec_Integer_Linear_Solvers;

package body Test_Integer_Linear_Solvers is

  procedure Random_Test_Standard_Solver
                 ( n,m,low,upp : in integer32 ) is

    use Standard_Integer_Vectors,Standard_Integer_Vectors_io;
    use Standard_Integer_Matrices;
    use Standard_Integer_Linear_Solvers;

    mat : constant Matrix(1..n,1..m)
        := Random_Matrix(natural32(n),natural32(m),low,upp);
    wrk : Matrix(1..n,1..m) := mat;
    sol : Vector(1..m);
    res : Vector(1..n);

  begin
    Upper_Triangulate(wrk);
    Solve0(wrk,sol);
    res := mat*sol;
    put("The residual : "); put(res); 
    put(" of solution : "); put(sol); new_line;
  end Random_Test_Standard_Solver;

  procedure Random_Test_Standard64_Solver
                 ( n,m : in integer32; low,upp : in integer64 ) is

    use Standard_Integer64_Vectors,Standard_Integer64_Vectors_io;
    use Standard_Integer64_Matrices;
    use Standard_Integer64_Linear_Solvers;

    mat : constant Matrix(1..n,1..m)
        := Random_Matrix(natural32(n),natural32(m),low,upp);
    wrk : Matrix(1..n,1..m) := mat;
    sol : Vector(1..m);
    res : Vector(1..n);

  begin
    Upper_Triangulate(wrk);
    Solve0(wrk,sol);
    res := mat*sol;
    put("The residual : "); put(res); 
    put(" of solution : "); put(sol); new_line;
  end Random_Test_Standard64_Solver;

  procedure Random_Test_Multprec_Solver
              ( n,m : in integer32; sz : in natural32 ) is

    use Multprec_Integer_Vectors,Multprec_Integer_Vectors_io;
    use Multprec_Integer_Matrices;
    use Multprec_Integer_Linear_Solvers;

    mat : constant Matrix(1..n,1..m)
        := Random_Matrix(natural32(n),natural32(m),sz);
    wrk : Matrix(1..n,1..m);
    sol : Vector(1..m);
    res : Vector(1..n);

  begin
    Copy(mat,wrk);
    Upper_Triangulate(wrk);
    Solve0(wrk,sol);
    res := mat*sol;
    put("The residual : "); put(res); 
    put(" of solution : "); put(sol); new_line;
  end Random_Test_Multprec_Solver;

  procedure Random_Test_Standard_Solvers is

    n,m,nb,low,upp : integer32 := 0;

  begin
    put("Give number of tests : "); get(nb);
    put("Give number of rows : "); get(n);
    put("Give number of columns : "); get(m);
    put("Give lower bound for numbers : "); get(low);
    put("Give upper bound for numbers : "); get(upp);
    for i in 1..nb loop
      Random_Test_Standard_Solver(n,m,low,upp);
    end loop;
  end Random_Test_Standard_Solvers;

  procedure Random_Test_Standard64_Solvers is

    n,m,nb : integer32 := 0;
    low,upp : integer64 := 0;

  begin
    put("Give number of tests : "); get(nb);
    put("Give number of rows : "); get(n);
    put("Give number of columns : "); get(m);
    put("Give lower bound for numbers : "); get(low);
    put("Give upper bound for numbers : "); get(upp);
    for i in 1..nb loop
      Random_Test_Standard64_Solver(n,m,low,upp);
    end loop;
  end Random_Test_Standard64_Solvers;

  procedure Random_Test_Multprec_Solvers is

    n,m : integer32 := 0;
    nb,sz : natural32 := 0;

  begin
    put("Give number of tests : "); get(nb);
    put("Give number of rows : "); get(n);
    put("Give number of columns : "); get(m);
    put("Give the size of the numbers : "); get(sz);
    for i in 1..nb loop
      Random_Test_Multprec_Solver(n,m,sz);
    end loop;
  end Random_Test_Multprec_Solvers;

  procedure Interactive_Test_Standard_Solvers is

    use Standard_Integer_Vectors,Standard_Integer_Vectors_io;
    use Standard_Integer_Matrices,Standard_Integer_Matrices_io;
    use Standard_Integer_Linear_Solvers;

    n,m : integer32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      sol : Vector(1..m) := (1..m => 0);
      res : Vector(1..n);
      mat,wrk : Matrix(1..n,1..m);
      l,p : Matrix(1..n,1..n);
      d : integer32 := 1;
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" integer matrix : "); get(mat);
      put_line("Your matrix : "); put(mat);
      wrk := mat;
      Upper_Triangulate(l,wrk);
      put_line("The matrix in upper triangular form : "); put(wrk);
      put_line("The transformation matrix T : "); put(l);
      d := Det(l);
      put("The determinant of T is "); put(d,1);
      if d*d /= 1
       then put_line(" BUG!!!!");
       else put_line(" okay");
      end if;
      p := l*mat;
      put_line("Product of T and original matrix : "); put(p);
      if Equal(p,wrk)
       then put_line("Product and triangulated matrix are equal, okay.");
       else put_line("Product and triangulated matrix are NOT equal, BUG!!");
      end if;
      if n = m then
        d := 1;
        for i in 1..n loop
          d := d*wrk(i,i);
        end loop;
        put("the determinant : "); put(d,1); new_line;
      end if;
      Solve0(wrk,sol);
      put_line("The solution of the homogeneous system : ");
      put(sol); new_line;
      res := mat*sol;
      put("The residual : "); put(res); new_line;
    end;
  end Interactive_Test_Standard_Solvers;

  procedure Interactive_Test_Standard64_Solvers is

    use Standard_Integer64_Vectors,Standard_Integer64_Vectors_io;
    use Standard_Integer64_Matrices,Standard_Integer64_Matrices_io;
    use Standard_Integer64_Linear_Solvers;

    n,m : integer32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      sol : Vector(1..m) := (1..m => 0);
      res : Vector(1..n);
      mat,wrk : Matrix(1..n,1..m);
      l,p : Matrix(1..n,1..n);
      d : integer64 := 1;
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" integer matrix : "); get(mat);
      put_line("Your matrix : "); put(mat);
      wrk := mat;
      Upper_Triangulate(l,wrk);
      put_line("The matrix in upper triangular form : "); put(wrk);
      put_line("The transformation matrix T : "); put(l);
      d := Det(l);
      put("The determinant of T is "); put(d,1);
      if d*d /= 1
       then put_line(" BUG!!!!");
       else put_line(" okay");
      end if;
      p := l*mat;
      put_line("Product of T and original matrix : "); put(p);
      if Equal(p,wrk)
       then put_line("Product and triangulated matrix are equal, okay.");
       else put_line("Product and triangulated matrix are NOT equal, BUG!!");
      end if;
      if n = m then
        d := 1;
        for i in 1..n loop
          d := d*wrk(i,i);
        end loop;
        put("the determinant : "); put(d,1); new_line;
      end if;
      Solve0(wrk,sol);
      put_line("The solution of the homogeneous system : ");
      put(sol); new_line;
      res := mat*sol;
      put("The residual : "); put(res); new_line;
    end;
  end Interactive_Test_Standard64_Solvers;

  procedure Show_Differences
              ( A,B : in Multprec_Integer_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if not Equal(A(i,j),B(i,j)) then
          put("A("); put(i,1); put(",");
          put(j,1); put(") = ");
          put(A(i,j)); new_line;
          put("B("); put(i,1); put(",");
          put(j,1); put(") = ");
          put(B(i,j)); new_line;
        end if;
      end loop;
    end loop;
  end Show_Differences;

  procedure Interactive_Test_Multprec_Solvers is

    use Multprec_Integer_Vectors,Multprec_Integer_Vectors_io;
    use Multprec_Integer_Matrices,Multprec_Integer_Matrices_io;
    use Multprec_Integer_Linear_Solvers;

    n,m : integer32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      sol : Vector(1..m);
      res : Vector(1..n);
      mat,wrk : Matrix(1..n,1..m);
      l,p : Matrix(1..n,1..n);
      d : Integer_Number := Create(integer(1));
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" integer matrix : "); get(mat);
      put_line("Your matrix : "); put(mat);
      Copy(mat,wrk);
      Upper_Triangulate(l,wrk);
      put_line("The matrix in upper triangular form : "); put(wrk);
      put_line("The transformation matrix T : "); put(l);
      d := Det(l);
      put("The determinant of T is "); put(d,1);
      if Equal(d,1) then
        put_line(" okay");
      else
        Min(d);
        if not Equal(d,1) 
         then put_line(" BUG!!!!");
         else put_line(" okay");
        end if;
      end if;
      p := l*mat;
      put_line("product of T and original matrix : "); put(p);
      if Equal(p,wrk)
       then put_line("Product and triangulated matrix are equal, okay.");
       else put_line("Product and triangulated matrix are NOT equal, BUG!!");
            Show_Differences(p,wrk);
      end if;
      if n = m then
        d := Create(integer(1));
        for i in 1..n loop
          Mul(d,wrk(i,i));
        end loop;
        put("the determinant : "); put(d); new_line;
      end if;
      Solve0(wrk,sol);
      put_line("The solution of the homogeneous system : ");
      put(sol); new_line;
      res := mat*sol;
      put("The residual : "); put(res); new_line;
    end;
  end Interactive_Test_Multprec_Solvers;

  procedure Main is

    ans,lng : character;

  begin
    new_line;
    put_line("Interactive testing of matrices of integer numbers");
    new_line;
    loop
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. solve given linear systems of standard numbers.");
      put_line("  2. solve random linear systems of standard numbers.");
      put_line("  3. solve given linear systems of multi-precision numbers.");
      put_line("  4. solve random linear systems of multi-precision numbers.");
      put("Type 0, 1, 2, 3, or 4 to choose: ");
      Ask_Alternative(ans,"01234");
      exit when (ans = '0');
      if ans = '1' or ans = '2'
       then put("Use 64-bit integers ? (y/n) "); Ask_Yes_or_No(lng);
      end if;
      new_line;
      case ans is
        when '1' => if lng = 'y'
                     then Interactive_Test_Standard64_Solvers;
                     else Interactive_Test_Standard_Solvers;
                    end if;
        when '2' => if lng = 'y'
                     then Random_Test_Standard64_Solvers;
                     else Random_Test_Standard_Solvers;
                    end if;
        when '3' => Interactive_Test_Multprec_Solvers;
        when '4' => Random_Test_Multprec_Solvers;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Integer_Linear_Solvers;
