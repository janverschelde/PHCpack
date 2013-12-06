with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Random_Matrices;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Standard_Floating_Norms_Equals;     use Standard_Floating_Norms_Equals;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Floating_Vectors;
with Multprec_Floating_Vectors_io;
with Multprec_Floating_Matrices;
with Multprec_Floating_Matrices_io;
with Multprec_Floating_Linear_Solvers;   use Multprec_Floating_Linear_Solvers;
with Multprec_Floating_Norms_Equals;     use Multprec_Floating_Norms_Equals;
with Test_LU_Decompositions;             use Test_LU_Decompositions;

procedure ts_fltlu is

-- DESCRIPTION :
--   Tests LU factorization of standard and multi-precision floats.

  function Vdm_Matrix ( v : Multprec_Floating_Vectors.Vector )
                      return Multprec_Floating_Matrices.Matrix is

    use Multprec_Floating_Matrices;

    n : constant integer32 := v'length;
    res : Matrix(1..n,1..n);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := v(i)**(j-1);
      end loop;
    end loop;
    return res;
  end Vdm_Matrix;

  function Vdm_Matrix ( v : Standard_Floating_Vectors.Vector )
                      return Standard_Floating_Matrices.Matrix is

    use Standard_Floating_Matrices;

    n : constant integer32 := v'length;
    res : Matrix(1..n,1..n);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := v(i)**integer(j-1);
      end loop;
    end loop;
    return res;
  end Vdm_Matrix;

--  procedure Write ( f : in Multprec_Floating_Numbers.Floating_Number ) is
--
--    frac : constant Integer_Number := Fraction(f);
--    expo : constant Integer_Number := Exponent(f);
--
--  begin
--    put("Fraction :");
--    for i in reverse 0..Size(frac) loop
--      put(Coefficient(frac,i));
--    end loop;
--    new_line;
--    put("Exponent :");
--    for i in reverse 0..Size(expo) loop
--      put(Coefficient(expo,i));
--    end loop;
--    new_line;
--  end Write;

 -- procedure Write ( v : in Multprec_Floating_Vectors.Vector ) is
 -- begin
 --   for i in v'range loop
 --     Write(v(i)); new_line;
 --   end loop;
 -- end Write;

 -- procedure Write ( m : in Multprec_Floating_Matrices.Matrix ) is
 -- begin
 --   for i in m'range(1) loop
 --     for j in m'range(2) loop
 --       Write(m(i,j)); new_line;
 --     end loop;
 --   end loop;
 -- end Write;

  procedure Set_Size ( v : in out Multprec_Floating_Vectors.Vector;
                       k : in natural32 ) is

  -- DESCRIPTION :
  --   Sets the size of the elements in v to k.

    use Multprec_Floating_Vectors;

  begin
    for i in v'range loop
      Set_Size(v(i),k);
    end loop;
    --Write(v);
  end Set_Size;

  procedure Set_Size ( m : in out Multprec_Floating_Matrices.Matrix;
                       k : in natural32 ) is

  -- DESCRIPTION :
  --   Sets the size of the elements in m to k.

    use Multprec_Floating_Matrices;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        Set_Size(m(i,j),k);
      end loop;
    end loop;
    --Write(m);
  end Set_Size;

  procedure lufco_Solve ( n : in integer32;
                          mat : in Standard_Floating_Matrices.Matrix;
                          rhs : in Standard_Floating_Vectors.Vector ) is

    use Standard_Floating_Vectors;
    use Standard_Floating_Vectors_io;
    use Standard_Floating_Matrices;

    wrk : Matrix(mat'range(1),mat'range(2)) := mat;
    piv : Standard_Integer_Vectors.Vector(mat'range(2));
    rcond,nrm : double_float;
    res,sol : Vector(rhs'range);

  begin
    put_line("Solving the linear system with lufco.");
    lufco(wrk,n,piv,rcond);
    put("inverse condition : "); put(rcond); new_line;
    sol := rhs;
    lusolve(wrk,n,piv,sol);
    put_line("The solution vector :"); put(sol); new_line;
    res := rhs - mat*sol;
    put_line("The residual : "); put(res); new_line;
    nrm := Max_Norm(res);
    put("Max norm of residual : "); put(nrm); new_line;
    nrm := Sum_Norm(res);
    put("Sum norm of residual : "); put(nrm); new_line;
  end lufco_Solve;

  procedure lufac_Solve ( n : in integer32;
                          mat : in Standard_Floating_Matrices.Matrix;
                          rhs : in Standard_Floating_Vectors.Vector ) is

    use Standard_Floating_Vectors;
    use Standard_Floating_Vectors_io;
    use Standard_Floating_Matrices;

    wrk : Matrix(mat'range(1),mat'range(2)) := mat;
    piv : Standard_Integer_Vectors.Vector(mat'range(2));
    info : integer32;
    res,sol : Vector(rhs'range);
    nrm : double_float;

  begin
    put_line("Solving the linear system with lufac.");
    lufac(wrk,n,piv,info);
    put("info : "); put(info,1); new_line;
    sol := rhs;
    lusolve(wrk,n,piv,sol);
    put_line("The solution vector :"); put(sol); new_line;
    res := rhs - mat*sol;
    put_line("The residual : "); put(res); new_line;
    nrm := Max_Norm(res);
    put("Max norm of residual : "); put(nrm); new_line;
    nrm := Sum_Norm(res);
    put("Sum norm of residual : "); put(nrm); new_line;
  end lufac_Solve;

  procedure Interactive_Test_Standard_Linear_Solvers is

    use Standard_Floating_Vectors;
    use Standard_Floating_Vectors_io;
    use Standard_Floating_Matrices;

    n : integer32 := 0;

  begin
    new_line;
    put_line("Interactive testing of solving standard linear systems.");
    new_line;
    put("Give the dimension : "); get(n);
    declare
      mat : Matrix(1..n,1..n);
      rhs : Vector(1..n);
    begin
      put("Give "); put(n,1); put("x"); put(n,1);
      put_line(" floating matrix : "); get(mat);
      put_line("-> the matrix : "); put(mat);
      put("Give "); put(n,1); put_line(" floating-numbers : "); get(rhs);
      put_line("-> right-hand side vector : "); put(rhs); new_line;
      lufac_Solve(n,mat,rhs);
      lufco_Solve(n,mat,rhs);
    end;
  end Interactive_Test_Standard_Linear_Solvers;

  procedure Random_Test_Standard_Linear_Solvers is

    use Standard_Floating_Vectors;
    use Standard_Floating_Vectors_io;
    use Standard_Floating_Matrices;

    n,nb : integer32 := 0;
    timer : Timing_Widget;

  begin
    new_line;
    put_line("Testing the solution of random standard linear systems.");
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the number of tests : "); get(nb);
    tstart(timer);
    for i in 1..nb loop
      declare
        mat : Matrix(1..n,1..n);
        rhs : Vector(1..n);
      begin
        mat := Vdm_Matrix(Random_Vector(1,n));
        rhs := Random_Vector(1,n);
        lufac_Solve(n,mat,rhs);
        lufco_Solve(n,mat,rhs);
      end;
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Random_Test_Standard_Linear_Solvers;

  procedure lufco_Solve ( n : in integer32;
                          mat : in Multprec_Floating_Matrices.Matrix;
                          rhs : in Multprec_Floating_Vectors.Vector ) is

    use Multprec_Floating_Vectors;
    use Multprec_Floating_Vectors_io;
    use Multprec_Floating_Matrices;
    use Multprec_Floating_Matrices_io;

    wrk : Matrix(mat'range(1),mat'range(2));
    piv : Standard_Integer_Vectors.Vector(mat'range(2));
    rcond,nrm : Floating_Number;
    res,sol,acc : Vector(rhs'range);

  begin
    put_line("Solving the linear system with lufco.");
    Copy(mat,wrk);
    lufco(wrk,n,piv,rcond);
    put("inverse condition : "); put(rcond); new_line;
    Copy(rhs,sol);
    lusolve(wrk,n,piv,sol);
   -- put_line("The solution vector :"); put(sol); new_line;
    acc := mat*sol;
    res := rhs - acc;
    put_line("The residual : "); put(res); new_line;
    nrm := Max_Norm(res);
    put("Max norm of residual : "); put(nrm); new_line;
    Clear(nrm);
    nrm := Sum_Norm(res);
    put("Sum norm of residual : "); put(nrm); new_line;
    Clear(nrm); Clear(rcond);
    Clear(wrk); Clear(acc); Clear(res); Clear(sol);
  end lufco_Solve;

  procedure lufac_Solve ( n : in integer32;
                          mat : in Multprec_Floating_Matrices.Matrix;
                          rhs : in Multprec_Floating_Vectors.Vector ) is

    use Multprec_Floating_Vectors;
    use Multprec_Floating_Vectors_io;
    use Multprec_Floating_Matrices;
    use Multprec_Floating_Matrices_io;

    wrk : Matrix(mat'range(1),mat'range(2));
    piv : Standard_Integer_Vectors.Vector(mat'range(2));
    info : integer32;
    res,sol,acc : Vector(rhs'range);
    nrm : Floating_Number;

  begin
    put_line("Solving the linear system with lufac.");
    Copy(mat,wrk);
    lufac(wrk,n,piv,info);
    put("info : "); put(info,1); new_line;
    Copy(rhs,sol);
    lusolve(wrk,n,piv,sol);
   -- put_line("The solution vector :"); put(sol); new_line;
    acc := mat*sol;
    res := rhs - acc;
    put_line("The residual : "); put(res); new_line;
    nrm := Max_Norm(res);
    put("Max norm of residual : "); put(nrm); new_line;
    Clear(nrm);
    nrm := Sum_Norm(res);
    put("Sum norm of residual : "); put(nrm); new_line;
    Clear(nrm); Clear(wrk); Clear(acc); Clear(res); Clear(sol);
  end lufac_Solve;

  procedure Interactive_Test_Multprec_Linear_Solvers is

    use Multprec_Floating_Vectors;
    use Multprec_Floating_Vectors_io;
    use Multprec_Floating_Matrices;
    use Multprec_Floating_Matrices_io;

    ans : character;
    n : integer32 := 0;

  begin
    put("Give the dimension : "); get(n);
    declare
      mat : Matrix(1..n,1..n);
      rhs : Vector(1..n);
      sz : integer32 := 0;
    begin
      put("Give "); put(n,1); put("x"); put(n,1);
      put_line(" floating matrix : "); get(mat);
      put_line("-> the matrix : "); put(mat);
      put("Give "); put(n,1); put_line(" floating-numbers : "); get(rhs);
      put_line("-> right-hand side vector : "); put(rhs); new_line;
      put("Give the size (-1 for default) : "); get(sz);
      if sz >= 0 then
        Set_Size(mat,natural32(sz));
        Set_Size(rhs,natural32(sz));
      end if;
      loop
        lufac_Solve(n,mat,rhs);
        lufco_Solve(n,mat,rhs);
        put("Do you want to resolve with other precision ? (y/n) "); get(ans);
        exit when ans /= 'y';
        put("Give the size : "); get(sz);
        Set_Size(mat,natural32(sz));
        Set_Size(rhs,natural32(sz));
      end loop;
    end;
  end Interactive_Test_Multprec_Linear_Solvers;

 procedure Random_Test_Multprec_Linear_Solvers is

    use Multprec_Floating_Vectors;
    use Multprec_Floating_Vectors_io;
    use Multprec_Floating_Matrices;
    use Multprec_Floating_Matrices_io;

    n : integer32 := 0;
    sz,nb : natural32 := 0;
    timer : Timing_Widget;

  begin
    new_line;
    put_line("Testing the solution of random multi-precision linear systems.");
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the size of the numbers : "); get(sz);
    put("Give the number of tests : "); get(nb);
    tstart(timer);
    for i in 1..nb loop
      declare
        mat : Matrix(1..n,1..n);
        rhs : Vector(1..n);
      begin
        rhs := Random_Vector(1,n,sz);
        mat := Vdm_Matrix(rhs); --Random_Matrix(n,n,sz);
        Clear(rhs);
        rhs := Random_Vector(1,n,sz);
       -- lufac_Solve(n,mat,rhs);
        lufco_Solve(n,mat,rhs);
        Clear(mat);
        Clear(rhs);
      end;
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Random_Test_Multprec_Linear_Solvers;

  procedure Standard_Accuracy_Test is

    n : integer32 := 0;

  begin
    new_line;
    put_line("Testing the accuracy of the the LU factorization.");
    new_line;
    put("Give the dimension : "); get(n);
    declare
      ans : character;
      A,LU : Standard_Floating_Matrices.Matrix(1..n,1..n);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      info : integer32;
      tol : constant double_float := 1.0E-8;
      maxerr : double_float;
      fail : boolean;
    begin
      put("Enter your own matrix (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("-> give "); put(n,1); put("-by-"); put(n,1);
        put_line(" matrix :"); get(A);
      else
        put("-> generating a random "); put(n,1); put("-by-"); put(n,1);
        put_line(" matrix ...");
        A := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n));
      end if;
      put_line("The  matrix A : "); put(A);
      LU := A;
      lufac(LU,n,ipvt,info);
      Test_Decomposition(n,A,LU,ipvt,tol,true,maxerr,fail);
      put("largest error :"); put(maxerr,3);
      put("  <"); put(tol,3); put(" ? ");
      if fail
       then put_line("bug!?");
       else put_line("okay.");
      end if;
    end;
  end Standard_Accuracy_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of matrices of floating numbers");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. interactive test on solving standard linear systems.");
      put_line("  2. test on solving random standard linear systems.");
      put_line("  3. interactive accuracy test on standard LU factorization.");
      put_line("  4. interactive test on solving multi-precision systems.");
      put_line("  5. test on solving random multi-precision systems.");
      put("Enter 0, 1, 2, 3, 4, or 5 to make your choice : ");
      Ask_Alternative(ans,"012345");
      exit when ans = '0';
      case ans is
        when '1' => Interactive_Test_Standard_Linear_Solvers;
        when '2' => Random_Test_Standard_Linear_Solvers;
        when '3' => Standard_Accuracy_Test;
        when '4' => Interactive_Test_Multprec_Linear_Solvers;
        when '5' => Random_Test_Multprec_Linear_Solvers;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_fltlu;
