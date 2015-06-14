with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Random_Matrices;           use Multprec_Random_Matrices;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;
with Multprec_Complex_Matrices;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Test_LU_Decompositions;

procedure ts_cmplu is

-- DESCRIPTION :
--   Tests LU factorization on standard and multi-precision complex numbers.

  function Vdm_Matrix ( v : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Matrices.Matrix is

    use Standard_Complex_Numbers;
    use Standard_Complex_Matrices;

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

  procedure lufac_Solve ( n : in integer32;
                          mat : in Standard_Complex_Matrices.Matrix;
                          rhs : in Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Vectors;
    use Standard_Complex_Vectors_io;
    use Standard_Complex_Matrices;
    use Test_LU_Decompositions;

    wrk : Matrix(mat'range(1),mat'range(2));
    piv : Standard_Integer_Vectors.Vector(mat'range(2));
    info : integer32;
    res,sol,acc : Vector(rhs'range);
    nrm : double_float;
    tol : constant double_float := 1.0E-8;
    maxerr : double_float;
    fail : boolean;

  begin
    new_line;
    put_line("Solving the linear system with lufac.");
    wrk := mat;
    lufac(wrk,n,piv,info);
    put("info : "); put(info,1); new_line;
    sol := rhs;
    lusolve(wrk,n,piv,sol);
    put_line("The solution vector :"); put(sol); new_line;
    acc := mat*sol;
    res := rhs - acc;
    put_line("The residual : "); put(res); new_line;
    nrm := Max_Norm(res);
    put("Max norm of residual : "); put(nrm); new_line;
    nrm := Sum_Norm(res);
    put("Sum norm of residual : "); put(nrm); new_line;
    new_line;
    put_line("Testing the LU factorization ...");
    new_line;
    Test_Decomposition(n,mat,wrk,piv,tol,true,maxerr,fail);
    put("largest error : "); put(maxerr,3);
    put(" < "); put(tol,3); put(" ? ");
    if fail
     then put_line("bug!?");
     else put_line("okay.");
    end if;
  end lufac_Solve;

  procedure lufco_Solve ( n : in integer32;
                          mat : in Standard_Complex_Matrices.Matrix;
                          rhs : in Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Vectors;
    use Standard_Complex_Vectors_io;
    use Standard_Complex_Matrices;
    use Test_LU_Decompositions;

    wrk : Matrix(mat'range(1),mat'range(2)) := mat;
    piv : Standard_Integer_Vectors.Vector(mat'range(2));
    rcond,nrm,maxerr : double_float;
    res,sol : Vector(rhs'range);
    tol : constant double_float := 1.0E-8;
    fail : boolean;

  begin
    new_line;
    put_line("Solving the linear system with lufco.");
    lufco(wrk,n,piv,rcond);
    put("inverse condition : "); put(rcond); new_line;
    declare
      anorm : constant double_float := Norm1(mat);
      wrk1 : Matrix(mat'range(1),mat'range(2)) := mat;
      ipvt : Standard_Integer_Vectors.Vector(mat'range(2));
      info : integer32;
    begin
      lufac(wrk1,n,ipvt,info);
      if info = 0
       then estco(wrk1,n,ipvt,anorm,rcond);
       else rcond := 0.0;
      end if;
      put_line("computed condition number again with estco");
      put("inverse condition : "); put(rcond); new_line;
    end;
    sol := rhs;
    lusolve(wrk,n,piv,sol);
    put_line("The solution vector :"); put(sol); new_line;
    res := rhs - mat*sol;
    put_line("The residual : "); put(res); new_line;
    nrm := Max_Norm(res);
    put("Max norm of residual : "); put(nrm); new_line;
    nrm := Sum_Norm(res);
    put("Sum norm of residual : "); put(nrm); new_line;
    new_line;
    put_line("Testing the LU factorization ...");
    new_line;
    Test_Decomposition(n,mat,wrk,piv,tol,true,maxerr,fail);
    put("largest error : "); put(maxerr,3);
    put(" < "); put(tol,3); put(" ? ");
    if fail
     then put_line("bug!?");
     else put_line("okay.");
    end if;
  end lufco_Solve;

  procedure Random_Test_Standard_Linear_Solvers is

    use Standard_Complex_Vectors;
    use Standard_Complex_Vectors_io;
    use Standard_Complex_Matrices;

    n,nb : integer32 := 0;
    timer : Timing_Widget;

  begin
    new_line;
    put_line("Testing of solving random standard linear systems.");
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

 -- function Vdm_Matrix ( v : Multprec_Complex_Vectors.Vector )
 --                     return Multprec_Complex_Matrices.Matrix is

  -- NOTE : this function causes a memory leak !!!

 --   use Multprec_Complex_Numbers;
 --   use Multprec_Complex_Matrices;

 --   n : constant natural := v'length;
 --   res : Matrix(1..n,1..n);

 -- begin
 --   for i in res'range(1) loop
 --     for j in res'range(2) loop
 --       res(i,j) := v(i)**(j-1);
 --     end loop;
 --   end loop;
 --   return res;
 -- end Vdm_Matrix;

 -- procedure lufac_Solve ( n : in natural;
 --                         mat : in Multprec_Complex_Matrices.Matrix;
 --                         rhs : in Multprec_Complex_Vectors.Vector ) is
 --
 --   use Multprec_Complex_Vectors;
 --   use Multprec_Complex_Vectors_io;
 --   use Multprec_Complex_Matrices;
 --   use Multprec_Complex_Matrices_io;
 --
 --   wrk : Matrix(mat'range(1),mat'range(2)) := +mat;
 --   piv : Standard_Integer_Vectors.Vector(mat'range(2));
 --   info : natural;
 --   res,sol,acc : Vector(rhs'range);
 --   nrm : Floating_Number;
 --
 -- begin
 --   put_line("Solving the linear system with lufac.");
 --   lufac(wrk,n,piv,info);
 --   put("info : "); put(info,1); new_line;
 --   sol := +rhs;
 --   lusolve(wrk,n,piv,sol);
 --   put_line("The solution vector :"); put(sol); new_line;
 --   acc := mat*sol;
 --   res := rhs - acc;
 --   put_line("The residual : "); put(res); new_line;
 --   nrm := Max_Norm(res);
 --   put("Max norm of residual : "); put(nrm); new_line;
 --   Clear(nrm);
 --   nrm := Sum_Norm(res);
 --   put("Sum norm of residual : "); put(nrm); new_line;
 --   Clear(nrm);
 --   Clear(wrk); Clear(acc); Clear(res); Clear(sol);
 -- end lufac_Solve;

  procedure lufco_Solve ( n : in integer32;
                          mat : in Multprec_Complex_Matrices.Matrix;
                          rhs : in Multprec_Complex_Vectors.Vector ) is

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Vectors_io;
    use Multprec_Complex_Matrices;

    wrk : Matrix(mat'range(1),mat'range(2)) := +mat;
    piv : Standard_Integer_Vectors.Vector(mat'range(2));
    rcond,nrm : Floating_Number;
    res,sol,acc : Vector(rhs'range);

  begin
    put_line("Solving the linear system with lufco.");
    lufco(wrk,n,piv,rcond);
    put("inverse condition : "); put(rcond); new_line;
    declare
      anorm : Floating_Number := Norm1(mat);
      wrk1 : Matrix(mat'range(1),mat'range(2)) := +mat;
      ipvt : Standard_Integer_Vectors.Vector(mat'range(2));
      info : integer32;
    begin
      lufac(wrk1,n,ipvt,info);
      if info = 0
       then estco(wrk1,n,ipvt,anorm,rcond);
       else rcond := Create(0.0);
      end if;
      put_line("computed condition number again with estco");
      put("inverse condition : "); put(rcond); new_line;
    end;
    sol := +rhs;
    lusolve(wrk,n,piv,sol);
    put_line("The solution vector :"); put(sol); new_line;
    acc := mat*sol;
    res := rhs - acc;
    put_line("The residual : "); put(res); new_line;
    nrm := Max_Norm(res);
    put("Max norm of residual : "); put(nrm); new_line;
    Clear(nrm);
    nrm := Sum_Norm(res);
    put("Sum norm of residual : "); put(nrm); new_line;
    Clear(wrk); Clear(acc); Clear(res); Clear(sol);
    Clear(nrm); Clear(rcond);
  end lufco_Solve;

  procedure Random_Test_Multprec_Linear_Solvers is

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Vectors_io;
    use Multprec_Complex_Matrices;

    n : integer32 := 0;
    sz,nb : natural32 := 0;
    timer : Timing_Widget;

  begin
    new_line;
    put_line("Testing of solving random multi-precision linear systems.");
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the size of the numbers : "); get(sz);
    put("Give the number of tests : "); get(nb);
    tstart(timer);
    for i in 1..nb loop
      declare
       -- rnd : Vector(1..n) := Random_Vector(1,n,sz);
        mat : Matrix(1..n,1..n)
            := Random_Matrix(natural32(n),natural32(n),sz); --Vdm_Matrix(rnd);
        rhs : Vector(1..n) := Random_Vector(1,n,sz);
      begin
       -- lufac_Solve(n,mat,rhs);
        lufco_Solve(n,mat,rhs);
       -- Clear(rnd);
        Clear(mat);
        Clear(rhs);
      end;
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Random_Test_Multprec_Linear_Solvers;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of LU factorization on complex numbers");
    loop
      new_line;
      put_line("Choose one of the following :                               ");
      put_line("  0. exit this program.                                     ");
      put_line("  1. test on solving random standard linear systems.        ");
      put_line("  2. test on solving random multi-precision complex numbers.");
      put("Type 0, 1, or 2 to make your choice : "); get(ans);
      exit when (ans = '0');
      case ans is
        when '1' => Random_Test_Standard_Linear_Solvers;
        when '2' => Random_Test_Multprec_Linear_Solvers;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_cmplu;
