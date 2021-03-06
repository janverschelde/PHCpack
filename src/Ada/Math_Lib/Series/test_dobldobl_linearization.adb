with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with DoblDobl_Complex_VecMats_io;        use DoblDobl_Complex_VecMats_io;
with DoblDobl_Complex_Singular_Values;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors_io; use DoblDobl_Complex_Series_Vectors_io;
with DoblDobl_Complex_Series_Matrices;
with DoblDobl_Complex_Vector_Series_io;  use DoblDobl_Complex_Vector_Series_io;
with DoblDobl_Complex_Matrix_Series;
with DoblDobl_Complex_Matrix_Series_io;  use DoblDobl_Complex_Matrix_Series_io;
with DoblDobl_Random_Series_Vectors;
with DoblDobl_Random_Series_Matrices;
with DoblDobl_Series_Matrix_Solvers;
with Series_Coefficient_Vectors;

package body Test_DoblDobl_Linearization is

  procedure Write_Difference
              ( x,y : in DoblDobl_Complex_Vector_Series.Vector ) is

    dim : constant integer32 := x.cff(0)'last;
    dif : DoblDobl_Complex_Vectors.Vector(1..dim);
    nrm,err : double_double := create(0.0);

    use DoblDobl_Complex_Vectors;

  begin
    for k in 0..x.deg loop
      dif := x.cff(k).all - y.cff(k).all;
      nrm := DoblDobl_Complex_Vector_Norms.Max_Norm(dif);
      put("Max norm of error at component "); put(k,1);
      put(" : "); put(nrm,3); new_line;
      if nrm > err
       then err := nrm;
      end if;
    end loop;
    put("Max norm of the error : "); put(err,3); new_line;
  end Write_Difference;

  procedure DoblDobl_Test ( n,m,d : in integer32 ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Series_Matrices;
    use DoblDobl_Series_Matrix_Solvers;

    sA : constant DoblDobl_Complex_Series_Matrices.Matrix(1..n,1..m)
       := DoblDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,d);
    As : constant DoblDobl_Complex_Matrix_Series.Matrix 
       := DoblDobl_Complex_Matrix_Series.Create(sA); 
    sx : constant DoblDobl_Complex_Series_Vectors.Vector(1..m)
       := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,m,d);
    xs : constant DoblDobl_Complex_Vector_Series.Vector(d)
       := DoblDobl_Complex_Vector_Series.Create(sx);
    sb : constant DoblDobl_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant DoblDobl_Complex_Vector_Series.Vector(d)
       := DoblDobl_Complex_Vector_Series.Create(sb);
    ys,ysn : DoblDobl_Complex_Vector_Series.Vector(d);
    ans : character;
    rcond : double_double;
    det : Complex_Number;
    info : integer32;

  begin
    put_line("The coefficients of the matrix series :"); put(As);
    put_line("The exact solution x :"); put_line(sx);
    put_line("The coefficients of the vector series x :"); put(xs);
    put_line("The right hand side vector b :"); put_line(sb);
    put_line("The coefficients of the vector series b :"); put(bs);
    new_line;
    put("Solve with Echelon form ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Echelon_Solve(As,bs,det,ys,ysn);
      put("n.deg : "); put(ysn.deg,1); 
      put("  det : "); put(det); new_line;
    else
      if n > m then
        new_line;
        put("Solve with SVD ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'y' then
          Solve_by_SVD(As,bs,info,rcond,ys);
          put("rcond : "); put(rcond,3); new_line;
        else
          Solve_by_QRLS(As,bs,info,ys);
          put("info : "); put(info,1); new_line;
        end if;
      else
        new_line;
        put("Condition number wanted ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'y' then
          Solve_by_lufco(As,bs,rcond,ys);
          put("rcond : "); put(rcond,3); new_line;
        else
          Solve_by_lufac(As,bs,info,ys);
          put("info : "); put(info,1); new_line;
        end if;
      end if;
    end if;
    put_line("The generated leading vector series of the solution :");
    put_line(xs.cff(0));
    put_line("The computed leading vector series of the solution :");
    put_line(ys.cff(0));
    for k in 1..bs.deg loop
      put("The generated term "); put(k,1);
      put_line(" of the vector series of the solution :");
      put_line(xs.cff(k));
      put("The computed term "); put(k,1);
      put_line(" of the vector series of the solution :");
      put_line(ys.cff(k));
    end loop;
    Write_Difference(xs,ys);
  end DoblDobl_Test;

  procedure DoblDobl_Timing ( n,m,d,f : in integer32 ) is

    sA : constant DoblDobl_Complex_Series_Matrices.Matrix(1..n,1..m)
       := DoblDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,d);
    As : constant DoblDobl_Complex_Matrix_Series.Matrix 
       := DoblDobl_Complex_Matrix_Series.Create(sA); 
    sb : constant DoblDobl_Complex_Series_Vectors.Vector(1..n)
       := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    bs : constant DoblDobl_Complex_Vector_Series.Vector(d)
       := DoblDobl_Complex_Vector_Series.Create(sb);
    xs : DoblDobl_Complex_Vector_Series.Vector(d);
    info : integer32;
    timer : Timing_Widget;

    use DoblDobl_Series_Matrix_Solvers;

  begin
    if n = m then
      tstart(timer);
      for k in 1..f loop
        Solve_by_lufac(As,bs,info,xs);
        DoblDobl_Complex_Vector_Series.Clear(xs); -- test on memory leaks
      end loop;
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"Solve by LUfac");
    else
      tstart(timer);
      for k in 1..f loop
        Solve_by_QRLS(As,bs,info,xs);
        DoblDobl_Complex_Vector_Series.Clear(xs); -- test on memory leaks
      end loop;
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"Solve by QRLS");
    end if;
  end DoblDobl_Timing;

  procedure DoblDobl_Coefficient_Test ( n,m,d : in integer32 ) is

    use DoblDobl_Complex_Series_Matrices;
    use DoblDobl_Series_Matrix_Solvers;

    sA : constant DoblDobl_Complex_Series_Matrices.Matrix(1..n,1..m)
       := DoblDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,d);
    As : constant DoblDobl_Complex_Matrix_Series.Matrix 
       := DoblDobl_Complex_Matrix_Series.Create(sA); 
    vm : constant DoblDobl_Complex_VecMats.VecMat(0..As.deg)
       := Series_Coefficient_Vectors.DoblDobl_Series_Coefficients(As);
    sx : constant DoblDobl_Complex_Series_Vectors.Vector(1..m)
       := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,m,d);
    xs : constant DoblDobl_Complex_Vector_Series.Vector(d)
       := DoblDobl_Complex_Vector_Series.Create(sx);
    sb : constant DoblDobl_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant DoblDobl_Complex_Vector_Series.Vector(d)
       := DoblDobl_Complex_Vector_Series.Create(sb);
    sbcff : constant DoblDobl_Complex_VecVecs.VecVec(1..n)
          := Series_Coefficient_Vectors.DoblDobl_Series_Coefficients(sb);
    bscff : constant DoblDobl_Complex_VecVecs.VecVec(0..bs.deg)
          := Series_Coefficient_Vectors.DoblDobl_Series_Coefficients(bs);
    ys : DoblDobl_Complex_Vector_Series.Vector(d);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    qraux,w1,w2,w3,w4,w5 : DoblDobl_Complex_Vectors.Vector(1..n);
    wrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..n);
    ewrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
         := new DoblDobl_Complex_Vectors.Vector(1..m);
    xsol : DoblDobl_Complex_VecVecs.VecVec(0..d);
    ans : character;
    rcond : double_double;
    info : integer32;
    mm : constant integer32
       := DoblDobl_Complex_Singular_Values.Min0(n,m);
    S : DoblDobl_Complex_Vectors.Vector(1..mm);
    U : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    V : DoblDobl_Complex_Matrices.Matrix(1..m,1..m);
    timer1,timer2 : Timing_Widget;

  begin
    put_line("The coefficients of the matrix series :"); put(As);
    put_line("The coefficient matrices : "); put(vm);
    put_line("The exact solution x :"); put_line(sx);
    put_line("The coefficients of the vector series x :"); put(xs);
    put_line("The right hand side vector b :"); put_line(sb);
    put_line("The coefficients of b : "); put_line(sbcff);
    put_line("The coefficients of the vector series b :"); put(bs);
    put_line("The coefficients of the vector series b :"); put_line(bscff);
    if n > m then
      put("Solve with SVD ? (y/n) "); Ask_Yes_or_No(ans);
      for i in xsol'range loop
        xsol(i) := new DoblDobl_Complex_Vectors.Vector(1..m);
      end loop;
      if ans = 'y' then
        tstart(timer1);
        Solve_by_SVD(As,bs,info,rcond,ys);
        tstop(timer1);
        put("rcond : "); put(rcond,3); new_line;
        tstart(timer2);
        Solve_by_SVD(vm,bscff,xsol,S,U,V,info,rcond,ewrk,wrk);
        tstop(timer2);
        put("rcond : "); put(rcond,3); new_line;
      else
        tstart(timer1);
        Solve_by_QRLS(As,bs,info,ys);
        tstop(timer1);
        put("info : "); put(info,1); new_line;
        tstart(timer2);
        Solve_by_QRLS(vm,bscff,xsol,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
        tstop(timer2);
        put("info : "); put(info,1); new_line;
      end if;
    else
      put("Condition number wanted ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        tstart(timer1);
        Solve_by_lufco(As,bs,rcond,ys);
        tstop(timer1);
        put("rcond : "); put(rcond,3); new_line;
        tstart(timer2);
        Solve_by_lufco(vm,bscff,ipvt,rcond,wrk);
        tstop(timer2);
        put("rcond : "); put(rcond,3); new_line;
      else
        tstart(timer1);
        Solve_by_lufac(As,bs,info,ys);
        tstop(timer1);
        put("info : "); put(info,1); new_line;
        tstart(timer2);
        Solve_by_lufac(vm,bscff,ipvt,info,wrk);
        tstop(timer2);
        put("info : "); put(info,1); new_line;
      end if;
    end if;
    put_line("The generated leading vector series of the solution :");
    put_line(xs.cff(0));
    put_line("The computed leading vector series of the solution :");
    put_line(ys.cff(0));
    if n > m then
      put_line("The computed leading vector series of the solution :");
      put_line(xsol(0));
    else
      put_line("The computed leading vector series of the solution :");
      put_line(bscff(0));
    end if;
    for k in 1..bs.deg loop
      put("The generated term "); put(k,1);
      put_line(" of the vector series of the solution :"); put_line(xs.cff(k));
      put("The computed term "); put(k,1);
      put_line(" of the vector series of the solution :"); put_line(ys.cff(k));
      if n > m then
        put("The computed term "); put(k,1);
        put_line(" of the vector series of the solution :"); put_line(xsol(k));
      else
        put("The computed term "); put(k,1);
        put_line(" of the vector series of the solution :"); put_line(bscff(k));
      end if;
    end loop;
    Test_DoblDobl_Linearization.Write_Difference(xs,ys);
    new_line;
    print_times(standard_output,timer1,"first series solver");
    new_line;
    print_times(standard_output,timer2,"second series solver");
  end DoblDobl_Coefficient_Test;

  procedure Main is

    neq,nvr,deg,frq : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing the linearization of systems of power series ...");
    put("  Give the number of equations in the system : "); get(neq);
    put("  Give the number of variables in the system : "); get(nvr);
    put("  Give the degree of the series : "); get(deg);
    put("Test linearization on coefficient vectors ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      DoblDobl_Coefficient_Test(neq,nvr,deg);
    else
      put("  Give frequency of testing (0 for interactive) : "); get(frq);
      new_line;
      if frq = 0
       then DoblDobl_Test(neq,nvr,deg);
       else DoblDobl_Timing(neq,nvr,deg,frq);
      end if;
    end if;
  end Main;

end Test_DoblDobl_Linearization;
