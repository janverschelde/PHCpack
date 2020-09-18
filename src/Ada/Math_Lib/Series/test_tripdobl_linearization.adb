with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Vectors_io;        use TripDobl_Complex_Vectors_io;
with TripDobl_Complex_Vector_Norms;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Series_Vectors_io; use TripDobl_Complex_Series_Vectors_io;
with TripDobl_Complex_Series_Matrices;
with TripDobl_Complex_Vector_Series_io;  use TripDobl_Complex_Vector_Series_io;
with TripDobl_Complex_Matrix_Series;
with TripDobl_Complex_Matrix_Series_io;  use TripDobl_Complex_Matrix_Series_io;
with TripDobl_Random_Series_Vectors;
with TripDobl_Random_Series_Matrices;
with TripDobl_Series_Matrix_Solvers;

package body Test_TripDobl_Linearization is

  procedure Write_Difference
              ( x,y : in TripDobl_Complex_Vector_Series.Vector ) is

    dim : constant integer32 := x.cff(0)'last;
    dif : TripDobl_Complex_Vectors.Vector(1..dim);
    nrm,err : triple_double := create(0.0);

    use TripDobl_Complex_Vectors;

  begin
    for k in 0..x.deg loop
      dif := x.cff(k).all - y.cff(k).all;
      nrm := TripDobl_Complex_Vector_Norms.Max_Norm(dif);
      put("Max norm of error at component "); put(k,1);
      put(" : "); put(nrm,3); new_line;
      if nrm > err
       then err := nrm;
      end if;
    end loop;
    put("Max norm of the error : "); put(err,3); new_line;
  end Write_Difference;

  procedure TripDobl_Test ( n,m,d : in integer32 ) is

    use TripDobl_Complex_Series_Matrices;
    use TripDobl_Series_Matrix_Solvers;

    sA : constant TripDobl_Complex_Series_Matrices.Matrix(1..n,1..m)
       := TripDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,d);
    As : constant TripDobl_Complex_Matrix_Series.Matrix 
       := TripDobl_Complex_Matrix_Series.Create(sA); 
    sx : constant TripDobl_Complex_Series_Vectors.Vector(1..m)
       := TripDobl_Random_Series_Vectors.Random_Series_Vector(1,m,d);
    xs : constant TripDobl_Complex_Vector_Series.Vector(d)
       := TripDobl_Complex_Vector_Series.Create(sx);
    sb : constant TripDobl_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant TripDobl_Complex_Vector_Series.Vector(d)
       := TripDobl_Complex_Vector_Series.Create(sb);
    ys : TripDobl_Complex_Vector_Series.Vector(d);
    ans : character;
    rcond : triple_double;
    info : integer32;

  begin
    put_line("The coefficients of the matrix series :"); put(As);
    put_line("The exact solution x :"); put_line(sx);
    put_line("The coefficients of the vector series x :"); put(xs);
    put_line("The right hand side vector b :"); put_line(sb);
    put_line("The coefficients of the vector series b :"); put(bs);
    new_line;
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
  end TripDobl_Test;

  procedure TripDobl_Timing ( n,m,d,f : in integer32 ) is

    sA : constant TripDobl_Complex_Series_Matrices.Matrix(1..n,1..m)
       := TripDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,d);
    As : constant TripDobl_Complex_Matrix_Series.Matrix 
       := TripDobl_Complex_Matrix_Series.Create(sA); 
    sb : constant TripDobl_Complex_Series_Vectors.Vector(1..n)
       := TripDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    bs : constant TripDobl_Complex_Vector_Series.Vector(d)
       := TripDobl_Complex_Vector_Series.Create(sb);
    xs : TripDobl_Complex_Vector_Series.Vector(d);
    info : integer32;
    timer : Timing_Widget;

    use TripDobl_Series_Matrix_Solvers;

  begin
    if n = m then
      tstart(timer);
      for k in 1..f loop
        Solve_by_lufac(As,bs,info,xs);
        TripDobl_Complex_Vector_Series.Clear(xs); -- test on memory leaks
      end loop;
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"Solve by LUfac");
    else
      tstart(timer);
      for k in 1..f loop
        Solve_by_QRLS(As,bs,info,xs);
        TripDobl_Complex_Vector_Series.Clear(xs); -- test on memory leaks
      end loop;
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"Solve by QRLS");
    end if;
  end TripDobl_Timing;

  procedure Main is

    neq,nvr,deg,frq : integer32 := 0;

  begin
    new_line;
    put_line("Testing the linearization of systems of power series ...");
    put("  Give the number of equations in the system : "); get(neq);
    put("  Give the number of variables in the system : "); get(nvr);
    put("  Give the degree of the series : "); get(deg);
    put("  Give frequency of testing (0 for interactive) : "); get(frq);
    new_line;
    if frq = 0
     then TripDobl_Test(neq,nvr,deg);
     else TripDobl_Timing(neq,nvr,deg,frq);
    end if;
  end Main;

end Test_TripDobl_Linearization;
