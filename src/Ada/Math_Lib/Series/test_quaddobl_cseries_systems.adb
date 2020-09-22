with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Symbol_Table;
with QuadDobl_Complex_Series;
with QuadDobl_Series_Linear_Solvers;    use QuadDobl_Series_Linear_Solvers;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_SysFun;
with QuadDobl_CSeries_Jaco_Matrices;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;

package body Test_QuadDobl_CSeries_Systems is

  procedure Read_Series_Vector
              ( v : out QuadDobl_Complex_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 ) is

    vx : integer32 := 1;

  begin
    if idx = 0 then
      Symbol_Table.Enlarge(1);
      vx := dim+1;
    end if;
    new_line;
    put("Reading a vector of "); put(dim,1); put_line(" series ...");
    Complex_Series_and_Polynomials_io.get(v,vx);
    put_line("The vector of series :");
    Complex_Series_and_Polynomials_io.put(v.all);
  end Read_Series_Vector;

  procedure Test_Evaluation
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

    use QuadDobl_CSeries_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : QuadDobl_Complex_Series_Vectors.Link_to_Vector;
    y : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    degree : integer32 := 0;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the degree of series : "); get(degree);
    Complex_Series_and_Polynomials.Set_Degree(x.all,degree);
    put_line("Evaluating the series ...");
    y := QuadDobl_CSeries_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Complex_Series_and_Polynomials_io.put(y);
  end Test_Evaluation;

  procedure Write ( A : QuadDobl_Complex_Series_Matrices.Matrix ) is

    use QuadDobl_Complex_Series;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A["); put(i,1); put(","); put(j,1); put_line("] :");
        if A(i,j) /= null then
          Complex_Series_and_Polynomials_io.put(A(i,j).all);
          new_line;
        end if;
      end loop;
    end loop;
  end Write;

  procedure Test_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

    use QuadDobl_CSeries_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : QuadDobl_Complex_Series_Vectors.Link_to_Vector;
    dx : QuadDobl_Complex_Series_Vectors.Vector(1..n);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jp : constant QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,1..n)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,1..n);
    degree : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    tol : constant double_float := 1.0e-20;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the degree of series : "); get(degree);
    Complex_Series_and_Polynomials.Set_Degree(x.all,degree);
    put_line("Evaluating the series ...");
    px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Complex_Series_and_Polynomials_io.put(px);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_Degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x.all);
    Complex_Series_and_Polynomials.Set_Degree(jm,degree);
    put_line("The Jacobian matrix : ");
    Write(jm);
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line("The update of the Newton step :");
      Complex_Series_and_Polynomials.Filter(dx,tol);
      Complex_Series_and_Polynomials_io.put(dx);
      QuadDobl_Complex_Series_Vectors.Add(x.all,dx);
      put_line("After adding the update to the current series :");
      Complex_Series_and_Polynomials.Filter(x.all,tol);
      Complex_Series_and_Polynomials_io.put(x.all);
      px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x.all);
      put_line("After evaluation in the original system :");
      Complex_Series_and_Polynomials.Filter(px,tol);
      Complex_Series_and_Polynomials_io.put(px);
    end if;
  end Test_Newton_Step;

  procedure Main is

    ls : QuadDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
    ix : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the index of the series variable : "); get(ix);
    new_line;
    Complex_Series_and_Polynomials_io.get(ls,ix);
    new_line;
    put_line("The polynomial system : ");
    Complex_Series_and_Polynomials_io.put(ls.all,ix);
    new_line;
    put_line("MENU to test systems of series polynomials :");
    put_line("  0. test evaluation in a given series;");
    put_line("  1. run one Newton step.");
    put("Type 0 or 1 to select a test : ");
    Ask_Alternative(ans,"01");
    case ans is
      when '0' => Test_Evaluation(ls.all,ix);
      when '1' => Test_Newton_Step(ls.all,ix);
      when others => null;
    end case;
  end Main;

end Test_QuadDobl_CSeries_Systems;
