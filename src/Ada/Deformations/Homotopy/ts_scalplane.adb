with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Projective_Transformations;         use Projective_Transformations;
with Homogenization;                     use Homogenization;
with Hyperplane_Solution_Scaling;        use Hyperplane_Solution_Scaling;

procedure ts_scalplane is

-- DESCRIPTION :
--   Test on readjusting the hyperplane in a projective transformation,
--   after scaling a solution.

  procedure Scale ( p : in out Poly_sys; evp : in Eval_Coeff_Poly_Sys;
                    cff : in out Standard_Complex_VecVecs.VecVec;
                    v : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Scales the vector v and evaluates the scaled v at p.

    y : Standard_Complex_Vectors.Vector(p'range);
    z : Standard_Complex_Vectors.Vector(p'range);

  begin
    y := Eval(p,v);
    z := Eval(evp,cff,v);
    put_line("Evaluation at p before coefficient adjustment :"); put_line(y);
    put_line("Evaluation at evp before coefficient adjustment :"); put_line(z);
    Sub(p(p'last),y(y'last));
    y := Eval(p,v);
    put_line("Evaluation before scaling : "); put_line(y);
    put_line("The vector before scaling : "); put_line(v);
    Scale(v);
    put_line("The vector after scaling : "); put_line(v);
    y := Eval(p,v);
    Sub(p(p'last),y(y'last));
    y := Eval(p,v);
    put_line("Evaluation at p after scaling : "); put_line(y);
    Adjust(cff(cff'last),v);
    z := Eval(evp,cff,v);  
    put_line("Evaluation at evp after scaling and coefficient adjustment :");
    put_line(z);
  end Scale;

  procedure Scale ( p : in out Poly_sys; evp : in Eval_Coeff_Poly_Sys;
                    cff : in out Standard_Complex_VecVecs.VecVec;
                    s : in out Solution_List ) is

  -- DESCRIPTION :
  --   Scales the solution vectors in s and checks their evaluation in p.

    tmp : Solution_List := s;
    ls : Link_to_Solution;
    ans : character;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Scale(p,evp,cff,ls.v);
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end Scale;

  procedure Projective_Transformation
              ( p : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Transform the system p and the solutions into projective coordinates.

    pt : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1);

  begin
    Projective_Transformation(sols);
    Projective_Transformation(p.all);
    pt := Add_Random_Hyperplanes(p.all,1,false);
    Standard_Complex_Poly_Systems.Clear(p);
    p := new Standard_Complex_Poly_Systems.Poly_Sys'(pt);
  end Projective_Transformation;

  procedure Test ( p : in out Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Runs a test on scaling the solutions and the last equation in p.

    evp : Eval_Coeff_Poly_Sys(p'range)
        := Standard_Complex_Poly_SysFun.Create(p);
    cff : Standard_Complex_VecVecs.VecVec(p'range)
        := Standard_Complex_Poly_SysFun.Coeff(p);
  begin
    put_line("The last polynomial in p : ");
    put(p(p'last)); new_line;
    put_line("The last coefficient vector :");
    put_line(cff(cff'last).all);
    Scale(p,evp,cff,sols);
    Standard_Complex_Poly_SysFun.Clear(evp);
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a system and its solutions.
  --   Then the test is done.

    p : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading a systems and its solutions ...");
    Standard_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    Projective_Transformation(p,sols);
    Test(p.all,sols);
  end Main;

begin
  Main;
end ts_scalplane;
