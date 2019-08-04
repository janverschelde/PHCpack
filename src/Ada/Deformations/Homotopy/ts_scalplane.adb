with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Projective_Transformations;         use Projective_Transformations;
with Homogenization;                     use Homogenization;
with Hyperplane_Solution_Scaling;        use Hyperplane_Solution_Scaling;

procedure ts_scalplane is

-- DESCRIPTION :
--   Test on readjusting the hyperplane in a projective transformation,
--   after scaling a solution.

  procedure Scale ( p : in out Standard_Complex_Poly_Systems.Poly_sys;
                    evp : in Standard_Complex_Poly_SysFun.Eval_Coeff_Poly_Sys;
                    cff : in out Standard_Complex_VecVecs.VecVec;
                    v : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Scales the vector v and evaluates the scaled v at p,
  --   in standard double precision.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;

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

  procedure Scale ( p : in out DoblDobl_Complex_Poly_Systems.Poly_sys;
                    evp : in DoblDobl_Complex_Poly_SysFun.Eval_Coeff_Poly_Sys;
                    cff : in out DoblDobl_Complex_VecVecs.VecVec;
                    v : in out DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Scales the vector v and evaluates the scaled v at p,
  --   in double double precision.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;

    y : DoblDobl_Complex_Vectors.Vector(p'range);
    z : DoblDobl_Complex_Vectors.Vector(p'range);

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

  procedure Scale ( p : in out QuadDobl_Complex_Poly_Systems.Poly_sys;
                    evp : in QuadDobl_Complex_Poly_SysFun.Eval_Coeff_Poly_Sys;
                    cff : in out QuadDobl_Complex_VecVecs.VecVec;
                    v : in out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Scales the vector v and evaluates the scaled v at p,
  --   in quad double precision.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;

    y : QuadDobl_Complex_Vectors.Vector(p'range);
    z : QuadDobl_Complex_Vectors.Vector(p'range);

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

  procedure Scale ( p : in out Standard_Complex_Poly_Systems.Poly_sys;
                    evp : in Standard_Complex_Poly_SysFun.Eval_Coeff_Poly_Sys;
                    cff : in out Standard_Complex_VecVecs.VecVec;
                    s : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Scales the solution vectors in s and checks their evaluation in p,
  --   in standard double precision.

    use Standard_Complex_Solutions;

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

  procedure Scale ( p : in out DoblDobl_Complex_Poly_Systems.Poly_sys;
                    evp : in DoblDobl_Complex_Poly_SysFun.Eval_Coeff_Poly_Sys;
                    cff : in out DoblDobl_Complex_VecVecs.VecVec;
                    s : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Scales the solution vectors in s and checks their evaluation in p,
  --   in double double precision.

    use DoblDobl_Complex_Solutions;

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

  procedure Scale ( p : in out QuadDobl_Complex_Poly_Systems.Poly_sys;
                    evp : in QuadDobl_Complex_Poly_SysFun.Eval_Coeff_Poly_Sys;
                    cff : in out QuadDobl_Complex_VecVecs.VecVec;
                    s : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Scales the solution vectors in s and checks their evaluation in p,
  --   in quad double precision.

    use QuadDobl_Complex_Solutions;

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
  --   Transform the system p and the solutions into projective coordinates,
  --   in standard double precision.

    pt : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1);

  begin
    Projective_Transformation(sols);
    Projective_Transformation(p.all);
    pt := Add_Random_Hyperplanes(p.all,1,false);
    Standard_Complex_Poly_Systems.Clear(p);
    p := new Standard_Complex_Poly_Systems.Poly_Sys'(pt);
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Transform the system p and the solutions into projective coordinates,
  --   in double double precision.

    pt : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1);

  begin
    Projective_Transformation(sols);
    Projective_Transformation(p.all);
    pt := Add_Random_Hyperplanes(p.all,1,false);
    DoblDobl_Complex_Poly_Systems.Clear(p);
    p := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(pt);
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Transform the system p and the solutions into projective coordinates,
  --   in double double precision.

    pt : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last+1);

  begin
    Projective_Transformation(sols);
    Projective_Transformation(p.all);
    pt := Add_Random_Hyperplanes(p.all,1,false);
    QuadDobl_Complex_Poly_Systems.Clear(p);
    p := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(pt);
  end Projective_Transformation;

  procedure Standard_Test
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs a test on scaling the solutions and the last equation in p,
  --   in standard double precision.

    use Standard_Complex_Poly_SysFun;

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
  end Standard_Test;

  procedure DoblDobl_Test
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs a test on scaling the solutions and the last equation in p,
  --   in double double precision.

    use DoblDobl_Complex_Poly_SysFun;

    evp : Eval_Coeff_Poly_Sys(p'range)
        := DoblDobl_Complex_Poly_SysFun.Create(p);
    cff : DoblDobl_Complex_VecVecs.VecVec(p'range)
        := DoblDobl_Complex_Poly_SysFun.Coeff(p);
  begin
    put_line("The last polynomial in p : ");
    put(p(p'last)); new_line;
    put_line("The last coefficient vector :");
    put_line(cff(cff'last).all);
    Scale(p,evp,cff,sols);
    DoblDobl_Complex_Poly_SysFun.Clear(evp);
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Runs a test on scaling the solutions and the last equation in p,
  --   in quad double precision.

    use QuadDobl_Complex_Poly_SysFun;

    evp : Eval_Coeff_Poly_Sys(p'range)
        := QuadDobl_Complex_Poly_SysFun.Create(p);
    cff : QuadDobl_Complex_VecVecs.VecVec(p'range)
        := QuadDobl_Complex_Poly_SysFun.Coeff(p);
  begin
    put_line("The last polynomial in p : ");
    put(p(p'last)); new_line;
    put_line("The last coefficient vector :");
    put_line(cff(cff'last).all);
    Scale(p,evp,cff,sols);
    QuadDobl_Complex_Poly_SysFun.Clear(evp);
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a system and its solutions.
  --   Then the test is done in standard double precision.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading a systems and its solutions ...");
    Standard_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    Projective_Transformation(p,sols);
    Standard_Test(p.all,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a system and its solutions.
  --   Then the test is done in double double precision.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading a systems and its solutions ...");
    DoblDobl_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    Projective_Transformation(p,sols);
    DoblDobl_Test(p.all,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a system and its solutions.
  --   Then the test is done in quad double precision.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading a systems and its solutions ...");
    QuadDobl_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    Projective_Transformation(p,sols);
    QuadDobl_Test(p.all,sols);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision,");
    put_line("  1. double double precision,");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_scalplane;
