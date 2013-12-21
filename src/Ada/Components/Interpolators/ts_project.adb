with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_VecVecs;
with Standard_Central_Projections;       use Standard_Central_Projections;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Central_Projections;       use Multprec_Central_Projections;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;
with Sample_Points;                      use Sample_Points;
with Projection_Operators;               use Projection_Operators;

procedure ts_project is

  procedure Standard_Random_Test_Single_Projector ( n,rd : in integer32 ) is

  -- DESCRIPTION :
  --   Does a test with the central projection using one base point.
  --   Random points are input and the output is tested whether to lie
  --   on the given hyperplane.  Note that this can only be checked
  --   when n = rd!!

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Complex_VecVecs;
    use Standard_Complex_Solutions;

    hyp : constant Vector(0..n) := Random_Vector(0,n);
    hyps : VecVec(1..rd);
    base : constant Vector(1..n) := Random_Vector(1,n);
    point : constant Vector(1..n) := Random_Vector(1,n);
    proj : Vector(1..rd) := Intersect(hyp,base,point,rd);
    eva : Complex_Number := hyp(0) + hyp(proj'range)*proj;
    m : integer32 := 0;
    sol : Solution(n);
    spt : Standard_Sample;
    p : Standard_Projector;

  begin
    for i in 1..rd loop
      hyps(i) := new Vector'(Random_Vector(0,n));
    end loop;
    p := Create(hyps);
    sol.v := base;
    spt := Create(sol,hyps);
    Update(p,spt,hyp);
    put_line("Projected point : "); put_line(proj);
    put("Evaluation : "); put(eva); new_line;
    proj := Project(p,point);
    put_line("Projected point : "); put_line(proj);
    put("Give number of points : "); get(m);
    declare
      points,prjpts : VecVec(1..m);
    begin
      for i in 1..m loop
        points(i) := new Vector'(Random_Vector(1,n));
      end loop;
      prjpts := Intersect(hyp,base,points,rd);
      for i in 1..m loop
        eva := hyp(0) + hyp(prjpts(i)'range)*prjpts(i).all;
        put("Evaluation : "); put(eva); new_line;
      end loop;
    end;
  end Standard_Random_Test_Single_Projector;

  procedure Standard_Random_Test_Multiple_Projector
              ( n,rd,m : in integer32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Complex_VecVecs;
    use Standard_Complex_Solutions;

    hyp,base : VecVec(1..m);
    hyps : VecVec(1..rd);
    point : constant Vector(1..n) := Random_Vector(1,n);
    projpt : Vector(1..rd);
    eva : Complex_Number;
    sol : Solution(n);
    spt : Standard_Sample;
    p : Standard_Projector;

  begin
    for i in 1..m loop
      hyp(i) := new Vector'(Random_Vector(0,n));
      base(i) := new Vector'(Random_Vector(1,n));
    end loop;
    for i in 1..rd loop
      hyps(i) := new Vector'(Random_Vector(0,n));
    end loop;
    p := Create(hyps);
    for i in 1..m loop
      sol.v := base(i).all;
      spt := Create(sol,hyps);
      Update(p,spt,hyp(i).all);
      Shallow_Clear(spt);
    end loop;
    put_line("Intersecting base points.");
    Intersect_Base_Points(hyp,base);
    for i in 1..m-1 loop
      put("Evaluating in hyperplane "); put(i,1); put_line(" :");
      for j in i+1..m loop
        put("  at base point "); put(j,1); put(" : ");
        eva := hyp(i)(0) + hyp(i)(base(j)'range)*base(j).all;
        put(eva,3); new_line;
      end loop;
    end loop;
    put_line("Projection of random point.");
    projpt := Intersect(hyp,base,point,rd);
    put_line("The projection : "); put_line(projpt);
    put_line("Evaluation of projected point :");
    for i in 1..m loop
      put("  in hyperplane "); put(i,1); put(" : ");
      eva := hyp(i)(0) + hyp(i)(projpt'range)*projpt;
      put(eva,3); new_line;
    end loop;
    projpt := Project(p,point);
    put_line("The projection : "); put_line(projpt);
  end Standard_Random_Test_Multiple_Projector;

  procedure Multprec_Random_Test_Single_Projector
              ( n,rd : in integer32; size : in natural32 ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_VecVecs;
    use Multprec_Complex_Solutions;

    hyp : constant Vector(0..n) := Random_Vector(0,n,size);
    hyps : VecVec(1..rd);
    base : constant Vector(1..n) := Random_Vector(1,n,size);
    point : constant Vector(1..n) := Random_Vector(1,n,size);
    proj : Vector(1..rd) := Intersect(hyp,base,point,rd);
    eva : Complex_Number := hyp(0) + hyp(proj'range)*proj;
    m : integer32 := 0;
    sol : Solution(n);
    spt : Multprec_Sample;
    p : Multprec_Projector;

  begin
    for i in 1..rd loop
      hyps(i) := new Vector'(Random_Vector(0,n,size));
    end loop;
    p := Create(hyps);
    put_line("The projection : "); put_line(proj);
    put("Evaluation : "); put(eva); new_line;
    sol.v := base;
    spt := Create(sol,hyps);
    Update(p,spt,hyp);
    Clear(proj);
    proj := Project(p,point);
    put_line("The projection : "); put_line(proj);
    put("Give number of points : "); get(m);
    declare
      points,prjpts : VecVec(1..m);
    begin
      for i in 1..m loop
        points(i) := new Vector'(Random_Vector(1,n,size));
      end loop;
      prjpts := Intersect(hyp,base,points,rd);
      for i in 1..m loop
        eva := hyp(prjpts(i)'range)*prjpts(i).all;
        Add(eva,hyp(0));
        put("Evaluation : "); put(eva); new_line;
        Clear(eva);
      end loop;
      Clear(points); Clear(prjpts);
    end;
  end Multprec_Random_Test_Single_Projector;

  procedure Multprec_Random_Test_Multiple_Projector
              ( n,rd : in integer32; size : in natural32;
                m : in integer32 ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_VecVecs;

    hyp,base : VecVec(1..m);
    point : constant Vector(1..n) := Random_Vector(1,n,size);
    projpt : Vector(1..rd);
    eva : Complex_Number;

  begin
    for i in 1..m loop
      hyp(i) := new Vector'(Random_Vector(0,n,size));
      base(i) := new Vector'(Random_Vector(1,n,size));
    end loop;
    put_line("Intersecting base points.");
    Intersect_Base_Points(hyp,base);
    for i in 1..m-1 loop
      put("Evaluating in hyperplane "); put(i,1); put_line(" :");
      for j in i+1..m loop
        put("  at base point "); put(j,1); put(" : ");
        eva := hyp(i)(base(j)'range)*base(j).all;
        Add(eva,hyp(i)(0));
        put(eva,3); new_line;
        Clear(eva);
      end loop;
    end loop;
    put_line("Projection of random point.");
    projpt := Intersect(hyp,base,point,rd);
    put_line("Evaluation of projected point :");
    for i in 1..m loop
      put("  in hyperplane "); put(i,1); put(" : ");
      eva := hyp(i)(projpt'range)*projpt;
      Add(eva,hyp(i)(0));
      put(eva,3); new_line;
      Clear(eva);
    end loop;
  end Multprec_Random_Test_Multiple_Projector;

  procedure Main is

    ans : character;
    d,n,rd : integer32 := 0;
    size : natural32 := 0;

  begin
    new_line;
    put_line("Random test on a central projector.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. use one base point with standard numbers.");
      put_line("  2. use multiple base points with standard numbers.");
      put_line("  3. use one base point with multi-precision numbers.");
      put_line("  4. use multiple base points with multi-precision numbers.");
      put("Type 0,1,2,3, or 4 to select your choice : ");
      Ask_Alternative(ans,"01234");
      exit when ans = '0';
      new_line;
      put("Give the origin dimension : "); get(d);
      put("Give the target dimension : "); get(rd);
      case ans is
        when '1' => Standard_Random_Test_Single_Projector(d,rd);
        when '2' => put("Give number of base points : "); get(n);
                    Standard_Random_Test_Multiple_Projector(d,rd,n);
        when '3' => put("Give the size of the numbers : "); get(size);
                    Multprec_Random_Test_Single_Projector(d,rd,size);
        when '4' => put("Give the number of base points : "); get(n);
                    put("Give the size of the numbers : "); get(size);
                    Multprec_Random_Test_Multiple_Projector(d,rd,size,n);
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_project;
