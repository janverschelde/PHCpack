with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Cyclic_Roots_System;
with Cyclic_Laurent_System;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with DEMiCs_Translated;
with DEMiCs_Output_Cells;

package body Test_DEMiCs_Translated is

  function Eigenvalue_Polynomial
             ( dim,idx : integer32; prb : double_float := 0.5 )
             return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the polynomial with index idx of an eigenvalue problem
  --   of dimension dim, so the polynomial has dim variables.
  --   The first variable plays the role of the eigenvalue.
  --   So, the polynomial starts with x(1)*x(idx+1).
  --   The other coefficients are one with probability prb.
  --   If prb equals one, then the eigenvalue problem is dense,
  --   and the smaller prb is, the sparser the polynomial.
  --   The constant term is added to each polynomial.

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    trm : Term;
    ldg : Standard_Natural_Vectors.Link_to_Vector
        := new Standard_Natural_Vectors.Vector'(1..dim+1 => 0);
    rnd : double_float;

  begin
    trm.cf := Create(1.0);
    trm.dg := Degrees(ldg);
    trm.dg(1) := 1;
    trm.dg(idx+1) := 1;
    res := Create(trm);
    trm.dg(1) := 0;
    trm.dg(idx+1) := 0;
    for i in 1..dim loop
      rnd := abs(Standard_Random_Numbers.Random);
      if rnd < prb then
        trm.dg(i+1) := 1;
        Add(res,trm);
        trm.dg(i+1) := 0;
      end if;
    end loop;
    Standard_Natural_Vectors.Clear(ldg);
    return res;
  end Eigenvalue_Polynomial;

  function Linear_Polynomial
             ( dim : integer32 ) return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns a the sum of the dim variables of an eigenvector,
  --   with also the eigenvalue variable added.

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    trm : Term;
    ldg : Standard_Natural_Vectors.Link_to_Vector
        := new Standard_Natural_Vectors.Vector'(1..dim+1 => 0);

  begin
    trm.cf := Create(1.0);
    trm.dg := Degrees(ldg);
    res := Create(trm);
    for i in 1..dim+1 loop
      trm.dg(i) := 1;
      Add(res,trm);
      trm.dg(i) := 0;
    end loop;
    Standard_Natural_Vectors.Clear(ldg);
    return res;
  end Linear_Polynomial;

  function Eigenvalue_System 
             ( dim : integer32; prb : double_float := 0.5 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns a polynomial system which represents an eigenvalue
  --   problem of dimension dim, with sparsity controlled by prb.
  --   The last polynomial, at position dim+1 is linear.

    use Standard_Complex_Poly_Systems;

    res : Poly_Sys(1..dim+1);

  begin
    for i in 1..dim loop
      res(i) := Eigenvalue_Polynomial(dim,i,prb);
    end loop;
    res(dim+1) := Linear_Polynomial(dim);
    return res;
  end Eigenvalue_System;

  procedure Test_Eigenvalue_Problem
              ( dim : in integer32; vrblvl : in integer32 := 0 ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..dim+1)
      := Eigenvalue_System(dim);
    mv : integer32;

  begin
    new_line;
    put("an eigenvalue problem of dimension ");
    put(dim,1); put_line(" : "); put(p);
    mv := DEMiCs_Translated.Mixed_Volume(p,0,false,false,vrblvl);
    put("mixed volume of a "); put(dim,1);
    put("-dimensional eigenproblem : "); put(mv,1); new_line;
  end Test_Eigenvalue_Problem;

  procedure Test_Cyclic ( dim : in integer32; vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Computes the mixed volume of the cyclic n-roots system
  --   where n equals dim.

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
      := Cyclic_Roots_System.Double_Cyclic_System(dim);
    mv : integer32;

  begin
    if vrblvl > 0 then
      new_line;
      put("the cyclic "); put(dim,1); put_line("-roots polynomials : ");
      put(p);
    end if;
    mv := DEMiCs_Translated.Mixed_Volume(p,0,false,false,vrblvl);
    put("mixed volume of cyclic "); put(dim,1);
    put("-roots : "); put(mv,1); new_line;
  end Test_Cyclic;

  procedure Test_Reformulated
              ( dim : in integer32; vrblvl : in integer32 := 0 ) is

    p : constant Standard_Complex_Laur_Systems.Laur_Sys(1..dim-1)
      := Cyclic_Laurent_System.Cyclic_System(natural32(dim));
    mv : integer32;

  begin
    if vrblvl > 0 then
      new_line;
      put("the reformulated cyclic "); put(dim,1);
      put_line("-roots polynomials : "); put(p);
    end if;
    mv := DEMiCs_Translated.Mixed_Volume(p,0,false,vrblvl);
    put("mixed volume of reformulated cyclic "); put(dim,1);
    put("-roots : "); put(mv,1); new_line;
  end Test_Reformulated;

  procedure Test_Labels ( dim : in integer32; vrblvl : in integer32 := 0 ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
      := Cyclic_Roots_System.Double_Cyclic_System(dim);
    mv : integer32;

  begin
    if vrblvl > 0 then
      new_line;
      put("the cyclic "); put(dim,1); put_line("-roots polynomials : ");
      put(p);
    end if;
    mv := DEMiCs_Translated.Mixed_Labels(p,true,0,false,false,vrblvl);
    put("mixed volume of cyclic "); put(dim,1);
    put("-roots : "); put(mv,1); new_line;
  end Test_Labels;

  procedure Test_Cells ( dim : in integer32; vrblvl : in integer32 := 0 ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
      := Cyclic_Roots_System.Double_Cyclic_System(dim);
    mv : integer32;
    mcc : Mixed_Subdivision;

  begin
    if vrblvl > 0 then
      new_line;
      put("the cyclic "); put(dim,1); put_line("-roots polynomials : ");
      put(p);
    end if;
    mv := DEMiCs_Translated.Mixed_Labels(p,true,0,false,false,vrblvl);
    put("mixed volume of cyclic "); put(dim,1);
    put("-roots : "); put(mv,1); new_line;
    mcc := DEMiCs_Translated.Mixed_Cells(vrblvl);
    put("number of mixed cells : ");
    put(integer32(Length_Of(mcc)),1); new_line;
  end Test_Cells;

  procedure Test_Cyclic_Roots ( vrblvl : in integer32 := 0 ) is
  begin
    new_line;
    put_line("-> running tests on the cyclic n-roots system ...");
    for dim in 3..11 loop
      Test_Cyclic(integer32(dim),vrblvl);
    end loop;
  end Test_Cyclic_Roots;

  procedure Test_Reformulated_Cyclic ( vrblvl : in integer32 := 0 ) is
  begin
    new_line;
    put_line("-> running tests on the reformulated cyclic n-roots system ...");
    for dim in 3..11 loop
      Test_Reformulated(integer32(dim),vrblvl);
    end loop;
  end Test_Reformulated_Cyclic;

  procedure Polynomial_User_Lifting ( vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Tests user lifting on polynomial system.

    c3 : constant Standard_Complex_Poly_Systems.Poly_Sys(1..3)
       := Cyclic_Roots_System.Double_Cyclic_System(3);
    mv : integer32;
    mcc : Mixed_Subdivision;
    mix : Standard_Integer_Vectors.Link_to_Vector;

  begin
    mv := DEMiCs_Translated.Mixed_Labels(c3,true,0,false,true,vrblvl);
    put("mixed volume of cyclic 3-roots : "); put(mv,1); new_line;
    mcc := DEMiCs_Translated.Mixed_Cells(vrblvl);
    mix := DEMiCs_Output_Cells.Get_Mixture;
    Floating_Mixed_Subdivisions_io.put(3,mix.all,mcc);
  end Polynomial_User_Lifting;

  procedure Laurent_User_Lifting ( vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Tests user lifting on Laurent polynomial system.

    c3 : constant Standard_Complex_Laur_Systems.Laur_Sys(1..2)
       := Cyclic_Laurent_System.Cyclic_System(3);
    mv : integer32;
    mcc : Mixed_Subdivision;
    mix : Standard_Integer_Vectors.Link_to_Vector;

  begin
    mv := DEMiCs_Translated.Mixed_Labels(c3,true,0,true,vrblvl);
    put("mixed volume of cyclic 3-roots : "); put(mv,1); new_line;
    mcc := DEMiCs_Translated.Mixed_Cells(vrblvl);
    mix := DEMiCs_Output_Cells.Get_Mixture;
    Floating_Mixed_Subdivisions_io.put(3,mix.all,mcc);
  end Laurent_User_Lifting;

  procedure Test_User_Lifting ( vrblvl : in integer32 := 0 ) is

    ans : character;

  begin
    new_line;
    put_line("-> testing user defined lifting ...");
    put("Test on Laurent system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Laurent_User_Lifting(vrblvl);
     else Polynomial_User_Lifting(vrblvl);
    end if;
  end Test_User_Lifting;

  procedure Main is

    vrblvl : integer32 := 99; -- default verbose level
    ans : character;
  
  begin
    put_line("Testing the DEMiCs algorithm ...");
    new_line;
    put("Intermediate output wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'n'
     then vrblvl := 0;
    end if;
    new_line;
    put_line("MENU for testing the translated DEMiCs :");
    put_line("  0. test a sparse eigenvalue problem");
    put_line("  1. run sequence of cyclic n-roots problems");
    put_line("  2. compute labels to points in the mixed cells");
    put_line("  3. convert labels into mixed cells");
    put_line("  4. test on the reformulated cyclic n-roots systems");
    put_line("  5. test user defined lifting");
    put("Type 0, 1, 2, 3, 4, or 5 to select a test : ");
    Ask_Alternative(ans,"012345");
    case ans is
      when '0' => Test_Eigenvalue_Problem(5,vrblvl);
      when '1' => Test_Cyclic_Roots(vrblvl);
      when '2' => Test_Labels(5,vrblvl);
      when '3' => Test_Cells(5,vrblvl);
      when '4' => Test_Reformulated_Cyclic(vrblvl);
      when '5' => Test_User_Lifting(vrblvl);
      when others => null;
    end case;
  end Main;

end Test_DEMiCs_Translated;
