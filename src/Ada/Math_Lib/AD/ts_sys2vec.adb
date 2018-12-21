with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Standard_Random_Vectors;
with Standard_Complex_Vector_Norms;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Matrices;
with Standard_Complex_Matrix_Norms;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrix_Norms;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrix_Norms;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Standard_Polynomial_Vectors;
with Standard_Polynomial_Vectors_io;    use Standard_Polynomial_Vectors_io;
with DoblDobl_Polynomial_Vectors;
with DoblDobl_Polynomial_Vectors_io;    use DoblDobl_Polynomial_Vectors_io;
with QuadDobl_Polynomial_Vectors;
with QuadDobl_Polynomial_Vectors_io;    use QuadDobl_Polynomial_Vectors_io;
with System_Vector_Convertors;          use System_Vector_Convertors;

procedure ts_sys2vec is

-- DESCRIPTION :
--   Development of the conversion from a polynomial system into
--   the polynomial vectors for efficient evaluation and differentiation.

  procedure Standard_Eval
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                q : in Standard_Polynomial_Vectors.System ) is

  -- DESCRIPTION :
  --   Evaluates p and q at a random point in standard double precision.

    nvr : constant natural32
        := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    dim : constant integer32 := integer32(nvr);
    xpt : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim);
    yvp : Standard_Complex_Vectors.Vector(p'range);
    yjm : Standard_Complex_Matrices.Matrix(p'range,xpt'range);
    zvp : Standard_Complex_Vectors.Vector(q.pols'range);
    zjm : Standard_Complex_Matrices.Matrix(p'range,xpt'range);
    jmt : Standard_Complex_Jaco_Matrices.Jaco_Mat(p'range,xpt'range)
        := Standard_Complex_Jaco_Matrices.Create(p);
    nrm : double_float;

  begin
    yvp := Standard_Complex_Poly_SysFun.Eval(p,xpt);
    put_line("Straightforward evaluation : "); put_line(yvp);
    Standard_Polynomial_Vectors.Speel(q,xpt,zvp,zjm);
    put_line("Speelpenning evaluation : "); put_line(zvp);
    Standard_Complex_Vectors.Sub(yvp,zvp);
    nrm := Standard_Complex_Vector_Norms.Max_Norm(yvp);
    put("Max norm of the evaluation difference :"); put(nrm,3); new_line;
    yjm := Standard_Complex_Jaco_Matrices.Eval(jmt,xpt);
    Standard_Complex_Matrices.Sub(yjm,zjm);
    nrm := Standard_Complex_Matrix_Norms.Max_Norm(yjm);
    put("Max norm of the differentiation difference :"); put(nrm,3); new_line;
  end Standard_Eval;

  procedure DoblDobl_Eval
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                q : in DoblDobl_Polynomial_Vectors.System ) is

  -- DESCRIPTION :
  --   Evaluates p and q at a random point in double double precision.

    nvr : constant natural32
        := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    dim : constant integer32 := integer32(nvr);
    xpt : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    yvp : DoblDobl_Complex_Vectors.Vector(p'range);
    yjm : DoblDobl_Complex_Matrices.Matrix(p'range,xpt'range);
    zvp : DoblDobl_Complex_Vectors.Vector(q.pols'range);
    zjm : DoblDobl_Complex_Matrices.Matrix(p'range,xpt'range);
    jmt : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(p'range,xpt'range)
        := DoblDobl_Complex_Jaco_Matrices.Create(p);
    nrm : double_double;

  begin
    yvp := DoblDobl_Complex_Poly_SysFun.Eval(p,xpt);
    put_line("Straightforward evaluation : "); put_line(yvp);
    DoblDobl_Polynomial_Vectors.Speel(q,xpt,zvp,zjm);
    put_line("Speelpenning evaluation : "); put_line(zvp);
    DoblDobl_Complex_Vectors.Sub(yvp,zvp);
    nrm := DoblDobl_Complex_Vector_Norms.Max_Norm(yvp);
    put("Max norm of the evaluation difference : "); put(nrm,3); new_line;
    yjm := DoblDobl_Complex_Jaco_Matrices.Eval(jmt,xpt);
    DoblDobl_Complex_Matrices.Sub(yjm,zjm);
    nrm := DoblDobl_Complex_Matrix_Norms.Max_Norm(yjm);
    put("Max norm of the differentiation difference : ");
    put(nrm,3); new_line;
  end DoblDobl_Eval;

  procedure QuadDobl_Eval
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                q : in QuadDobl_Polynomial_Vectors.System ) is

  -- DESCRIPTION :
  --   Evaluates p and q at a random point in quad double precision.

    nvr : constant natural32
        := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    dim : constant integer32 := integer32(nvr);
    xpt : constant QuadDobl_Complex_Vectors.Vector(1..dim)
        := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    yvp : QuadDobl_Complex_Vectors.Vector(p'range);
    yjm : QuadDobl_Complex_Matrices.Matrix(p'range,xpt'range);
    zvp : QuadDobl_Complex_Vectors.Vector(q.pols'range);
    zjm : QuadDobl_Complex_Matrices.Matrix(p'range,xpt'range);
    jmt : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(p'range,xpt'range)
        := QuadDobl_Complex_Jaco_Matrices.Create(p);
    nrm : quad_double;

  begin
    yvp := QuadDobl_Complex_Poly_SysFun.Eval(p,xpt);
    put_line("Straightforward evaluation : "); put_line(yvp);
    QuadDobl_Polynomial_Vectors.Speel(q,xpt,zvp,zjm);
    put_line("Speelpenning evaluation : "); put_line(zvp);
    QuadDobl_Complex_Vectors.Sub(yvp,zvp);
    nrm := QuadDobl_Complex_Vector_Norms.Max_Norm(yvp);
    put("Max norm of the evaluation difference : "); put(nrm,3); new_line;
    yjm := QuadDobl_Complex_Jaco_Matrices.Eval(jmt,xpt);
    QuadDobl_Complex_Matrices.Sub(yjm,zjm);
    nrm := QuadDobl_Complex_Matrix_Norms.Max_Norm(yjm);
    put("Max norm of the differentiation difference : ");
    put(nrm,3); new_line;
  end QuadDobl_Eval;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then converts it
  --   into polynomial vectors, in standard double precision.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    pv : Standard_Polynomial_Vectors.Link_to_System;
    ans : character;
 
  begin
    put_line("Reading a polynomial system ..."); get(lp);
    pv := Convert(lp);
    put_line("The converted system :"); put(pv);
    put("Test evaluation ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Eval(lp.all,pv.all);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then converts it
  --   into polynomial vectors, in double double precision.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    pv : DoblDobl_Polynomial_Vectors.Link_to_System;
    ans : character;
 
  begin
    put_line("Reading a polynomial system ..."); get(lp);
    pv := Convert(lp);
    put_line("The converted system :"); put(pv);
    put("Test evaluation ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then DoblDobl_Eval(lp.all,pv.all);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then converts it
  --   into polynomial vectors, in quad double precision.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    pv : QuadDobl_Polynomial_Vectors.Link_to_System;
    ans : character;
 
  begin
    put_line("Reading a polynomial system ..."); get(lp);
    pv := Convert(lp);
    put_line("The converted system :"); put(pv);
    put("Test evaluation ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then QuadDobl_Eval(lp.all,pv.all);
    end if;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a precision and then launches
  --   the corresponding test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sys2vec;
