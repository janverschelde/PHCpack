with text_io;                           use text_io;
with Timing_Package;                    use Timing_Package;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Random_Vectors;
with Standard_Random_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Lexicographical_Supports;
with Standard_Polynomial_Flatteners;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Random_Vectors;
with DoblDobl_Random_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Polynomial_Flatteners;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Random_Vectors;
with QuadDobl_Random_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Polynomial_Flatteners;
with Multitasking_Polynomial_Functions; use Multitasking_Polynomial_Functions;

procedure ts_mtpolval is

-- DESCRIPTION :
--   Interactive development of multitasking polynomial evaluation.

  procedure Prompt_to_Fix_Seed is

  -- DESCRIPTION :
  --   Prompts the user if the seed of the random number generator
  --   needs to be fixed or not.  Allows the user to give a seed.

    ans : character;
    seed : integer32 := 0;

  begin
    put("Set the seed for the random number generator ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give a value for the seed : "); get(seed);
      put("Setting the seed to "); put(seed,1); put_line(".");
      Standard_Random_Numbers.Set_Seed(natural32(seed));
    end if;
  end Prompt_to_Fix_Seed;

  function Standard_Random_System 
             ( n,d,m : integer32)
             return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Generates a random polynomial system in n variables
  --   of degrees at most d and with m terms (0 for dense).
  --   The system on return also contains all partial derivatives.

    use Standard_Complex_Polynomials;
    use Standard_Random_Polynomials;
    use Standard_Complex_Poly_Systems;

    res : Poly_Sys(1..(n+1)*n);

  begin
    Prompt_to_Fix_Seed;
    for i in 1..n loop
      if m = 0
       then res(i) := Random_Dense_Poly(natural32(n),natural32(d),0);
       else res(i) := Random_Sparse_Poly
                        (natural32(n),natural32(d),natural32(m),0);
      end if;
    end loop;
    for i in 1..n loop
      for j in 1..n loop
        res(i*n+j) := Diff(res(i),j);
      end loop;
    end loop;
    return res;
  end Standard_Random_System;

  function DoblDobl_Random_System 
             ( n,d,m : integer32)
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Generates a random polynomial system in n variables
  --   of degrees at most d and with m terms (0 for dense).
  --   The system on return also contains all partial derivatives.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Random_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    res : Poly_Sys(1..(n+1)*n);

  begin
    Prompt_to_Fix_Seed;
    for i in 1..n loop
      if m = 0
       then res(i) := Random_Dense_Poly(natural32(n),natural32(d),0);
       else res(i) := Random_Sparse_Poly
                        (natural32(n),natural32(d),natural32(m),0);
      end if;
    end loop;
    for i in 1..n loop
      for j in 1..n loop
        res(i*n+j) := Diff(res(i),j);
      end loop;
    end loop;
    return res;
  end DoblDobl_Random_System;

  function QuadDobl_Random_System 
             ( n,d,m : integer32)
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Generates a random polynomial system in n variables
  --   of degrees at most d and with m terms (0 for dense).
  --   The system on return also contains all partial derivatives.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Random_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    res : Poly_Sys(1..(n+1)*n);

  begin
    Prompt_to_Fix_Seed;
    for i in 1..n loop
      if m = 0
       then res(i) := Random_Dense_Poly(natural32(n),natural32(d),0);
       else res(i) := Random_Sparse_Poly
                        (natural32(n),natural32(d),natural32(m),0);
      end if;
    end loop;
    for i in 1..n loop
      for j in 1..n loop
        res(i*n+j) := Diff(res(i),j);
      end loop;
    end loop;
    return res;
  end QuadDobl_Random_System;

  procedure Standard_Dense_Test
              ( nt,n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests evaluation at a random system, assumed to be dense.

    v : Standard_Integer_VecVecs.Link_to_VecVec;
    c : Standard_Complex_Matrices.Link_to_Matrix;
    x : Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    y : Standard_Complex_Vectors.Vector(p'range)
      := Standard_Complex_Poly_SysFun.Eval(p,x);
    z : Standard_Complex_Vectors.Vector(p'range);

  begin
    Standard_Polynomial_Flatteners.Flatten(p,v,c);
    z := Standard_Polynomial_Flatteners.Eval(c.all,v.all,x);
    put_line("y : "); put_line(y);
    put_line("z : "); put_line(z);
    put("Number of monomials : "); put(v'last,1); new_line;
    declare
      vy : Standard_Complex_Vectors.Vector(v'range)
         := Standard_Polynomial_Flatteners.Eval(v.all,x);
      vz : Standard_Complex_Vectors.Vector(v'range);
      z2 : Standard_Complex_Vectors.Vector(p'range);
      ans : character;
    begin
      put("Output during multitasking ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then vz := Reporting_Eval(nt,v.all,x);
            z2 := Reporting_Eval(nt,c.all,v.all,x);
       else vz := Silent_Eval(nt,v.all,x);
            z2 := Silent_Eval(nt,c.all,v.all,x);
      end if;
      put_line("vy : "); put_line(vy);
      put_line("vz : "); put_line(vz);
      put_line("z1 : "); put_line(z);
      put_line("z2 : "); put_line(z2);
    end;
  end Standard_Dense_Test;

  procedure Standard_Sparse_Test
              ( nt,n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests evaluation at a random system, assumed to be sparse.

    sup : List := Standard_Polynomial_Flatteners.Distinct_Supports(p);
    lsp : List := Lexicographical_Supports.Sort(sup);
    v : constant Standard_Integer_VecVecs.VecVec 
      := Lists_of_Integer_Vectors.Shallow_Create(lsp);
    c : Standard_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    x : Standard_Complex_Vectors.Vector(1..n);
    y,z1,z2,z3 : Standard_Complex_Vectors.Vector(p'range);
    ans : character;

  begin
    Standard_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      y := Standard_Complex_Poly_SysFun.Eval(p,x);
      z1 := Standard_Polynomial_Flatteners.Eval(c,v,k,x);
      put("Monitor progress of multitasked evaluation ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then z2 := Reporting_Eval(nt,c,v,k,x);
            z3 := Reporting_Looping_Eval(nt,c,v,k,x);
       else z2 := Silent_Eval(nt,c,v,k,x);
            z3 := Silent_Looping_Eval(nt,c,v,k,x);
      end if;
      put_line("y :"); put_line(y);
      put_line("z1 :"); put_line(z1);
      put_line("z2 :"); put_line(z2);
      put_line("z3 :"); put_line(z3);
      put("More samples ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Sparse_Test;

  procedure Standard_Test ( nt,n,d,m : in integer32 ) is

  -- DESCRIPTION :
  --   Tests evaluation at a random system.

    p : constant Standard_Complex_Poly_Systems.Poly_Sys
      := Standard_Random_System(n,d,m);
    ans : character;

  begin
    put_line("The system : "); put_line(p);
    put("Is the system dense ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Dense_Test(nt,n,p);
     else Standard_Sparse_Test(nt,n,p);
    end if;
  end Standard_Test;

  procedure Standard_Performance_Test ( nt,n,d,m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of dimension n, with largest degree d
  --   and at most m monomials in every equation.
  --   After prompting the user for the number of evaluations,
  --   nt tasks will evaluate the system at randomly generated points.

    p : constant Standard_Complex_Poly_Systems.Poly_Sys
      := Standard_Random_System(n,d,m);
    ntimes : integer32;
    timer : Timing_Widget;
    sup : List := Standard_Polynomial_Flatteners.Distinct_Supports(p);
    lsp : List := Lexicographical_Supports.Sort(sup);
    v : constant Standard_Integer_VecVecs.VecVec 
      := Lists_of_Integer_Vectors.Shallow_Create(lsp);
    c : Standard_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    x : Standard_Complex_Vectors.Vector(1..n);
    y : Standard_Complex_Vectors.Vector(p'range);

  begin
    put("Number of distinct monomials : "); put(v'last,1); new_line;
    put("Give the number of evaluations : "); get(ntimes);
    put_line("... flattening the system ...");
    Standard_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    put_line("... starting the evaluation ...");
    tstart(timer);
    for i in 1..ntimes loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      y := Silent_Eval(nt,c,v,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"multitasked polynomial evaluation");
  end Standard_Performance_Test;

  procedure DoblDobl_Dense_Test
              ( nt,n : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests evaluation at a random system, assuming it is dense.

    v : Standard_Integer_VecVecs.Link_to_VecVec;
    c : DoblDobl_Complex_Matrices.Link_to_Matrix;
    x : DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    y : DoblDobl_Complex_Vectors.Vector(p'range)
      := DoblDobl_Complex_Poly_SysFun.Eval(p,x);
    z : DoblDobl_Complex_Vectors.Vector(p'range);

  begin
    put_line("The system : "); put_line(p);
    DoblDobl_Polynomial_Flatteners.Flatten(p,v,c);
    z := DoblDobl_Polynomial_Flatteners.Eval(c.all,v.all,x);
    put_line("y : "); put_line(y);
    put_line("z : "); put_line(z);
    put("Number of monomials : "); put(v'last,1); new_line;
    declare
      vy : DoblDobl_Complex_Vectors.Vector(v'range)
         := DoblDobl_Polynomial_Flatteners.Eval(v.all,x);
      vz : DoblDobl_Complex_Vectors.Vector(v'range);
      z2 : DoblDobl_Complex_Vectors.Vector(p'range);
      ans : character;
    begin
      put("Output during multitasking ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then vz := Reporting_Eval(nt,v.all,x);
            z2 := Reporting_Eval(nt,c.all,v.all,x);
       else vz := Silent_Eval(nt,v.all,x);
            z2 := Silent_Eval(nt,c.all,v.all,x);
      end if;
      put_line("vy : "); put_line(vy);
      put_line("vz : "); put_line(vz);
      put_line("z1 : "); put_line(z);
      put_line("z2 : "); put_line(z2);
    end;
  end DoblDobl_Dense_Test;

  procedure DoblDobl_Sparse_Test
              ( nt,n : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests evaluation at a random system, assumed to be sparse.

    sup : List := DoblDobl_Polynomial_Flatteners.Distinct_Supports(p);
    lsp : List := Lexicographical_Supports.Sort(sup);
    v : constant Standard_Integer_VecVecs.VecVec 
      := Lists_of_Integer_Vectors.Shallow_Create(lsp);
    c : DoblDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    x : DoblDobl_Complex_Vectors.Vector(1..n);
    y,z1,z2,z3 : DoblDobl_Complex_Vectors.Vector(p'range);
    ans : character;

  begin
    DoblDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      y := DoblDobl_Complex_Poly_SysFun.Eval(p,x);
      z1 := DoblDobl_Polynomial_Flatteners.Eval(c,v,k,x);
      put("Monitor progress of multitasked evaluation ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then z2 := Reporting_Eval(nt,c,v,k,x);
            z3 := Reporting_Looping_Eval(nt,c,v,k,x);
       else z2 := Silent_Eval(nt,c,v,k,x);
            z3 := Silent_Looping_Eval(nt,c,v,k,x);
      end if;
      put_line("y :"); put_line(y);
      put_line("z1 :"); put_line(z1);
      put_line("z2 :"); put_line(z2);
      put_line("z3 :"); put_line(z3);
      put("More samples ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Sparse_Test;

  procedure DoblDobl_Test ( nt,n,d,m : in integer32 ) is

  -- DESCRIPTION :
  --   Tests evaluation at a random system.

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
      := DoblDobl_Random_System(n,d,m);
    ans : character;

  begin
    put_line("The system : "); put_line(p);
    put("Is the system dense ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then DoblDobl_Dense_Test(nt,n,p);
     else DoblDobl_Sparse_Test(nt,n,p);
    end if;
  end DoblDobl_Test;

  procedure DoblDobl_Performance_Test ( nt,n,d,m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of dimension n, with largest degree d
  --   and at most m monomials in every equation.
  --   After prompting the user for the number of evaluations,
  --   nt tasks will evaluate the system at randomly generated points.

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
      := DoblDobl_Random_System(n,d,m);
    ntimes : integer32;
    timer : Timing_Widget;
    sup : List := DoblDobl_Polynomial_Flatteners.Distinct_Supports(p);
    lsp : List := Lexicographical_Supports.Sort(sup);
    v : constant Standard_Integer_VecVecs.VecVec 
      := Lists_of_Integer_Vectors.Shallow_Create(lsp);
    c : DoblDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    x : DoblDobl_Complex_Vectors.Vector(1..n);
    y : DoblDobl_Complex_Vectors.Vector(p'range);

  begin
    put("Number of distinct monomials : "); put(v'last,1); new_line;
    put("Give the number of evaluations : "); get(ntimes);
    put_line("... flattening the system ...");
    DoblDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    put_line("... starting the evaluation ...");
    tstart(timer);
    for i in 1..ntimes loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      y := Silent_Eval(nt,c,v,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"multitasked polynomial evaluation");
  end DoblDobl_Performance_Test;

  procedure QuadDobl_Dense_Test
              ( nt,n : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests evaluation at a random system, assuming it is dense.

    v : Standard_Integer_VecVecs.Link_to_VecVec;
    c : QuadDobl_Complex_Matrices.Link_to_Matrix;
    x : QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    y : QuadDobl_Complex_Vectors.Vector(p'range)
      := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    z : QuadDobl_Complex_Vectors.Vector(p'range);

  begin
    put_line("The system : "); put_line(p);
    QuadDobl_Polynomial_Flatteners.Flatten(p,v,c);
    z := QuadDobl_Polynomial_Flatteners.Eval(c.all,v.all,x);
    put_line("y : "); put_line(y);
    put_line("z : "); put_line(z);
    put("Number of monomials : "); put(v'last,1); new_line;
    declare
      vy : QuadDobl_Complex_Vectors.Vector(v'range)
         := QuadDobl_Polynomial_Flatteners.Eval(v.all,x);
      vz : QuadDobl_Complex_Vectors.Vector(v'range);
      z2 : QuadDobl_Complex_Vectors.Vector(p'range);
      ans : character;
    begin
      put("Output during multitasking ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then vz := Reporting_Eval(nt,v.all,x);
            z2 := Reporting_Eval(nt,c.all,v.all,x);
       else vz := Silent_Eval(nt,v.all,x);
            z2 := Silent_Eval(nt,c.all,v.all,x);
      end if;
      put_line("vy : "); put_line(vy);
      put_line("vz : "); put_line(vz);
      put_line("z1 : "); put_line(z);
      put_line("z2 : "); put_line(z2);
    end;
  end QuadDobl_Dense_Test;

  procedure QuadDobl_Sparse_Test
              ( nt,n : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests evaluation at a random system, assumed to be sparse.

    sup : List := QuadDobl_Polynomial_Flatteners.Distinct_Supports(p);
    lsp : List := Lexicographical_Supports.Sort(sup);
    v : constant Standard_Integer_VecVecs.VecVec 
      := Lists_of_Integer_Vectors.Shallow_Create(lsp);
    c : QuadDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    x : QuadDobl_Complex_Vectors.Vector(1..n);
    y,z1,z2,z3 : QuadDobl_Complex_Vectors.Vector(p'range);
    ans : character;

  begin
    QuadDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    loop
      x := QuadDobl_Random_Vectors.Random_Vector(1,n);
      y := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
      z1 := QuadDobl_Polynomial_Flatteners.Eval(c,v,k,x);
      put("Monitor progress of multitasked evaluation ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then z2 := Reporting_Eval(nt,c,v,k,x);
            z3 := Reporting_Looping_Eval(nt,c,v,k,x);
       else z2 := Silent_Eval(nt,c,v,k,x);
            z3 := Silent_Looping_Eval(nt,c,v,k,x);
      end if;
      put_line("y :"); put_line(y);
      put_line("z1 :"); put_line(z1);
      put_line("z2 :"); put_line(z2);
      put_line("z3 :"); put_line(z3);
      put("More samples ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Sparse_Test;

  procedure QuadDobl_Test ( nt,n,d,m : in integer32 ) is

  -- DESCRIPTION :
  --   Tests evaluation at a random system.

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
      := QuadDobl_Random_System(n,d,m);
    ans : character;

  begin
    put_line("The system : "); put_line(p);
    put("Is the system dense ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then QuadDobl_Dense_Test(nt,n,p);
     else QuadDobl_Sparse_Test(nt,n,p);
    end if;
  end QuadDobl_Test;

  procedure QuadDobl_Performance_Test ( nt,n,d,m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of dimension n, with largest degree d
  --   and at most m monomials in every equation.
  --   After prompting the user for the number of evaluations,
  --   nt tasks will evaluate the system at randomly generated points.

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
      := QuadDobl_Random_System(n,d,m);
    ntimes : integer32;
    timer : Timing_Widget;
    sup : List := QuadDobl_Polynomial_Flatteners.Distinct_Supports(p);
    lsp : List := Lexicographical_Supports.Sort(sup);
    v : constant Standard_Integer_VecVecs.VecVec 
      := Lists_of_Integer_Vectors.Shallow_Create(lsp);
    c : QuadDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    x : QuadDobl_Complex_Vectors.Vector(1..n);
    y : QuadDobl_Complex_Vectors.Vector(p'range);
    looping : character;

  begin
    put("Number of distinct monomials : "); put(v'last,1); new_line;
    put("Give the number of evaluations : "); get(ntimes);
    put("Version with looping workers ? (y/n) "); Ask_Yes_or_No(looping);
    put_line("... flattening the system ...");
    QuadDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    put_line("... starting the evaluation ...");
    if looping = 'y' then
      tstart(timer);
      for i in 1..ntimes loop
        x := QuadDobl_Random_Vectors.Random_Vector(1,n);
        y := Silent_Looping_Eval(nt,c,v,k,x);
      end loop;
      tstop(timer);
    else
      tstart(timer);
      for i in 1..ntimes loop
        x := QuadDobl_Random_Vectors.Random_Vector(1,n);
        y := Silent_Eval(nt,c,v,k,x);
      end loop;
      tstop(timer);
    end if;
    new_line;
    print_times(standard_output,timer,"multitasked polynomial evaluation");
  end QuadDobl_Performance_Test;

  procedure Main is

    n,d,m,nt : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing multitasked polynomial system evaluation");
    new_line;
    put_line("Reading parameters for random polynomial system ...");
    put("Give dimension : "); get(n);
    put("Give largest degree : "); get(d);
    put("Give number of terms (0 for dense) : "); get(m);
    new_line;
    put_line("MENU for type of arithmetic : ");
    put_line("  1. evaluate polynomials in standard complex arithmetic;");
    put_line("  2. evaluate polynomials with complex double doubles;");
    put_line("  3. evaluate polynomials with complex quad doubles;");
    put_line("  4. performance test with standard complex arithmetic;");
    put_line("  5. performance test with double double complex arithmetic;");
    put_line("  6. performance test with quad double complex arithmetic.");
    put("Type 1, 2, 3, 4, 5, or 6 to choose : ");
    Ask_Alternative(ans,"123456");
    new_line;
    put("Give number of tasks : "); get(nt);
    put_line("... generating a random system ...");
    case ans is 
      when '1' => Standard_Test(nt,n,d,m);
      when '2' => DoblDobl_Test(nt,n,d,m);
      when '3' => QuadDobl_Test(nt,n,d,m);
      when '4' => Standard_Performance_Test(nt,n,d,m);
      when '5' => DoblDobl_Performance_Test(nt,n,d,m);
      when '6' => QuadDobl_Performance_Test(nt,n,d,m);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtpolval;
