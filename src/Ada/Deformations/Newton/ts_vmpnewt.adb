with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Random_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vector_Tools;
with Multprec_Random_Vectors;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Strings;      use Multprec_Complex_Poly_Strings;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_System_and_Solutions_io;
with Varbprec_Complex_Newton_Steps;      use Varbprec_Complex_Newton_Steps;
with Random_Conditioned_Root_Problems;   use Random_Conditioned_Root_Problems;
with Verification_of_Solutions;          use Verification_of_Solutions;

procedure ts_vmpnewt is

-- DESCRIPTION :
--   Development of variable precision Newton's method.

  function Standard_Initial_Approximation
              ( n : integer32 ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user for n coordinates for an initial approximation
  --   to run Newton's method or generates a random vector.

    res : Standard_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    new_line;
    put("Use random vector in Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      res := Standard_Random_Vectors.Random_Vector(1,n);
    else
      put("Reading "); put(n,1); put(" complex numbers ...");
      for i in 1..n loop
        put("x("); put(i,1); put(") : "); get(res(i));
      end loop;
    end if;
    return res;
  end Standard_Initial_Approximation;

  function DoblDobl_Initial_Approximation
              ( n : integer32 ) return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user for n coordinates for an initial approximation
  --   to run Newton's method or generates a random vector.

    res : DoblDobl_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    new_line;
    put("Use random vector in Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      res := DoblDobl_Random_Vectors.Random_Vector(1,n);
    else
      put("Reading "); put(n,1); put(" complex numbers ...");
      for i in 1..n loop
        put("x("); put(i,1); put(") : "); get(res(i));
      end loop;
    end if;
    return res;
  end DoblDobl_Initial_Approximation;

  function QuadDobl_Initial_Approximation
              ( n : integer32 ) return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user for n coordinates for an initial approximation
  --   to run Newton's method or generates a random vector.

    res : QuadDobl_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    new_line;
    put("Use random vector in Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      res := QuadDobl_Random_Vectors.Random_Vector(1,n);
    else
      put("Reading "); put(n,1); put(" complex numbers ...");
      for i in 1..n loop
        put("x("); put(i,1); put(") : "); get(res(i));
      end loop;
    end if;
    return res;
  end QuadDobl_Initial_Approximation;

  function Multprec_Initial_Approximation
              ( n : integer32; size : natural32 )
              return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user for n coordinates for an initial approximation
  --   to run Newton's method or generates a random vector.

    res : Multprec_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    new_line;
    put("Use random vector in Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      res := Multprec_Random_Vectors.Random_Vector(1,n,size);
    else
      put("Reading "); put(n,1); put(" complex numbers ...");
      for i in 1..n loop
        put("x("); put(i,1); put(") : "); get(res(i));
      end loop;
    end if;
    Multprec_Complex_Vector_Tools.Set_Size(res,size);
    return res;
  end Multprec_Initial_Approximation;

  procedure Standard_Test
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps on p,
  --   in standard double precision, prompting the user for
  --   an initial approximation for a root, 
  --   or generating just a random vector.

    z : Standard_Complex_Vectors.Vector(p'range)
      := Standard_Initial_Approximation(p'last);

  begin
    Standard_Test(p,z);
  end Standard_Test;

  procedure DoblDobl_Test
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in double double precision, prompting the user for an initial
  --   approximation for a solution, or generating a random vector.

    z : DoblDobl_Complex_Vectors.Vector(p'range)
      := DoblDobl_Initial_Approximation(p'last);

  begin
    DoblDobl_Test(p,z);
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in quad double precision, prompting the user for an initial
  --   approximation for a solution or generating a random vector.

    z : QuadDobl_Complex_Vectors.Vector(p'range)
      := QuadDobl_Initial_Approximation(p'last);

  begin
    QuadDobl_Test(p,z);
  end QuadDobl_Test;

  procedure Multprec_Test
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                size : in natural32 ) is

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in arbitrary multiprecision with numbers of the given size,
  --   prompting the user to provide an initial approximation for
  --   a solution, or generating a random vector.

    z : Multprec_Complex_Vectors.Vector(p'range)
      := Multprec_Initial_Approximation(p'last,size);

  begin
    Multprec_Test(p,z);
  end Multprec_Test;

  procedure Standard_Test_on_Given_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and runs the variable precision Newton steps
  --   in standard double precision.

    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    Standard_Test(lp.all);
  end Standard_Test_on_Given_System;

  procedure DoblDobl_Test_on_Given_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and runs the variable precision Newton steps
  --   in double double precision.

    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    DoblDobl_Test(lp.all);
  end DoblDobl_Test_on_Given_System;

  procedure QuadDobl_Test_on_Given_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and runs the variable precision Newton steps
  --   in quad double precision.

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    QuadDobl_Test(lp.all);
  end QuadDobl_Test_on_Given_System;

  procedure Multprec_Test_on_Given_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and runs the variable precision Newton steps
  --   in arbitrary multiprecision.

    use Multprec_Complex_Poly_Systems;
    file : file_type;
    n,m : natural32;
    ls : Link_to_Array_of_Strings;
    deci,size : natural32 := 0;

  begin
    new_line;
    put("Give the number of decimal places : "); get(deci); skip_line;
    size := Decimal_to_Size(deci);
    new_line;
    put_line("Reading a polynomial system ...");
    Read_Name_and_Open_File(file);
    get(file,integer(n),integer(m),ls);
    if Symbol_Table.Number < m
     then Symbol_Table.Init(m);
    end if;
    declare
      p : Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(n))
        := Parse(m,size,ls.all);
    begin
      Multprec_Test(p,size);
      Multprec_Complex_Poly_Systems.Clear(p);
    end;
  end Multprec_Test_on_Given_System;

  procedure Test_in_Fixed_Precision is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then either calls the procedure to generate a random
  --   conditioned root problems or calls the procedure that will 
  --   test variable precision Newton's method on a given system.

    precision,random : character;

  begin
    new_line;
    put_line("Testing variable precision Newton steps ...");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision;");
    put_line("  3. arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(precision,"0123");
    new_line;
    put("Generate random conditioned root problem ? (y/n) ");
    Ask_Yes_or_No(random);
    if random = 'y' then
      new_line;
      put_line("Generating a random conditioned root problem ...");
      Random_Conditioned_Root_Problem(precision);
    else
      case precision is 
        when '0' => Standard_Test_on_Given_System;
        when '1' => DoblDobl_Test_on_Given_System;
        when '2' => QuadDobl_Test_on_Given_System;
        when '3' => Multprec_Test_on_Given_System;
        when others => null;
      end case;
    end if;
  end Test_in_Fixed_Precision;

  procedure Interactive_Test_Variable_Precision is

  -- DESCRIPTION :
  --   Generates a conditioned root problem and then estimates
  --   the loss of decimal places.

    f : Link_to_Array_of_Strings;   
    z : Link_to_String;
    loss,want,accu : integer32 := 0;
    precision : natural32;
    err,rco,res : double_float;
    ans : character;

  begin
    Random_Conditioned_Root_Problem(f,z);
   -- put_line("The polynomial in the system :");
   -- for i in f'range loop
   --   put_line(f(i).all);
   -- end loop;
   -- put_line("The initial approximation for a root :");
   -- put_line(z.all);
    loop
      loss := Estimate_Loss_for_Polynomial_System(f.all,z.all,1000);
      put("Estimated loss of decimal places : "); put(loss,1); new_line;
      new_line;
      put("Give the wanted number of accurate decimal places : ");
      get(want);
      precision := natural32(-loss) + natural32(want);
      put("-> the working precision should be ");
      put(precision,1); put_line(" decimal places.");
      do_Newton_Step_on_Polynomial_System(f.all,z,loss,want,err,rco,res);
      accu := integer32(log10(err));
      put("err :"); put(err,3);
      put("  rco :"); put(rco,3);
      put("  res :"); put(res,3); new_line;
      put("-> number of accurate decimal places : ");
      put(abs(accu),1); new_line;
      put("Do another Newton step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Interactive_Test_Variable_Precision;

  procedure Test_Variable_Precision_Newton_Steps is

  -- DESCRIPTION :
  --   Generates a conditioned root problem and then estimates
  --   the loss of decimal places.

    f : Link_to_Array_of_Strings;   
    z : Link_to_String;
    want,loss : integer32 := 0;
    nbr : natural32 := 0;
    err,rco,res : double_float;
    ans : character;
    laurent : boolean;

  begin
    Random_Conditioned_Root_Problem(f,z);
    new_line;
    put("Consider as Laurent polynomial system ? (y/n) ");
    Ask_Yes_or_No(ans);
    laurent := (ans = 'y');
    new_line;
    loop
      put("Give the wanted number of accurate decimal places : "); get(want);
      put("Give maximum number of steps : "); get(nbr);
      if laurent then
        Newton_Steps_on_Laurent_Polynomials
          (standard_output,f.all,z,want,1000,nbr,loss,err,rco,res);
      else
        Newton_Steps_on_Polynomial_System
          (standard_output,f.all,z,want,1000,nbr,loss,err,rco,res);
      end if;
      put("More steps ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Variable_Precision_Newton_Steps;

  procedure Test_Verification_of_Solutions is

  -- DESCRIPTION :
  --   Calls the operations of the package Verification_of_Solutions.

    wanted,maxitr,maxprc : natural32;
    verbose : boolean;
    p : Link_to_Array_of_Strings;
    sols : Multprec_Complex_Solutions.Solution_List;
    nq,nv,ns,len : natural32;
    file : file_type;
    ans : character;
    laurent : boolean;

  begin
    new_line;
    put_line("Reading a system with solutions from file.");
    Multprec_System_and_Solutions_io.get(nq,nv,p,sols);
    new_line;
    len := Multprec_Complex_Solutions.Length_Of(sols);
    put("Read list of "); put(len,1); put_line(" solutions from file.");
    if len > 0 then
      ns := Symbol_Table.Number;
      put("Number of symbols read : "); put(ns,1); new_line;
      if ns > 0 then
        put("Symbols in the table :");
        Symbol_Table_io.Write; new_line;
      end if;
      Menu_to_set_Parameters(wanted,maxitr,maxprc,verbose);
      new_line;
      put("Consider system as Laurent polynomial system ? (y/n) ");
      Ask_Yes_or_No(ans);
      laurent := (ans = 'y');
      if verbose then
        new_line;
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(file);
        Write_Parameters(file,wanted,maxitr,maxprc,true);
        if laurent then
          Verify_Solutions_of_Laurent_Polynomials
            (file,p.all,sols,wanted,maxitr,maxprc);
        else
          Verify_Solutions_of_Polynomial_System
            (file,p.all,sols,wanted,maxitr,maxprc);
        end if;
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,len,nv,sols);
      else 
        new_line;
        Write_Parameters(standard_output,wanted,maxitr,maxprc,verbose);
        if laurent then
          Verify_Solutions_of_Laurent_Polynomials
            (p.all,sols,wanted,maxitr,maxprc);
        else
          Verify_Solutions_of_Polynomial_System
            (p.all,sols,wanted,maxitr,maxprc);
        end if;
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,len,nv,sols);
      end if;
    end if;
  end Test_Verification_of_Solutions;

  procedure Main is

  -- DESCRIPTION :
  --   Offers the user to select precision first, before generating
  --   the problem, or later, after generating the problem.

    ans : character;

  begin
    new_line;
    put_line("MENU to test variable precision Newton's method :");
    put_line("  1. test in a fixed level of precision;");
    put_line("  2. fix precision after estimating condition numbers;");
    put_line("  3. run sequence of Newton steps in variable precision;");
    put_line("  4. verify solutions of a system on file.");
    put("Type 1, 2, 3, or 4 to select test : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Test_in_Fixed_Precision;
      when '2' => Interactive_Test_Variable_Precision;
      when '3' => Test_Variable_Precision_Newton_Steps;
      when '4' => Test_Verification_of_Solutions;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_vmpnewt;
