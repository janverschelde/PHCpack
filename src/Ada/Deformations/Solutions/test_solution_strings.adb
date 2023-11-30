with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Standard_Solution_Strings;
with DoblDobl_Solution_Strings;
with QuadDobl_Solution_Strings;
with Multprec_Solution_Strings;
with Solution_String_Splitters;          use Solution_String_Splitters;

package body Test_Solution_Strings is

  procedure Standard_Test_Write
              ( ls : in Standard_Complex_Solutions.Link_to_Solution ) is

    use Standard_Complex_Solutions;
    s : constant string := Standard_Solution_Strings.Write(ls.all);
    ind : integer := s'first;
    sol : Solution(ls.n);
    fail : boolean;

  begin
    put_line("The solution written as a string : ");
    put_line(s);
    put_line("The solution vector part : ");
    put_line(Coordinates(s).all);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(Standard_Solution_Strings.Length(ls.all),1);
    new_line;
    Standard_Solution_Strings.Parse(s,ind,natural32(ls.n),sol,fail);
    put_line("The solution parsed from string : ");
    put(sol); new_line;
  end Standard_Test_Write;

  procedure DoblDobl_Test_Write
              ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution ) is

    use DoblDobl_Complex_Solutions;
    s : constant string := DoblDobl_Solution_Strings.Write(ls.all);
    ind : integer := s'first;
    sol : Solution(ls.n);
    fail : boolean;

  begin
    put_line("The solution written as a string : ");
    put_line(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(DoblDobl_Solution_Strings.Length(ls.all),1);
    new_line;
    DoblDobl_Solution_Strings.Parse(s,ind,natural32(ls.n),sol,fail);
    put_line("The solution parsed from string : ");
    put(sol); new_line;
  end DoblDobl_Test_Write;

  procedure QuadDobl_Test_Write
              ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution ) is

    use QuadDobl_Complex_Solutions;
    s : constant string := QuadDobl_Solution_Strings.Write(ls.all);
    ind : integer := s'first;
    sol : Solution(ls.n);
    fail : boolean;

  begin
    put_line("The solution written as a string : ");
    put_line(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(QuadDobl_Solution_Strings.Length(ls.all),1);
    new_line;
    QuadDobl_Solution_Strings.Parse(s,ind,natural32(ls.n),sol,fail);
    put_line("The solution parsed from string : ");
    put(sol); new_line;
  end QuadDobl_Test_Write;

  procedure Multprec_Test_Write
              ( ls : in Multprec_Complex_Solutions.Link_to_Solution ) is

    use Multprec_Complex_Solutions;
    s : constant string := Multprec_Solution_Strings.Write(ls.all);
    ind : integer := s'first;
    sol : Solution(ls.n);
    fail : boolean;

  begin
    put_line("The solution written as a string : ");
    put_line(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(Multprec_Solution_Strings.Length(ls.all),1);
    new_line;
    Multprec_Solution_Strings.Parse(s,ind,natural32(ls.n),sol,fail);
    put_line("The solution parsed from string : ");
    put(sol); new_line;
  end Multprec_Test_Write;

  procedure Standard_Test_Write_Intro
              ( ls : in Standard_Complex_Solutions.Link_to_Solution ) is

    use Standard_Complex_Numbers;
    s : constant string := Standard_Solution_Strings.Write_Intro(ls.all);
    t : Complex_Number;
    m : integer32;
    ind : integer := s'first;
    fail : boolean;

  begin
    put_line("Introduction written as string :");
    put(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(Standard_Solution_Strings.Length_Intro(ls.all),1);
    new_line;
    Standard_Solution_Strings.Parse_Intro(s,ind,t,m,fail);
    if fail then
      put_line("an exception was raised when parsing the intro");
    else
      put_line("after parsing the string :");
      put(" t : "); put(t); new_line;
      put(" m : "); put(m,1); new_line;
    end if;
  end Standard_Test_Write_Intro;

  procedure DoblDobl_Test_Write_Intro
              ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution ) is

    use DoblDobl_Complex_Numbers;
    s : constant string := DoblDobl_Solution_Strings.Write_Intro(ls.all);
    t : Complex_Number;
    m : integer32;
    ind : integer := s'first;
    fail : boolean;

  begin
    put_line("Introduction written as string :");
    put(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(DoblDobl_Solution_Strings.Length_Intro(ls.all),1);
    new_line;
    DoblDobl_Solution_Strings.Parse_Intro(s,ind,t,m,fail);
    if fail then
      put_line("an exception was raised when parsing the intro");
    else
      put_line("after parsing the string :");
      put(" t : "); put(t); new_line;
      put(" m : "); put(m,1); new_line;
    end if;
  end DoblDobl_Test_Write_Intro;

  procedure QuadDobl_Test_Write_Intro
              ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution ) is

    use QuadDobl_Complex_Numbers;
    s : constant string := QuadDobl_Solution_Strings.Write_Intro(ls.all);
    t : Complex_Number;
    m : integer32;
    ind : integer := s'first;
    fail : boolean;

  begin
    put_line("Introduction written as string :");
    put(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(QuadDobl_Solution_Strings.Length_Intro(ls.all),1);
    new_line;
    QuadDobl_Solution_Strings.Parse_Intro(s,ind,t,m,fail);
    if fail then
      put_line("an exception was raised when parsing the intro");
    else
      put_line("after parsing the string :");
      put(" t : "); put(t); new_line;
      put(" m : "); put(m,1); new_line;
    end if;
  end QuadDobl_Test_Write_Intro;

  procedure Multprec_Test_Write_Intro
              ( ls : in Multprec_Complex_Solutions.Link_to_Solution ) is

    use Multprec_Complex_Numbers;
    s : constant string := Multprec_Solution_Strings.Write_Intro(ls.all);
    t : Complex_Number;
    m : integer32;
    ind : integer := s'first;
    fail : boolean;

  begin
    put_line("Introduction written as string :");
    put(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(Multprec_Solution_Strings.Length_Intro(ls.all),1);
    new_line;
    Multprec_Solution_Strings.Parse_Intro(s,ind,t,m,fail);
    if fail then
      put_line("an exception was raised when parsing the intro");
    else
      put_line("after parsing the string :");
      put(" t : "); put(t); new_line;
      put(" m : "); put(m,1); new_line;
    end if;
  end Multprec_Test_Write_Intro;

  procedure Standard_Test_Write_Vector
              ( ls : in Standard_Complex_Solutions.Link_to_Solution ) is

    use Standard_Complex_Solutions;
    s : constant string := Standard_Solution_Strings.Write_Vector(ls.all);
    v : Standard_Complex_Vectors.Vector(ls.v'range);
    fail : boolean;
    ns : Solution(ls.n);
    ind : integer := s'first;

  begin
    put_line("Solution vector written as string :");
    put(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(Standard_Solution_Strings.Length_Vector(ls.all),1);
    new_line;
    Standard_Solution_Strings.Parse_Vector(s,ind,natural32(ls.n),v,fail);
    if fail then
      put_line("an exception was raised when parsing the vector");
    else
      ns.v := v;
      put_line("The solution vector after parsing : ");
      put_vector(ns);
    end if;
  end Standard_Test_Write_Vector;

  procedure DoblDobl_Test_Write_Vector
              ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution ) is

    use DoblDobl_Complex_Solutions;
    s : constant string := DoblDobl_Solution_Strings.Write_Vector(ls.all);
    v : DoblDobl_Complex_Vectors.Vector(ls.v'range);
    fail : boolean;
    ns : Solution(ls.n);
    ind : integer := s'first;

  begin
    put_line("Solution vector written as string :");
    put(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(DoblDobl_Solution_Strings.Length_Vector(ls.all),1);
    new_line;
    DoblDobl_Solution_Strings.Parse_Vector(s,ind,natural32(ls.n),v,fail);
    if fail then
      put_line("an exception was raised when parsing the vector");
    else
      ns.v := v;
      put_line("The solution vector after parsing : ");
      put_vector(ns);
    end if;
  end DoblDobl_Test_Write_Vector;

  procedure QuadDobl_Test_Write_Vector
              ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution ) is

    use QuadDobl_Complex_Solutions;
    s : constant string := QuadDobl_Solution_Strings.Write_Vector(ls.all);
    v : QuadDobl_Complex_Vectors.Vector(ls.v'range);
    fail : boolean;
    ns : Solution(ls.n);
    ind : integer := s'first;

  begin
    put_line("Solution vector written as string :");
    put(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(QuadDobl_Solution_Strings.Length_Vector(ls.all),1);
    new_line;
    QuadDobl_Solution_Strings.Parse_Vector(s,ind,natural32(ls.n),v,fail);
    if fail then
      put_line("an exception was raised when parsing the vector");
    else
      ns.v := v;
      put_line("The solution vector after parsing : ");
      put_vector(ns);
    end if;
  end QuadDobl_Test_Write_Vector;

  procedure Multprec_Test_Write_Vector
              ( ls : in Multprec_Complex_Solutions.Link_to_Solution ) is

    use Multprec_Complex_Solutions;
    s : constant string := Multprec_Solution_Strings.Write_Vector(ls.all);
    v : Multprec_Complex_Vectors.Vector(ls.v'range);
    fail : boolean;
    ns : Solution(ls.n);
    ind : integer := s'first;

  begin
    put_line("Solution vector written as string :");
    put(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(Multprec_Solution_Strings.Length_Vector(ls.all),1);
    new_line;
    Multprec_Solution_Strings.Parse_Vector(s,ind,natural32(ls.n),v,fail);
    if fail then
      put_line("an exception was raised when parsing the vector");
    else
      ns.v := v;
      put_line("The solution vector after parsing : ");
      put_vector(ns);
    end if;
  end Multprec_Test_Write_Vector;

  procedure Standard_Test_Write_Diagnostics
              ( ls : in Standard_Complex_Solutions.Link_to_Solution ) is

    s : constant string
      := Standard_Solution_Strings.Write_Diagnostics(ls.all);
    err,rco,res : double_float;
    ind : integer := s'first;
    fail : boolean;

  begin
    put_line("Diagnostics written as string :");
    put_line(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(Standard_Solution_Strings.Length_Diagnostics,1);
    new_line;
    put_line("parsing the diagnostics ...");
    Standard_Solution_Strings.Parse_Diagnostics(s,ind,err,rco,res,fail);
    if fail then
      put_line("an exception occurred when parsing diagnostics");
    else
      put("  err : "); put(err,3);
      put("  rco : "); put(rco,3);
      put("  res : "); put(res,3); new_line;
    end if;
  end Standard_Test_Write_Diagnostics;

  procedure DoblDobl_Test_Write_Diagnostics
              ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution ) is

    s : constant string
      := DoblDobl_Solution_Strings.Write_Diagnostics(ls.all);
    err,rco,res : double_double;
    ind : integer := s'first;
    fail : boolean;

  begin
    put_line("Diagnostics written as string :");
    put_line(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(DoblDobl_Solution_Strings.Length_Diagnostics,1);
    new_line;
    put_line("parsing the diagnostics ...");
    DoblDobl_Solution_Strings.Parse_Diagnostics(s,ind,err,rco,res,fail);
    if fail then
      put_line("an exception occurred when parsing diagnostics");
    else
      put("  err : "); put(err,3);
      put("  rco : "); put(rco,3);
      put("  res : "); put(res,3); new_line;
    end if;
  end DoblDobl_Test_Write_Diagnostics;

  procedure QuadDobl_Test_Write_Diagnostics
              ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution ) is

    s : constant string
      := QuadDobl_Solution_Strings.Write_Diagnostics(ls.all);
    err,rco,res : quad_double;
    ind : integer := s'first;
    fail : boolean;

  begin
    put_line("Diagnostics written as string :");
    put_line(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(QuadDobl_Solution_Strings.Length_Diagnostics,1);
    new_line;
    put_line("parsing the diagnostics ...");
    QuadDobl_Solution_Strings.Parse_Diagnostics(s,ind,err,rco,res,fail);
    if fail then
      put_line("an exception occurred when parsing diagnostics");
    else
      put("  err : "); put(err,3);
      put("  rco : "); put(rco,3);
      put("  res : "); put(res,3); new_line;
    end if;
  end QuadDobl_Test_Write_Diagnostics;

  procedure Multprec_Test_Write_Diagnostics
              ( ls : in Multprec_Complex_Solutions.Link_to_Solution ) is

    s : constant string
      := Multprec_Solution_Strings.Write_Diagnostics(ls.all);
    err,rco,res : Floating_Number;
    ind : integer := s'first;
    fail : boolean;

  begin
    put_line("Diagnostics written as string :");
    put_line(s);
    put("Length of string : "); put(natural32(s'last),1); put(" = ");
    put(Multprec_Solution_Strings.Length_Diagnostics,1);
    new_line;
    put_line("parsing the diagnostics ...");
    Multprec_Solution_Strings.Parse_Diagnostics(s,ind,err,rco,res,fail);
    if fail then
      put_line("an exception occurred when parsing diagnostics");
    else
      put("  err : "); put(err,3);
      put("  rco : "); put(rco,3);
      put("  res : "); put(res,3); new_line;
    end if;
  end Multprec_Test_Write_Diagnostics;

  procedure Main is

    st_sols : Standard_Complex_Solutions.Solution_List;
    dd_sols : DoblDobl_Complex_Solutions.Solution_List;
    qd_sols : QuadDobl_Complex_Solutions.Solution_List;
    mp_sols : Multprec_Complex_Solutions.Solution_List;
    st_ls : Standard_Complex_Solutions.Link_to_Solution;
    dd_ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    qd_ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    mp_ls : Multprec_Complex_Solutions.Link_to_Solution;
    len : natural32;
    ans : character;

  begin
    new_line;
    put_line("MENU to test the writing and parsing of solutions :");
    put_line("  1. test with standard double complex numbers;");
    put_line("  2. test with double double complex numbers;");
    put_line("  3. test with quad double complex numbers;");
    put_line("  4. test with multiprecision double complex numbers.");
    put("Type 1, 2, 3, or 4 to select precision : ");
    Ask_Alternative(ans,"1234");
    new_line;
    case ans is 
      when '1' =>
        Read(st_sols);
        len := Standard_Complex_Solutions.Length_Of(st_sols);
        new_line;
        if len = 0 then
          put("No solutions read.");
        else
          put("Read "); put(len,1); put_line(" solutions.");
          st_ls := Standard_Complex_Solutions.Head_Of(st_sols);
          put_line("The first solution : "); put(st_ls.all); new_line;
          Standard_Test_Write(st_ls);
          Standard_Test_Write_Intro(st_ls);
          Standard_Test_Write_Vector(st_ls);
          Standard_Test_Write_Diagnostics(st_ls);
        end if;
      when '2' =>
        Read(dd_sols);
        len := DoblDobl_Complex_Solutions.Length_Of(dd_sols);
        new_line;
        if len = 0 then
          put("No solutions read.");
        else
          put("Read "); put(len,1); put_line(" solutions.");
          dd_ls := DoblDobl_Complex_Solutions.Head_Of(dd_sols);
          put_line("The first solution : "); put(dd_ls.all); new_line;
          DoblDobl_Test_Write(dd_ls);
          DoblDobl_Test_Write_Intro(dd_ls);
          DoblDobl_Test_Write_Vector(dd_ls);
          DoblDobl_Test_Write_Diagnostics(dd_ls);
        end if;
      when '3' =>
        Read(qd_sols);
        len := QuadDobl_Complex_Solutions.Length_Of(qd_sols);
        new_line;
        if len = 0 then
          put("No solutions read.");
        else
          put("Read "); put(len,1); put_line(" solutions.");
          qd_ls := QuadDobl_Complex_Solutions.Head_Of(qd_sols);
          put_line("The first solution : "); put(qd_ls.all); new_line;
          QuadDobl_Test_Write(qd_ls);
          QuadDobl_Test_Write_Intro(qd_ls);
          QuadDobl_Test_Write_Vector(qd_ls);
          QuadDobl_Test_Write_Diagnostics(qd_ls);
        end if;
      when '4' =>
        Read(mp_sols);
        len := Multprec_Complex_Solutions.Length_Of(mp_sols);
        new_line;
        if len = 0 then
          put("No solutions read.");
        else
          put("Read "); put(len,1); put_line(" solutions.");
          mp_ls := Multprec_Complex_Solutions.Head_Of(mp_sols);
          put_line("The first solution : "); put(mp_ls.all); new_line;
          Multprec_Test_Write(mp_ls);
          Multprec_Test_Write_Intro(mp_ls);
          Multprec_Test_Write_Vector(mp_ls);
          Multprec_Test_Write_Diagnostics(mp_ls);
        end if;
      when others => null;
    end case;
  end Main;

end Test_Solution_Strings;
