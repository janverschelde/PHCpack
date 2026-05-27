with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Random_Vectors;
with Double_Real_Powered_Series;
with Real_Powered_Series_IO;

package body Test_Real_Powered_Series is

  procedure Write ( cff : in Standard_Complex_Vectors.Vector;
                    pwt : in Standard_Floating_Vectors.Vector;
                    t : in character := 't' ) is
  begin
    put(cff(0)); new_line;
    for k in pwt'range loop
      put(cff(k)); put(" * "); put(t); put("**");
      put(pwt(k),1,14,3); new_line;
    end loop;
  end Write;

  procedure Random_Series
              ( size : in integer32;
                cff : out Standard_Complex_Vectors.Vector;
                pwt : out Standard_Floating_Vectors.Vector ) is
  begin
    cff := Standard_Random_Vectors.Random_Vector(0,size);
    pwt := Standard_Random_Vectors.Random_Vector(1,size);
    for i in pwt'range loop
      pwt(i) := abs(pwt(i));
    end loop;
    Double_Real_Powered_Series.sort(pwt,cff);
  end Random_Series;

  procedure Test_String_Series ( size : in integer32 ) is

    cff : Standard_Complex_Vectors.Vector(0..size);
    pwt : Standard_Floating_Vectors.Vector(1..size);
    ptr_cff : Standard_Complex_Vectors.Link_to_Vector;
    ptr_pwt : Standard_Floating_Vectors.Link_to_Vector;
    tol : constant double_float := 1.0e-14;
    difference : Complex_Number;
    fail : boolean := false;

  begin
    Random_Series(size,cff,pwt);
    put_line("A random series : "); Write(cff,pwt);
    declare
      s : constant string := Real_Powered_Series_IO.to_string(cff,pwt);
    begin
      new_line;
      put_line("The string format of the series :");
      put_line(s);
      new_line;
      put_line("Parsing the string for the series ...");
      Real_Powered_Series_IO.parse_string(s,ptr_cff,ptr_pwt,vrblvl=>1);
      put_line("The series parsed from string :");
      Write(ptr_cff.all,ptr_pwt.all);
      difference := ptr_cff(0) - cff(0);
      if AbsVal(difference) > tol then
        put_line("Constants are not equal!");
        put("difference : "); put(difference); new_line;
        fail := true;
      end if;
      for i in pwt'range loop
        difference := ptr_cff(i) - cff(i);
        if AbsVal(difference) > tol then
          put("coefficients "); put(i,1); put_line(" are not equal!");
          put("difference : "); put(difference); new_line;
          fail := true;
        end if;
        exit when fail;
        if abs(ptr_pwt(i) - pwt(i)) > tol then
          put("powers "); put(i,1); put_line(" are not equal!");
          put("difference : "); put(ptr_pwt(i)-pwt(i)); new_line;
          fail := true;
        end if;
        exit when fail;
      end loop;
      if fail
       then put_line("Parsing of string failed!");
       else put_line("Parsing of string succeeded.");
      end if;
    end;
  end Test_String_Series;

  procedure Write_Series_to_File
              ( cff : in Standard_Complex_Vectors.Vector;
                pwt : in Standard_Floating_Vectors.Vector;
                newline : in boolean := true ) is

  -- DESCRPTION :
  --   Prompts for a file name and writes the series to file,
  --   using the new line format if newline.

    file : file_type;
 
  begin
    Communications_with_User.Read_Name_and_Create_File(file);
    if newline
     then Real_Powered_Series_IO.put_line(file,cff,pwt);
     else Real_Powered_Series_IO.put(file,cff,pwt); new_line(file);
    end if;
    close(file);
  end Write_Series_to_File;

  procedure Read_Series_from_File 
              ( cff : out Standard_Complex_Vectors.Vector;
                pwt : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Prompts for a file name and reads a series from file,
  --   returned in cff and pwt.

    file : file_type;

  begin
    new_line;
    put_line("Reading the series from file ...");
    Communications_with_User.Read_Name_and_Open_File(file);
    Real_Powered_Series_IO.get(file,cff,pwt,vrblvl=>1);
    close(file);
  end Read_Series_from_File;

  procedure Test_Series_IO ( size : in integer32 ) is

    cff1,cff2 : Standard_Complex_Vectors.Vector(0..size);
    pwt1,pwt2 : Standard_Floating_Vectors.Vector(1..size);
    ans : character;
    tol : constant double_float := 1.0E-14;
    fail : boolean := false;
    difference : complex_number;

  begin
    Random_Series(size,cff1,pwt1);
    put_line("A random series : "); Write(cff1,pwt1);
    new_line;
    put("Write in new_line format ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y' then
      Real_Powered_Series_IO.put_line(cff1,pwt1);
      Write_Series_to_File(cff1,pwt1);
    else
      Real_Powered_Series_IO.put(cff1,pwt1); new_line;
      Write_Series_to_File(cff1,pwt1,false);
    end if;
    Read_Series_from_File(cff2,pwt2);
    new_line;
    put_line("The series read from file :"); Write(cff2,pwt2);
    difference := cff2(0) - cff1(0);
    if AbsVal(difference) > tol then
      put_line("Constants are not equal!");
      put("difference : "); put(difference); new_line;
      fail := true;
    end if;
    for i in pwt1'range loop
      difference := cff2(i) - cff1(i);
      if AbsVal(difference) > tol then
        put("coefficients "); put(i,1); put_line(" are not equal!");
        put("difference : "); put(difference); new_line;
        fail := true;
      end if;
      exit when fail;
      if abs(pwt2(i) - pwt1(i)) > tol then
        put("powers "); put(i,1); put_line(" are not equal!");
        put("difference : "); put(pwt2(i)-pwt1(i)); new_line;
        fail := true;
      end if;
      exit when fail;
    end loop;
    if fail
     then put_line("Reading from file failed!");
     else put_line("Reading from file succeeded.");
    end if;
  end Test_Series_IO;

  procedure Test_Random_Series ( size : in integer32 ) is

    cff : Standard_Complex_Vectors.Vector(0..size);
    pwt : Standard_Floating_Vectors.Vector(1..size);

  begin
    put("-> generating a series of size "); put(size,1); put_line(" ...");
    Random_Series(size,cff,pwt);
    Write(cff,pwt);
  end Test_Random_Series;

  procedure Test_Addition ( size : in integer32 ) is

    acf,bcf : Standard_Complex_Vectors.Vector(0..size);
    apw,bpw : Standard_Floating_Vectors.Vector(1..size);
    sumcf : Standard_Complex_Vectors.Vector(0..2*size);
    sumpw : Standard_Floating_Vectors.Vector(1..2*size);
    difcf : Standard_Complex_Vectors.Vector(0..3*size);
    difpw : Standard_Floating_Vectors.Vector(1..3*size);
    equ : boolean;

  begin
    put("-> adding two series of size "); put(size,1); put_line(" ...");
    Random_Series(size,acf,apw);
    Random_Series(size,bcf,bpw);
    put_line("The sum of "); Write(acf,apw);
    put_line("and"); Write(bcf,bpw);
    put_line("is");
    Double_Real_Powered_Series.Add(acf,bcf,apw,bpw,sumcf,sumpw);
    Write(sumcf,sumpw);
    Double_Real_Powered_Series.Sub(sumcf,bcf,sumpw,bpw,difcf,difpw);
    put_line("After subtracting the second series from the sum :");
    Write(difcf(0..size),difpw(1..size));
    put_line("the first series :");
    Write(acf,apw);
    equ := Double_Real_Powered_Series.Equal
             (difcf(0..size),acf,difpw(1..size),apw);
    put("Component wise equal ? ");
    if equ
     then put_line("Yes.");
     else put_line("No, bug!?");
    end if;
  end Test_Addition;

  procedure Test_Multiplication ( size : in integer32 ) is

    acf,bcf : Standard_Complex_Vectors.Vector(0..size);
    apw,bpw : Standard_Floating_Vectors.Vector(1..size);
    prdsize : constant integer32 := (size+1)*(size+1) - 1;
    prdcf : Standard_Complex_Vectors.Vector(0..prdsize);
    prdpw : Standard_Floating_Vectors.Vector(1..prdsize);

  begin
    put("-> multiplying two series of size "); put(size,1); put_line(" ...");
    Random_Series(size,acf,apw);
    Random_Series(size,bcf,bpw);
    put_line("The product of "); Write(acf,apw);
    put_line("and"); Write(bcf,bpw);
    put_line("is");
    Double_Real_Powered_Series.Mul(acf,bcf,apw,bpw,prdcf,prdpw);
    Write(prdcf,prdpw);
  end Test_Multiplication;

  procedure Test_Inverse ( size : in integer32 ) is

    acf,bcf : Standard_Complex_Vectors.Vector(0..size);
    apw,bpw : Standard_Floating_Vectors.Vector(1..size);
    prdsize : constant integer32 := (size+1)*(size+1) - 1;
    prdcf : Standard_Complex_Vectors.Vector(0..prdsize);
    prdpw : Standard_Floating_Vectors.Vector(1..prdsize);

  begin
    Random_Series(size,acf,apw);
    for i in apw'range loop
      apw(i) := apw(i) + 1.0;
    end loop;
    put("-> the inverse of a random series of size "); put(size,1);
    put_line(" :"); Write(acf,apw);
    put_line("is");
    Double_Real_Powered_Series.Inv(acf,apw,bcf,bpw);
    Write(bcf,bpw);
    put_line("the product of a random series with its inverse :");
    Double_Real_Powered_Series.Mul(acf,bcf,apw,bpw,prdcf,prdpw);
    Write(prdcf(0..size),prdpw(1..size));
  end Test_Inverse;

  procedure Test_Division ( size : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random series of the given size,
  --   tests their division.

    acf,bcf : Standard_Complex_Vectors.Vector(0..size);
    apw,bpw : Standard_Floating_Vectors.Vector(1..size);
    qsize : constant integer32 := (size+1)*(size+1) - 1;
    qcf : Standard_Complex_Vectors.Vector(0..qsize);
    qpw : Standard_Floating_Vectors.Vector(1..qsize);
    psize : constant integer32 := (size+1)*(qsize+1) - 1;
    pcf : Standard_Complex_Vectors.Vector(0..psize);
    ppw : Standard_Floating_Vectors.Vector(1..psize);
    equ : boolean;

  begin
    put("-> dividing two series of size "); put(size,1); put_line(" ...");
    Random_Series(size,acf,apw);
    Random_Series(size,bcf,bpw);
    for i in apw'range loop
      apw(i) := apw(i) + 1.0;
      bpw(i) := bpw(i) + 1.0;
    end loop;
    put_line("The division of "); Write(acf,apw);
    put_line("by"); Write(bcf,bpw);
    put_line("is");
    Double_Real_Powered_Series.Div(acf,bcf,apw,bpw,qcf,qpw);
    Write(qcf(0..size),qpw(1..size));
    put_line("multiplying the quotient with the second series ...");
    Double_Real_Powered_Series.Mul(bcf,qcf,bpw,qpw,pcf,ppw);
    Write(pcf(0..size),ppw(1..size));
    put_line("The first series :");
    Write(acf,apw);
    equ := Double_Real_Powered_Series.Equal(pcf(0..size),acf,ppw(1..size),apw);
    put("Component wise equal ? ");
    if equ
     then put_line("Yes.");
     else put_line("No, bug!?");
    end if;
  end Test_Division;

  procedure Test_Series_IO is

    deg : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the size of the series : "); get(deg);
    new_line;
    put("-> generating a random series of size "); put(deg,1);
    put_line(" ...");
    new_line;
    put("Test output to string ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_String_Series(deg);
     else Test_Series_IO(deg);
    end if;
  end Test_Series_IO;

  procedure Test_Series_Arithmetic is

    ans : character;
    size : integer32 := 0;

  begin
    new_line;
    put_line("MENU to test real powered series arithmetic ...");
    put_line("  0. generate a random series");
    put_line("  1. test addition and subtraction");
    put_line("  2. test multiplication");
    put_line("  3. test multiplicative inverse");
    put_line("  4. test division");
    put("Type 0, 1, 2, 3, or 4 to select a test : ");
    Communications_with_User.Ask_Alternative(ans,"01234");
    new_line;
    put("Give the size of the series : "); get(size);
    case ans is
      when '0' => Test_Random_Series(size);
      when '1' => Test_Addition(size);
      when '2' => Test_Multiplication(size);
      when '3' => Test_Inverse(size);
      when '4' => Test_Division(size);
      when others => null;
    end case;
  end Test_Series_Arithmetic;

  procedure Main is
  begin
    new_line;
    put_line("Testing the input/output procedures ...");
    Test_Series_IO;
    new_line;
    put_line("Testing the arithmetic ...");
    Test_Series_Arithmetic;
  end Main;
 
end Test_Real_Powered_Series;
