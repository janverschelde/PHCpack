with Ada.Text_IO;                       use Ada.Text_IO;
with String_Splitters;
with Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Symbol_Table;
with Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;
with Standard_Random_Laurentials;       use Standard_Random_Laurentials;
with Real_Powered_Series_IO;
with Test_Real_Powered_Series;
with Real_Powered_Homotopy_io;

package body Test_Real_Powered_Homotopy is

  procedure Random_Series_Vector
              ( size : in integer32;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwt : out Standard_Floating_VecVecs.VecVec ) is
  begin
    for i in cff'range loop
      declare
        icf : Standard_Complex_Vectors.Vector(0..size);
        ipw : Standard_Floating_Vectors.Vector(1..size);
      begin
        Test_Real_Powered_Series.Random_Series(size,icf,ipw);
        cff(i) := new Standard_Complex_Vectors.Vector'(icf);
        pwt(i) := new Standard_Floating_Vectors.Vector'(ipw);
      end;
    end loop;
  end Random_Series_Vector;

  procedure Test_Random_Polynomial ( nbr,nvr,size : in integer32 ) is

    q : Standard_Complex_Laurentials.Poly;
    cff : Standard_Complex_VecVecs.VecVec(1..nbr);
    pwt : Standard_Floating_VecVecs.VecVec(1..nbr);
    ans : character;
    file : file_type;
    name : String_Splitters.Link_to_String;
  
  begin
    new_line;
    put("-> generating "); put(nbr,1); put(" monomials in ");
    put(nvr,1); put_line(" variables ...");
    q := Random_Laurent_Polynomial(natural32(nvr),natural32(nbr),-9,+9);
    Standard_Complex_Laurentials_io.put(q); new_line;
    Random_Series_Vector(size,cff,pwt);
    for i in cff'range loop
      put("random series "); put(i,1); put_line(" :");
      Real_Powered_Series_IO.put_line(cff(i).all,pwt(i).all);
    end loop;
    new_line;
    put("new line format ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y'
     then Real_Powered_Homotopy_IO.put_line(q,cff,pwt);
     else Real_Powered_Homotopy_IO.put(q,cff,pwt);
    end if;
    new_line;
    put_line("Reading file name for output ...");
    Communications_with_User.Read_Name_and_Create_File(file,name);
    if ans = 'y'
     then Real_Powered_Homotopy_IO.put_line(file,q,cff,pwt);
     else Real_Powered_Homotopy_IO.put(file,q,cff,pwt);
    end if;
    Close(file);
    new_line;
    put_line("Closed file.  Reopening again for reading ...");
    new_line;
    Communications_with_User.Open_Input_File(file,name.all);
    declare
      q2 : Standard_Complex_Laurentials.Poly;
      c2 : Standard_Complex_VecVecs.VecVec(1..nbr);
      p2 : Standard_Floating_VecVecs.VecVec(1..nbr);
    begin
      Real_Powered_Homotopy_IO.get(file,nvr,size,q2,c2,p2,vrblvl=>2);
      for i in 1..nbr loop
        put("-> power series "); put(i,1); put_line(" :");
        Real_Powered_Series_IO.put_line(cff(i).all,pwt(i).all);
      end loop;
    end;
  end Test_Random_Polynomial;

  procedure Test_String_Polynomial ( nbr,nvr,size : in integer32 ) is

    q : Standard_Complex_Laurentials.Poly;
    cff : Standard_Complex_VecVecs.VecVec(1..nbr);
    pwt : Standard_Floating_VecVecs.VecVec(1..nbr);
    nsyms : constant integer32 := integer32(Symbol_Table.number);
  
  begin
    new_line;
    put("-> generating "); put(nbr,1); put(" monomials in ");
    put(nvr,1); put_line(" variables ...");
    q := Random_Laurent_Polynomial(natural32(nvr),natural32(nbr),-9,+9);
    Standard_Complex_Laurentials_io.put(q); new_line;
    Random_Series_Vector(size,cff,pwt);
    for i in cff'range loop
      put("random series "); put(i,1); put_line(" :");
      Real_Powered_Series_IO.put_line(cff(i).all,pwt(i).all);
    end loop;
    put("Number of symbols : "); put(nsyms,1); new_line;
    if nsyms = 0 then
      declare
        asb : constant Symbol_Table.Array_of_Symbols(1..nvr)
            := Symbol_Table.Standard_Symbols(nvr);
      begin
        Symbol_Table.Init(asb);
      end;
    end if;
    declare
      s : constant string
        := Real_Powered_Homotopy_IO.to_string(q,cff,pwt,vrblvl=>1);
    begin
      new_line;
      put_line("The string representation of the polynomial :");
      put_line(s);
      new_line;
      put_line("Parsing the string for the data ...");
      declare
        m : constant integer32
          := Real_Powered_Homotopy_IO.number_of_terms(s,2);
      begin
        put("number of terms in string representation : ");
        put(m,1);
        if m /= nbr then
          put(" /= "); put(nbr,1); put_line(", bug!");
        else
          put_line(", okay.");
          declare
            q2 : Standard_Complex_Laurentials.Poly;
            c2 : Standard_Complex_VecVecs.VecVec(1..nbr);
            p2 : Standard_Floating_VecVecs.VecVec(1..nbr);
          begin
            Real_Powered_Homotopy_IO.parse_string(s,nvr,q2,c2,p2,vrblvl=>2);
            for i in cff'range loop
              put("parsed series "); put(i,1); put_line(" :");
              Real_Powered_Series_IO.put_line(c2(i).all,p2(i).all);
            end loop;
            put_line("parsed Laurent polynomial :");
            Standard_Complex_Laurentials_io.put(q2); new_line;
          end;
        end if;
      end;
    end;
  end Test_String_Polynomial;

  procedure Main is

    nbr,nvr,size : integer32 := 0;
    ans : character;

  begin
    put("Give the number of variables : "); get(nvr);
    put("Give the number of monomials : "); get(nbr);
    put("Give the size of each series : "); get(size);
    put("Test output to string ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_String_Polynomial(nbr,nvr,size);
     else Test_Random_Polynomial(nbr,nvr,size);
    end if;
  end Main;

end Test_Real_Powered_Homotopy;
