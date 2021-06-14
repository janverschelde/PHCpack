with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Interfaces.C;
with C_integer_io;                       use C_Integer_io;
with C_Integer_Arrays;                   use C_Integer_Arrays;
with C_Integer_Arrays_io;                use C_Integer_Arrays_io;
with C_Double_io;                        use C_Double_io;
with C_Double_Arrays;                    use C_Double_Arrays;
with C_Double_Arrays_io;                 use C_Double_Arrays_io;
with C_to_Ada_Arrays;                    use C_to_Ada_Arrays;

procedure ts_arrays is

-- DESCRIPTION :
--   This package allows interactive testing of the conversions
--   between the various array types in C and Ada.

  procedure Test_C_to_Ada_Integer_Array is

  -- DESCRIPTION :
  --   Reads in an array of C integers and writes the converted Ada array.

    n : integer32 := 0;

  begin
    put("Give the number of elements in the array : "); get(n);
    declare
      ca : C_Integer_Array(0..Interfaces.C.size_T(n-1));
    begin
      put("Give "); put(n,1); put(" integers : ");
      for i in ca'range loop
        get(ca(i));
      end loop;
      put("Your integer array : "); put(integer(n),ca); new_line;
      declare
        aa : constant Standard_Integer_Vectors.Vector := Convert(ca);
      begin
        put("The Ada array : "); put(aa); new_line;
      end;
    end;
  end Test_C_to_Ada_Integer_Array;

  procedure Test_C_Double_to_Ada_Complex_Array is

  -- DESCRIPTION :
  --   Reads in an array of C doubles and writes the converted Ada
  --   array of complex numbers.

    n : integer32 := 0;

  begin
    loop
      put("Give the number of elements in the array : "); get(n);
      exit when n mod 2 = 0;
      put("  "); put(n,1);
      put_line(" is not a multiple of 2.  Please try again...");
    end loop;
    declare
      ca : C_Double_Array(0..Interfaces.C.size_T(n-1));
    begin
      put("Give "); put(n,1); put(" doubles : ");
      for i in ca'range loop
        get(ca(i));
      end loop;
      put_line("Your double array : "); put(integer(n),ca);
      declare
        aa : constant Standard_Complex_Vectors.Vector := Convert(ca);
      begin
        put_line("The Ada complex array : ");
        put_line(aa);
      end;
    end;
  end Test_C_Double_to_Ada_Complex_Array;

  procedure Test_Ada_to_C_Integer_Array is

    a,b : integer32 := 0;

  begin
    put("Give lower index : "); get(a);
    put("Give upper index : "); get(b);
    declare
      aa : Standard_Integer_Vectors.Vector(a..b);
    begin
      put("Give "); put(b-a+1,1); put(" integers : ");
      for i in aa'range loop
        aa(i) := 0; get(aa(i));
      end loop;
      put("Your integer array : "); put(aa); new_line;
      declare
        ca : constant C_Integer_Array := Convert(aa);
      begin
        put("The C integer array : "); put(ca'length,ca); new_line;
      end;
    end;
  end Test_Ada_to_C_Integer_Array;

  procedure Test_Ada_Complex_to_C_Double_Array is

    a,b : integer32 := 0;

  begin
    put("Give lower index : "); get(a);
    put("Give upper index : "); get(b);
    declare
      aa : Standard_Complex_Vectors.Vector(a..b);
    begin
      put("Give "); put(b-a+1,1); put(" complex numbers : ");
      for i in aa'range loop
        get(aa(i));
      end loop;
      put_line("Your complex array : "); put_line(aa); 
      declare
        ca : constant C_Double_Array := Convert(aa);
      begin
        put_line("The C double array : "); put(ca'length,ca);
      end;
    end;
  end Test_Ada_Complex_to_C_Double_Array;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Conversion between C and Ada arrays.");
    new_line;
    put_line("MENU to convert arrays from C to Ada and from Ada to C :");
    put_line("  1. convert C integer array into Ada integer array;");
    put_line("  2. convert Ada integer array into C integer array;");
    put_line("  3. convert C double array into Ada complex array;");
    put_line("  4. convert Ada complex array into C double array;");
    put("Type 1, 2, 3 or 4 to choose : "); Ask_Alternative(ans,"1234");
    new_line;
    case ans is
      when '1' => Test_C_to_Ada_Integer_Array;
      when '2' => Test_Ada_to_C_Integer_Array;
      when '3' => Test_C_Double_to_Ada_Complex_Array;
      when '4' => Test_Ada_Complex_to_C_Double_Array;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_arrays;
