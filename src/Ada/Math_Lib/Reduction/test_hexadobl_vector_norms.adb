with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Hexa_Double_Numbers_io;             use Hexa_Double_Numbers_io;
with HexaDobl_Complex_Numbers;           use HexaDobl_Complex_Numbers;
with Hexa_Double_Vectors;
with HexaDobl_Complex_Vectors;
with HexaDobl_Random_Vectors;
with Hexa_Double_Vector_Norms;           use Hexa_Double_Vector_Norms;
with HexaDobl_Complex_Vector_Norms;      use HexaDobl_Complex_Vector_Norms;

package body Test_HexaDobl_Vector_Norms is

  procedure Test_Real_Two_Norm is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    declare   
      v : constant Hexa_Double_Vectors.Vector(1..dim)
        := HexaDobl_Random_Vectors.Random_Vector(1,dim);
      vnrm : constant hexa_double := Norm2(v);
      w : Hexa_Double_Vectors.Vector(1..dim);
      wnrm : hexa_double;
    begin
      put("2-norm of a random vector : "); put(vnrm,3); new_line;
      for k in 1..dim loop
        w(k) := v(k)/vnrm;
      end loop;
      wnrm := Norm2(w);
      put_line("2-norm of the normalized vector :");
      put(wnrm); new_line;
    end;
  end Test_Real_Two_Norm;

  procedure Test_Real_Sum_Norm is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    declare   
      v : constant Hexa_Double_Vectors.Vector(1..dim)
        := HexaDobl_Random_Vectors.Random_Vector(1,dim);
      vnrm : constant hexa_double := Sum_Norm(v);
      w : Hexa_Double_Vectors.Vector(1..dim);
      wnrm : hexa_double;
    begin
      put("sum norm of a random vector : "); put(vnrm,3); new_line;
      for k in 1..dim loop
        w(k) := v(k)/vnrm;
      end loop;
      wnrm := Sum_Norm(w);
      put_line("sum norm of the normalized vector :");
      put(wnrm); new_line;
    end;
  end Test_Real_Sum_Norm;

  procedure Test_Real_Max_Norm is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    declare   
      v : constant Hexa_Double_Vectors.Vector(1..dim)
        := HexaDobl_Random_Vectors.Random_Vector(1,dim);
      vnrm : constant hexa_double := Max_Norm(v);
      w : Hexa_Double_Vectors.Vector(1..dim);
      wnrm : hexa_double;
    begin
      put("max norm of a random vector : "); put(vnrm,3); new_line;
      for k in 1..dim loop
        w(k) := v(k)/vnrm;
      end loop;
      wnrm := Max_Norm(w);
      put_line("max norm of the normalized vector :");
      put(wnrm); new_line;
    end;
  end Test_Real_Max_Norm;

  procedure Test_Complex_Two_Norm is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    declare   
      v : constant HexaDobl_Complex_Vectors.Vector(1..dim)
        := HexaDobl_Random_Vectors.Random_Vector(1,dim);
      vnrm : constant hexa_double := Norm2(v);
      w : HexaDobl_Complex_Vectors.Vector(1..dim);
      wnrm : hexa_double;
    begin
      put("2-norm of a random vector : "); put(vnrm,3); new_line;
      for k in 1..dim loop
        w(k) := v(k)/vnrm;
      end loop;
      wnrm := Norm2(w);
      put_line("2-norm of the normalized vector :");
      put(wnrm); new_line;
    end;
  end Test_Complex_Two_Norm;

  procedure Test_Complex_Sum_Norm is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    declare   
      v : constant HexaDobl_Complex_Vectors.Vector(1..dim)
        := HexaDobl_Random_Vectors.Random_Vector(1,dim);
      vnrm : constant hexa_double := Sum_Norm(v);
      w : HexaDobl_Complex_Vectors.Vector(1..dim);
      wnrm : hexa_double;
    begin
      put("sum norm of a random vector : "); put(vnrm,3); new_line;
      for k in 1..dim loop
        w(k) := v(k)/vnrm;
      end loop;
      wnrm := Sum_Norm(w);
      put_line("sum norm of the normalized vector :");
      put(wnrm); new_line;
    end;
  end Test_Complex_Sum_Norm;

  procedure Test_Complex_Max_Norm is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    declare   
      v : constant HexaDobl_Complex_Vectors.Vector(1..dim)
        := HexaDobl_Random_Vectors.Random_Vector(1,dim);
      vnrm : constant hexa_double := Max_Norm(v);
      w : HexaDobl_Complex_Vectors.Vector(1..dim);
      wnrm : hexa_double;
    begin
      put("max norm of a random vector : "); put(vnrm,3); new_line;
      for k in 1..dim loop
        w(k) := v(k)/vnrm;
      end loop;
      wnrm := Max_Norm(w);
      put_line("max norm of the normalized vector :");
      put(wnrm); new_line;
    end;
  end Test_Complex_Max_Norm;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test norms of hexa double vectors :");
    put_line("  1. 2-norm and normalization of a real random vector");
    put_line("  2. sum norm and normalization of a real random vector");
    put_line("  3. max norm and normalization of a real random vector");
    put_line("  4. 2-norm and normalization of a complex random vector");
    put_line("  5. sum norm and normalization of a complex random vector");
    put_line("  6. max norm and normalization of a complex random vector");
    put("Type 1, 2, 3, 4, 5, or 6 to select a test : ");
    Ask_Alternative(ans,"123456");
    new_line;
    case ans is
      when '1' => Test_Real_Two_Norm;
      when '2' => Test_Real_Sum_Norm;
      when '3' => Test_Real_Max_Norm;
      when '4' => Test_Complex_Two_Norm;
      when '5' => Test_Complex_Sum_Norm;
      when '6' => Test_Complex_Max_Norm;
      when others => null;
    end case;
  end Main;

end Test_HexaDobl_Vector_Norms;
