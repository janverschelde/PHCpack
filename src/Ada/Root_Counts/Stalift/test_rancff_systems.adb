with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Lists_of_Floating_Vectors_io;       use Lists_of_Floating_Vectors_io;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Arrays_of_Floating_Vector_Lists_io; use Arrays_of_Floating_Vector_Lists_io;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Floating_Integer_Convertors;        use Floating_Integer_Convertors;
with Supports_of_Polynomial_Systems;
with Random_Coefficient_Systems;

package body Test_RanCff_Systems is

  procedure Write_to_File ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    new_line;
    put("Write the random coefficient system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,p);
      close(file);
    end if;
  end Write_to_File;

  procedure Write_to_File ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    new_line;
    put("Write the random coefficient system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,p);
      close(file);
    end if;
  end Write_to_File;

  procedure Write_to_File ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    new_line;
    put("Write the random coefficient system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,p);
      close(file);
    end if;
  end Write_to_File;

  procedure Standard_Creation_of_Random_Polynomials is

    n,k : natural32 := 0;
    support : Lists_of_Integer_Vectors.List;
    fltsupp : Lists_of_Floating_Vectors.List;
    p : Standard_Complex_Polynomials.Poly;

  begin
    put("Give the number of variables : "); get(n);
    put("Give the number of monomials : "); get(k);
    put("Give "); put(k,1); put(" natural vectors of length ");
    put(n,1); put_line(" : "); get(n,k,support);
    p := Random_Coefficient_Systems.Create(n,support);
    put("The random polynomial :"); put_line(p);
    fltsupp := Convert(support);
    put_line("The floating point support list : "); put(fltsupp);
    Standard_Complex_Polynomials.Clear(p);
    p := Random_Coefficient_Systems.Create(n,fltsupp);
    put("The random polynomial created from a floating support list :");
    put_line(p);
  end Standard_Creation_of_Random_Polynomials;

  procedure DoblDobl_Creation_of_Random_Polynomials is

    n,k : natural32 := 0;
    support : Lists_of_Integer_Vectors.List;
    fltsupp : Lists_of_Floating_Vectors.List;
    p : DoblDobl_Complex_Polynomials.Poly;

  begin
    put("Give the number of variables : "); get(n);
    put("Give the number of monomials : "); get(k);
    put("Give "); put(k,1); put(" natural vectors of length ");
    put(n,1); put_line(" : "); get(n,k,support);
    p := Random_Coefficient_Systems.Create(n,support);
    put("The random polynomial :"); put_line(p);
    fltsupp := Convert(support);
    put_line("The floating point support list : "); put(fltsupp);
    DoblDobl_Complex_Polynomials.Clear(p);
    p := Random_Coefficient_Systems.Create(n,fltsupp);
    put("The random polynomial created from a floating support list :");
    put_line(p);
  end DoblDobl_Creation_of_Random_Polynomials;

  procedure QuadDobl_Creation_of_Random_Polynomials is

    n,k : natural32 := 0;
    support : Lists_of_Integer_Vectors.List;
    fltsupp : Lists_of_Floating_Vectors.List;
    p : QuadDobl_Complex_Polynomials.Poly;

  begin
    put("Give the number of variables : "); get(n);
    put("Give the number of monomials : "); get(k);
    put("Give "); put(k,1); put(" natural vectors of length ");
    put(n,1); put_line(" : "); get(n,k,support);
    p := Random_Coefficient_Systems.Create(n,support);
    put("The random polynomial :"); put_line(p);
    fltsupp := Convert(support);
    put_line("The floating point support list : "); put(fltsupp);
    QuadDobl_Complex_Polynomials.Clear(p);
    p := Random_Coefficient_Systems.Create(n,fltsupp);
    put("The random polynomial created from a floating support list :");
    put_line(p);
  end QuadDobl_Creation_of_Random_Polynomials;

  procedure Read_Supports
              ( dim : in natural32;
                s : out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is             
    n,k : natural32 := 0;

  begin
    for i in s'range loop
      put("Give the cardinality of support "); put(i,1); put(" : ");
      get(k);
      put("Give "); put(k,1); put(" natural vectors of length ");
      n := dim;
      put(dim,1); put_line(" : "); get(n,k,s(i));
    end loop;
  end Read_Supports;

  procedure Standard_Random_Polynomial_System ( n,r : integer32 ) is

    mix : Standard_Integer_Vectors.Vector(1..r);
    int_supports : Arrays_of_Integer_Vector_Lists.Array_of_Lists(1..r);
    flt_supports : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r);
    p : Standard_Complex_Poly_Systems.Poly_Sys(1..n);

  begin
    put("Give the type of mixture : "); get(mix);
    put("The type of mixture is "); put(mix); new_line;
    Read_Supports(natural32(n),int_supports);
    put_line("The supports : "); put(int_supports);
    p := Random_Coefficient_Systems.Create(natural32(n),mix,int_supports);
    put_line("The random coefficient system : "); put_line(p);
    flt_supports := Convert(int_supports);
    put_line("The floating supports : "); put(flt_supports);
    p := Random_Coefficient_Systems.Create(natural32(n),mix,flt_supports);
    put_line("The random coefficient system : "); put_line(p);
    Write_to_File(p);
  end Standard_Random_Polynomial_System;

  procedure DoblDobl_Random_Polynomial_System ( n,r : integer32 ) is

    mix : Standard_Integer_Vectors.Vector(1..r);
    int_supports : Arrays_of_Integer_Vector_Lists.Array_of_Lists(1..r);
    flt_supports : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r);
    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..n);

  begin
    put("Give the type of mixture : "); get(mix);
    put("The type of mixture is "); put(mix); new_line;
    Read_Supports(natural32(n),int_supports);
    put_line("The supports : "); put(int_supports);
    p := Random_Coefficient_Systems.Create(natural32(n),mix,int_supports);
    put_line("The random coefficient system : "); put_line(p);
    flt_supports := Convert(int_supports);
    put_line("The floating supports : "); put(flt_supports);
    p := Random_Coefficient_Systems.Create(natural32(n),mix,flt_supports);
    put_line("The random coefficient system : "); put_line(p);
    Write_to_File(p);
  end DoblDobl_Random_Polynomial_System;

  procedure QuadDobl_Random_Polynomial_System ( n,r : integer32 ) is

    mix : Standard_Integer_Vectors.Vector(1..r);
    int_supports : Arrays_of_Integer_Vector_Lists.Array_of_Lists(1..r);
    flt_supports : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r);
    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..n);

  begin
    put("Give the type of mixture : "); get(mix);
    put("The type of mixture is "); put(mix); new_line;
    Read_Supports(natural32(n),int_supports);
    put_line("The supports : "); put(int_supports);
    p := Random_Coefficient_Systems.Create(natural32(n),mix,int_supports);
    put_line("The random coefficient system : "); put_line(p);
    flt_supports := Convert(int_supports);
    put_line("The floating supports : "); put(flt_supports);
    p := Random_Coefficient_Systems.Create(natural32(n),mix,flt_supports);
    put_line("The random coefficient system : "); put_line(p);
    Write_to_File(p);
  end QuadDobl_Random_Polynomial_System;

  procedure Random_Polynomial_System ( precision : in character ) is

    n,r : integer32 := 0;

  begin
    put("Give the number of variables : "); get(n);
    put("Give the number of different supports : "); get(r);
    case precision is
      when '0' => Standard_Random_Polynomial_System(n,r);
      when '1' => DoblDobl_Random_Polynomial_System(n,r);
      when '2' => QuadDobl_Random_Polynomial_System(n,r);
      when others => null;
    end case;
  end Random_Polynomial_System;

  procedure Standard_Random_Coefficient_System is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(p);
    declare
      s : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p.all);
      n : constant integer32 := p'last;
      q : constant Standard_Complex_Poly_Systems.Poly_Sys(p'range)
        := Random_Coefficient_Systems.Create(natural32(n),s);
    begin
      put_line(q);
      Write_to_File(q);
    end;
  end Standard_Random_Coefficient_System;

  procedure DoblDobl_Random_Coefficient_System is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(p);
    declare
      s : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p.all);
      n : constant integer32 := p'last;
      q : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range)
        := Random_Coefficient_Systems.Create(natural32(n),s);
    begin
      put_line(q);
      Write_to_File(q);
    end;
  end DoblDobl_Random_Coefficient_System;

  procedure QuadDobl_Random_Coefficient_System is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(p);
    declare
      s : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p.all);
      n : constant integer32 := p'last;
      q : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range)
        := Random_Coefficient_Systems.Create(natural32(n),s);
    begin
      put_line(q);
      Write_to_File(q);
    end;
  end QuadDobl_Random_Coefficient_System;

  procedure Main is

    ans,dqd : character;

  begin
    new_line;
    put_line("Generation of random coefficient polynomial systems :");
    put_line("  1. create one random coefficient polynomial; or ");
    put_line("  2. generate a system of random coefficient polynomials;");
    put_line("  3. generate random coefficients for given system.");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    new_line;
    put_line("MENU for the precision : ");
    put_line("  0. standard double complex numbers;");
    put_line("  1. double double complex numbers;");
    put_line("  2. quad double complex numbers.");
    put("Type 0, 2, or 2 to make your choice : ");
    Ask_Alternative(dqd,"012");
    new_line;
    case ans is
      when '1' =>
        case dqd is
          when '0' => Standard_Creation_of_Random_Polynomials;
          when '1' => DoblDobl_Creation_of_Random_Polynomials;
          when '2' => QuadDobl_Creation_of_Random_Polynomials;
          when others => null;
        end case;
      when '2' => Random_Polynomial_System(dqd);
      when '3' =>
        case dqd is
          when '0' => Standard_Random_Coefficient_System;
          when '1' => DoblDobl_Random_Coefficient_System;
          when '2' => QuadDobl_Random_Coefficient_System;
          when others => null;
        end case;
      when others => null;
    end case;
  end Main;

end Test_RanCff_Systems;
