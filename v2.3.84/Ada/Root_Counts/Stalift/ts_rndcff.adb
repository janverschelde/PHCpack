with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Lists_of_Floating_Vectors;
with Lists_of_Floating_Vectors_io;       use Lists_of_Floating_Vectors_io;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Arrays_of_Floating_Vector_Lists;
with Arrays_of_Floating_Vector_Lists_io; use Arrays_of_Floating_Vector_Lists_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Integer_Convertors;        use Floating_Integer_Convertors;
with Supports_of_Polynomial_Systems;
with Random_Coefficient_Systems;

procedure ts_rndcff is

-- DESCRIPTION :
--   Test on the generation of polynomials and systems with given
--   supports and with random coefficients.

  procedure Write_to_File ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Gives the user the opportunity to write the system q to file,
  --   prompting for a file name.

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

  procedure Test_Creation_of_Random_Polynomials is

    n,k : natural32 := 0;
    support : Lists_of_Integer_Vectors.List;
    fltsupp : Lists_of_Floating_Vectors.List;
    p : Poly;

  begin
    put("Give the number of variables : "); get(n);
    put("Give the number of monomials : "); get(k);
    put("Give "); put(k,1); put(" natural vectors of length ");
    put(n,1); put_line(" : "); get(n,k,support);
    p := Random_Coefficient_Systems.Create(n,support);
    put("The random polynomial :"); put_line(p);
    fltsupp := Convert(support);
    put_line("The floating point support list : "); put(fltsupp);
    Clear(p);
    p := Random_Coefficient_Systems.Create(n,fltsupp);
    put("The random polynomial created from a floating support list :");
    put_line(p);
  end Test_Creation_of_Random_Polynomials;

  procedure Test_Creation_of_Random_Polynomial_System 
              ( n,r : integer32 ) is

    mix : Standard_Integer_Vectors.Vector(1..r);
    int_supports : Arrays_of_Integer_Vector_Lists.Array_of_Lists(1..r);
    flt_supports : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r);
    k : natural32 := 0;
    p : Poly_Sys(1..n);

  begin
    put("Give the type of mixture : "); get(mix);
    put("The type of mixture is "); put(mix); new_line;
    for i in 1..r loop
      put("Give the cardinality of support "); put(i,1); put(" : ");
      get(k);
      put("Give "); put(k,1); put(" natural vectors of length ");
      put(n,1); put_line(" : "); get(natural32(n),k,int_supports(i));
    end loop;
    put_line("The supports : "); put(int_supports);
    p := Random_Coefficient_Systems.Create(natural32(n),mix,int_supports);
    put_line("The random coefficient system : "); put_line(p);
    flt_supports := Convert(int_supports);
    put_line("The floating supports : "); put(flt_supports);
    p := Random_Coefficient_Systems.Create(natural32(n),mix,flt_supports);
    put_line("The random coefficient system : "); put_line(p);
    Write_to_File(p);
  end Test_Creation_of_Random_Polynomial_System;

  procedure Test_Creation_of_Random_Polynomial_Systems is

    n,r : integer32 := 0;

  begin
    put("Give the number of variables : "); get(n);
    put("Give the number of different supports : "); get(r);
    Test_Creation_of_Random_Polynomial_System(n,r);
  end Test_Creation_of_Random_Polynomial_Systems;

  procedure Test_Random_Coefficient_System is

    p : Link_to_Poly_Sys;

  begin
    get(p);
    declare
      s : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p.all);
      n : constant integer32 := p'last;
      q : constant Poly_Sys(p'range)
        := Random_Coefficient_Systems.Create(natural32(n),s);
    begin
      put_line(q);
      Write_to_File(q);
    end;
  end Test_Random_Coefficient_System;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Generation of random coefficient polynomial systems :");
    put_line("  1. create one random coefficient polynomial; or ");
    put_line("  2. generate a system of random coefficient polynomials;");
    put_line("  3. generate random coefficients for given system.");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Test_Creation_of_Random_Polynomials;
      when '2' => Test_Creation_of_Random_Polynomial_Systems;
      when '3' => Test_Random_Coefficient_System;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_rndcff;
