with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Total_Degree_Start_Systems;
with Standard_Complex_Solutions_io;
with Standard_Fabry_on_Homotopy;
with DoblDobl_Fabry_on_Homotopy;
with TripDobl_Fabry_on_Homotopy;
with QuadDobl_Fabry_on_Homotopy;
with PentDobl_Fabry_on_Homotopy;
with OctoDobl_Fabry_on_Homotopy;
with DecaDobl_Fabry_on_Homotopy;
with HexaDobl_Fabry_on_Homotopy;

package body Newton_Fabry_on_Homotopy is

  function Prompt_for_Precision return character is
    
    res : character;
   
  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. triple double precision");
    put_line("  3. quad double precision");
    put_line("  4. penta double precision");
    put_line("  5. octo double precision");
    put_line("  6. deca double precision");
    put_line("  7. hexa double precision");
    put("Type 0, 1, 2, 3, 4, 5, 6, or 7 to select a precision : ");
    Ask_Alternative(res,"01234567");
    return res;
  end Prompt_for_Precision;

  procedure Run_Newton_Fabry
              ( nbtasks : in natural32; precision : in character;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_fabry_on_homotopy.Run_Newton_Fabry ...");
    end if;
    case precision is
      when '0' => Standard_Fabry_on_Homotopy.Main(nbtasks,vrblvl-1);
      when '1' => DoblDobl_Fabry_on_Homotopy.Main(nbtasks,vrblvl-1);
      when '2' => TripDobl_Fabry_on_Homotopy.Main(nbtasks,vrblvl-1);
      when '3' => QuadDobl_Fabry_on_Homotopy.Main(nbtasks,vrblvl-1);
      when '4' => PentDobl_Fabry_on_Homotopy.Main(nbtasks,vrblvl-1);
      when '5' => OctoDobl_Fabry_on_Homotopy.Main(nbtasks,vrblvl-1);
      when '6' => DecaDobl_Fabry_on_Homotopy.Main(nbtasks,vrblvl-1);
      when '7' => HexaDobl_Fabry_on_Homotopy.Main(nbtasks,vrblvl-1);
      when others => null;
    end case;
  end Run_Newton_Fabry;

  procedure Main is

    prc : constant character := Prompt_for_Precision;

  begin
    Run_Newton_Fabry(0,prc);
  end Main;

  procedure Generate_Homotopy
              ( file : in file_type;
                dim : in integer32; tvalue : in double_float ) is

    d : constant Standard_Natural_Vectors.Vector(1..dim) := (1..dim => 2);
    ones : constant Standard_Complex_Vectors.Vector(1..dim)
         := (1..dim => Standard_Complex_Numbers.Create(1.0));
    q : Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
      := Total_Degree_Start_Systems.Start_System(d,ones);
    qsol : Standard_Complex_Solutions.Solution(dim+1);
    sols : Standard_Complex_Solutions.Solution_List;
    s : constant Symbol_Table.Array_of_Symbols(1..dim)
      := Symbol_Table.Standard_Symbols(dim);

  begin
    qsol.t := Standard_Complex_Numbers.Create(0.0);
    qsol.m := 1;
    qsol.v := (1..dim+1 => Standard_Complex_Numbers.Create(1.0));
    qsol.v(dim+1) := Standard_Complex_Numbers.Create(0.0);
    qsol.err := 0.0;
    qsol.rco := 1.0;
    qsol.res := 0.0;
    Standard_Complex_Solutions.Add(sols,qsol);
    put(file,dim,1); put(file," "); put(file,dim+1,1); new_line(file);
    put(file," x1");
    for i in 2..dim loop                 -- write x1 + x2 + ...
      put(file," + x"); put(file,i,1);
    end loop;
    for i in 1..dim loop                 -- write - x1 - x2 - ...
      put(file," - x"); put(file,i,1);
    end loop;
    put_line(file," + ");                -- to initialize symbol table
    for i in 1..dim loop
      put(file," (1-t)*(t -");
      put(file,tvalue); put(file,")*(");
      Standard_Complex_Polynomials_io.put_terms(file,q(i));
      put_line(file,")"); 
      put(file," + t*(t -"); put(file,tvalue); put(file,")*(");
      put(file,"x"); put(file,i,1); put_line(file,"^2);");
    end loop;
    Symbol_Table.Init(s);
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add_String("t");
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    Standard_Complex_Solutions_io.put(file,1,natural32(dim)+1,sols);
    Standard_Complex_Poly_Systems.Clear(q);
  end Generate_Homotopy;

  procedure Test is

    dim : integer32 := 0;
    tvalue : double_float := 0.0;
    file : file_type;

  begin
    new_line;
    put_line("Generating a test homotopy ...");
    put("Give the dimension : "); get(dim);
    put("Give a nonzero critical value for t : "); get(tvalue);
    new_line;
    put_line("Reading the name of the output file ..."); skip_line;
    Read_Name_and_Create_File(file);
    Generate_Homotopy(file,dim,tvalue);
    close(file);
  end Test;

end Newton_Fabry_on_Homotopy;
