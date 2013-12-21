with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Witness_Sets_io;
with Hypersurface_Witness_Sets_io;      use Hypersurface_Witness_Sets_io;

procedure ts_hypwit is

-- DESCRIPTION :
--   Standalone testing routine to compute a witness set for one polynomial.

  procedure Call_Driver ( p : in Poly ) is

    file : file_type;
    ans : character;
    output,fail : boolean;
    eps : double_float := 1.0E-12;

  begin
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    loop
      put("Current tolerance : "); put(eps,3); new_line;
      put("Do you want to change it ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
      put("Give new tolerance : "); get(eps); new_line;
    end loop;
    new_line;
    put_line("calling the root finder, see output file ...");
    new_line;
    Call_Root_Finder(file,p,output,eps,fail);
  end Call_Driver;

  procedure Blackbox_Root_Finder ( p : in Poly ) is

    n : constant natural32 := Number_of_Unknowns(p);
    eps : constant double_float := 1.0E-12;
    fail : boolean;
    e : Link_to_Poly_Sys;
    esols : Solution_List;

  begin
    Silent_Root_Finder(p,eps,fail,e,esols);
    Witness_Sets_io.Add_Embed_Symbols(n-1);
    new_line;
    if fail
     then put_line("A failure occurred!");
     else put_line("No failure occurred.");
    end if;
    new_line;
    put_line("The embedded system : "); put(e.all);
    new_line;
    put_line("The witness points :");
    put(standard_output,Length_Of(esols),natural32(Head_Of(esols).n),esols);
  end Blackbox_Root_Finder;

  procedure Main is

    lp : Link_to_Poly_Sys;
    ans : character;

  begin
    new_line;
    put_line("Computing a witness set for one polynomial...");
    new_line;
    put_line("Reading one polynomial in several variables...");
    get(lp);
    new_line;
    put_line("MENU to test the drivers :");
    put_line("  1. fully interactive driver;");
    put_line("  2. second driver without interaction;");
    put_line("  3. a blackbox root finder.");
    put("Type 1, 2, or 3 to choose : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Call_Root_Finder(lp(lp'first));
      when '2' => Call_Driver(lp(lp'first));
      when others => Blackbox_Root_Finder(lp(lp'first));
    end case;
  end Main;

begin
  Main;
end ts_hypwit;
