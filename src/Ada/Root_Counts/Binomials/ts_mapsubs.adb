with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Monomial_Maps;             use Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;          use Standard_Monomial_Maps_io;
with Standard_Monomial_Map_Substitutors; use Standard_Monomial_Map_Substitutors;

procedure ts_mapsubs is

-- DESCRIPTION :
--   Prompts the user for a general polynomial system and a binomial system
--   with its solution maps.  Substitutes the maps into the polynomial system.

  procedure Substitute_Map
              ( p : in Laur_Sys; map : in Monomial_Map ) is

  -- DESCRIPTION :
  --   Substitutes the map into the polynomial system p.

    s : Laur_Sys(p'range);

  begin
    put_line("The map : "); put(map);
    put("number of variables : "); put(map.n,1); new_line;
    put("dimension of the map : "); put(map.d,1); new_line;
    put_line("The configuration of tropisms : ");
    put(Tropism_Configuration(map));
    for i in p'range loop
      s(i) := Subs(p(i),map);
    end loop;
    put_line("After substitution : "); put(s);
    declare
      tol : constant double_float := 1.0E-8;
      f : constant Laur_Sys := Filter(s,tol);
      ans : character;
      file : file_type;
    begin
      put_line("After substitution and filtering : "); 
      put(f'last,1); new_line;
      put(f);
      new_line;
      put("Write filtered system to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of the file ...");
        Read_Name_and_Create_File(file);
        put(file,f'last);
        put(file," ");
        put(file,Number_of_Unknowns(f(f'first)),1);
        new_line(file);
        put(file,f);
        close(file);
      end if;
    end;
  end Substitute_Map;

  procedure Substitute_Maps
              ( p : in Laur_Sys; maps : in Monomial_Map_List ) is

  -- DESCRIPTION :
  --   Substitutes the maps into the polynomial system p.

    tmp : Monomial_Map_List := maps;
    map : Link_to_Monomial_Map;

  begin
    while not Is_Null(tmp) loop
      map := Head_Of(tmp);
      Substitute_Map(p,map.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Substitute_Maps;

  procedure Main is

    p,q : Link_to_Laur_Sys;
    maps : Monomial_Map_List;

  begin
    new_line;
    put_line("Substituting a monomial map into a system ...");
    new_line;
    put_line("Reading a polynomial system ..."); get(p);
    new_line;
    put_line("Reading system and its monomial maps ...");
    Read_System_and_Maps(q,maps);
    new_line;
    put("number of solution maps : "); put(Length_Of(maps),1); new_line;
    Substitute_Maps(p.all,maps);
  end Main;

begin
  Main;
end ts_mapsubs;
