with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with String_Splitters;                  use String_Splitters;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Integer_Matrices;
-- with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;       use Lists_of_Integer_Vectors_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Binomial_Varieties;       use Standard_Binomial_Varieties;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;         use Standard_Monomial_Maps_io;
with Standard_Monomial_Map_Solvers;     use Standard_Monomial_Map_Solvers;
with Standard_Monomial_Map_Circuits;    use Standard_Monomial_Map_Circuits;
with Standard_Monomial_Map_Ideals;      use Standard_Monomial_Map_Ideals;
with Standard_Monomial_Map_Filters;     use Standard_Monomial_Map_Filters;
with Standard_Permanent_Factors;
with Black_Box_Binomial_Solvers;        use Black_Box_Binomial_Solvers;
-- with DoblDobl_Monomial_Maps;
-- with QuadDobl_Monomial_Maps;

procedure ts_monmap is

-- DESCRIPTION :
--   Development and testing of monomial maps.

  procedure Write_to_File ( maps : in Monomial_Map_Array ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and writes the maps to file.

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of an output file to write the maps ...");
    Read_Name_and_Create_File(file);
    put(file,maps'last(1),1);
   -- put(file," ");
   -- put(file,Top_Dimension(maps),1);
    new_line(file);
    put(file,maps);
    close(file);
  end Write_to_File;

  procedure Write_to_File ( maps : in Array_of_Monomial_Map_Lists ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and writes the maps to file.

    file : file_type;
    len : constant Standard_Natural_Vectors.Vector := Lengths(maps);
    sum : constant natural32 := Standard_Natural_Vectors.Sum(len);
    tmp : Monomial_Map_List;
    link_to_map : Link_to_Monomial_Map;

  begin
    new_line;
    put_line("Reading the name of an output file to write the maps ...");
    Read_Name_and_Create_File(file);
    put(file,sum,1); new_line(file);
    for i in reverse maps'range loop
      if not Is_Null(maps(i)) then
        tmp := maps(i);
        while not Is_Null(tmp) loop
          link_to_map := Head_Of(tmp);
          put(file,link_to_map);
          tmp := Tail_Of(tmp);
        end loop;
      end if;
    end loop;
    close(file);
  end Write_to_File;

  procedure Solve ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Solves the binomial system defined by p
  --   and stores the solutions as monomial maps.

    d : integer32;
    M : Standard_Integer_Matrices.Link_to_Matrix;
    c : Solution_List;
    fail : boolean;
    ans : character;
    maps : Link_to_Monomial_Map_Array := Toric_Solve(p);
    list_of_maps : Monomial_Map_List;

  begin
    Black_Box_Solver(p,fail,d,M,c);
    if maps /= null then
      Insert_Parameter_Symbols(natural32(d));
      put_line("the solution maps : "); put(maps.all);
      new_line;
      put("Write solution maps to file (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Write_to_File(maps.all);
      end if;
      list_of_maps := Create(maps.all);
      Clear(maps);
    end if;
    put_line("the list of solution maps : "); put(list_of_maps);
  end Solve;

  procedure Show_Free_Maps ( maps : in Monomial_Map_List ) is

  -- DESCRIPTION :
  --   Runs through the maps and shows the ones that are free.
 
    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;
 
  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      if Is_Free(link_to_map.all)
       then put(link_to_map);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Show_Free_Maps;

  procedure Show_Free_Maps
              ( p : in Laur_Sys; maps : in Monomial_Map_List ) is

  -- DESCRIPTION :
  --   Runs through the maps and shows the ones that are free.
 
    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;
   -- cnt : natural := 0;
   -- ans : character;
 
  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
     -- cnt := cnt + 1;
     -- put("checking map "); put(cnt,1); put_line(" :");
     -- Show_Ideal(p,link_to_map.all);
      if Is_Free(link_to_map.all)
       then Show_Ideal(p,link_to_map.all);
           -- put_line("the map is free!");
      end if;
     -- put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
     -- exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end Show_Free_Maps;
 
  procedure Show_Free_Maps ( maps : in Array_of_Monomial_Map_Lists ) is

  -- DESCRIPTION :
  --   Runs through the maps and shows the ones that are free.

  begin
    for i in reverse maps'range loop
      put("the free maps at dimension "); put(i,1); put_line(" :");
      Show_Free_Maps(maps(i));
    end loop;
  end Show_Free_Maps;

  procedure Show_Free_Maps
              ( p : in Laur_Sys; maps : in Array_of_Monomial_Map_Lists ) is

  -- DESCRIPTION :
  --   Runs through the maps and shows the ones that are free.

  begin
    for i in reverse maps'range loop
      put("the free maps at dimension "); put(i,1); put_line(" :");
      Show_Free_Maps(p,maps(i));
    end loop;
  end Show_Free_Maps;

  procedure Show_Maps ( maps : in Array_of_Monomial_Map_Lists ) is

  -- DESCRIPTION :
  --   Shows the maps to the user.

    len : constant Standard_Natural_Vectors.Vector := Lengths(maps);
    ans : character;

  begin
    put("Lengths of components : "); put(len); new_line;
    for i in reverse maps'range loop
      if len(i) > 0 then
        put("Do you want to see the "); put(len(i),1);
        put(" maps of dimension ");
        put(i,1); put(" ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans = 'n');
        put(maps(i));
      end if;
    end loop;
  end Show_Maps;

  procedure Process_Solution_Maps ( maps : in Monomial_Map_List ) is

    td : constant integer32 := integer32(Top_Dimension(maps));
    components : Array_of_Monomial_Map_Lists(0..td)
               := Pure_Dimensional_Maps(maps);
    ans : character;

  begin
    put("The top dimension : "); put(td,1); new_line;
    if Is_Pure_Dimensional(maps)
     then put_line("The list of maps is pure dimensional.");
     else put_line("The list of maps is not pure dimensional.");
    end if;
    Show_Maps(components);
    Show_Free_Maps(components);
    put_line("filtering the free submaps out ...");
    Filter_Free_Submaps(components);
    put_line("after filtering :");
    Show_Maps(components);
    new_line;
    put("Write the solution maps to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Write_to_File(components);
    end if;
    Clear(components);
  end Process_Solution_Maps;

  procedure Read_Solution_Maps is

  -- DESCRIPTION :
  --   Prompts the user for a file name and reads solution maps from file.

    file : file_type;
    maps_list : Monomial_Map_List;
    maps_array : Link_to_Monomial_Map_Array;

  begin
    put_line("Reading a file name for the solution maps ...");
    Read_Name_and_Open_File(file);
   -- get(file,maps_array);
    get(file,maps_list);
    put_line("the list : "); put(maps_list);
    close(file);
    maps_array := new Monomial_Map_Array'(Create(maps_list));
   -- clear(maps_list);
    put("Read "); put(maps_array'last,1); put_line(" monomial maps : ");
    put(maps_array.all);
    Process_Solution_Maps(maps_list);
  end Read_Solution_Maps;

  procedure Check_Maps ( p : in Laur_Sys; maps : in Monomial_Map_List ) is

  -- DESCRIPTION :
  --   Checks for free maps to belong to toric maps.

    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      if Is_Free(link_to_map.all) then
        put(link_to_map.all);
        put("-> monomial ideal : ");
        Show_Ideal(p,link_to_map.all);
      else
        if Has_Zeroes(link_to_map.all) then
          cnt := cnt + 1;
          put("-> affine map "); put(cnt,1); put_line(" with zeroes :");
          put(link_to_map.all);
          put("-> I : "); Show_Ideal(p,link_to_map.all);
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put("Found "); put(cnt,1); put_line(" affine maps.");
  end Check_Maps;

  procedure Write_Lengths ( c : in Array_of_Monomial_Map_Lists ) is

    len : constant Standard_Natural_Vectors.Vector := Lengths(c);
    sum : constant natural32 := Standard_Natural_Vectors.Sum(len);

  begin
    put("#components per dimension :"); put(len);
    put(" sum : "); put(sum,1); new_line;
  end Write_Lengths;

  procedure Filter_Solution_Maps
              ( p : in Laur_Sys; maps : in Monomial_Map_List ) is

    td : constant integer32 := integer32(Top_Dimension(maps));
    components : Array_of_Monomial_Map_Lists(0..td)
               := Pure_Dimensional_Maps(maps);
    ans : character;

  begin
    put("The top dimension : "); put(td,1); new_line;
    if Is_Pure_Dimensional(maps)
     then put_line("The list of maps is pure dimensional.");
     else put_line("The list of maps is not pure dimensional.");
    end if;
   -- Show_Maps(components);
    put_line("Representations as ideals :"); Show_Ideals(p,components);
    put_line("The free maps :"); Show_Free_Maps(p,components);
    put_line("filtering the free submaps out ...");
    Filter_Free_Submaps(components);
    put_line("after filtering :"); Show_Ideals(p,components);
    Write_Lengths(components);
   -- Show_Maps(components);
    put_line("filtering free maps as submaps of affine maps ...");
    Filter_Free_of_Affine_Submaps(p,components);
    put_line("after filtering :"); Show_Ideals(p,components);
    Write_Lengths(components);
    new_line;
    put_line("filtering maps as affine submaps of affine maps ...");
    Filter_Affine_Submaps(p,components);
    put_line("after filtering :"); Show_Ideals(p,components);
    Write_Lengths(components);
    new_line;
    put("Write the solution maps to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Write_to_File(components);
    end if;
    Clear(components);
  end Filter_Solution_Maps;

  procedure Read_System_and_Maps is

  -- DESCRIPTION :
  --   Prompts the user for a file name and reads binomial system
  --   and solutions from file.

    file : file_type;
    lp : Link_to_Laur_Sys;
    found : boolean;
    maps : Monomial_Map_List;
    ans : character;

  begin
    put_line("Reading a file name for the system ...");
    Read_Name_and_Open_File(file);
    get(file,lp);
    put("number of equations : "); put(lp'last,1); new_line;
    put("number of variables : ");
    put(Number_of_Unknowns(lp(lp'first)),1); new_line;
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    if not found then
      put_line("Found no solutions on file.");
    else
      get(file,maps);
      new_line;
      put("Read "); put(Length_Of(maps),1); put_line(" solution maps.");
     -- Show_Ideals(lp.all,maps);
     -- Show_Free_Maps(lp.all,maps);
      put("Continue filtering the maps ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Filter_Solution_Maps(lp.all,maps);
      end if;
      put("Continue checking the maps ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Check_Maps(lp.all,maps);
      end if;
    end if;
  end Read_System_and_Maps;

  procedure Solve_and_Filter
              ( outfilename : in string; p : in Laur_Sys;
                append_sols : in boolean ) is

  -- DESCRIPTION :
  --   Generates the affine solution maps of p and filters the submaps.
  --   The solution maps are written to file, or appended to the file
  --   with given name is the append_sols is true.

    maps : Monomial_Map_List;
    components : Link_to_Array_of_Monomial_Map_Lists;
    fail : boolean;
    ans : character;
    use Standard_Permanent_Factors;

  begin
    Interactive_Affine_Solutions_with_Iterator(p,maps,fail);
    Reporting_Filter(standard_output,p,maps,components);
    new_line;
    put("Write the solution maps to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      if not append_sols then
        Write_to_File(components.all);
      else
        put_line("Appending solutions to the file " & outfilename);
        Append(outfilename,components.all);
      end if;
    end if;
  end Solve_and_Filter;

  function Product_of_Degrees ( p : Laur_Sys ) return natural32 is

  -- DESCRIPTION :
  --   Returns the product of the degrees of the polynomial in p.

    res : natural32 := 1;

  begin
    for i in p'range loop
      res := res*natural32(Degree(p(i)));
    end loop;
    return res;
  end Product_of_Degrees;

  procedure Call_Solve_and_Filter is

    file : file_type;
    lp : Link_to_Laur_Sys;
    name : Link_to_String;
    found,append_sols : boolean;

  begin
    put_line("Reading a binomial system ...");
    Read_Name_and_Open_File(file,name);
    get(file,lp);
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    append_sols := not found;
    close(file);
    new_line;
    put_line("solving the system ...");
    new_line;
    Solve_and_Filter(name.all,lp.all,append_sols);
  end Call_Solve_and_Filter;

  procedure Compute_Degrees_of_Maps is

    lp : Link_to_Laur_Sys;
    sols : Monomial_Map_List;

  begin
    Read_System_and_Maps(lp,sols);
    new_line;
    put("Product of the degrees in the system : ");
    put(Product_of_Degrees(lp.all),1); put_line(".");
    new_line;
    put("Read "); put(Length_Of(sols),1);
    put_line(" monomial maps from file.");
    new_line;
    Show_Degrees(lp.all,sols);
  end Compute_Degrees_of_Maps;

  procedure Compute_Circuit_Equations_for_Maps is

    lp : Link_to_Laur_Sys;
    sols,tmp : Monomial_Map_List;
    lmp : Link_to_Monomial_Map;
    dim : integer32;

  begin
    Read_System_and_Maps(lp,sols);
    new_line;
   -- put("Give the dimension of the circuits : "); get(dim);
    if not Is_Null(sols) then
      Remove_Parameter_Symbols(natural32(Head_Of(sols).n));
      tmp := sols;
      while not Is_Null(tmp) loop
        lmp := Head_Of(tmp);
        dim := lmp.d;
        declare
          c : List := Circuits(lmp.all,dim);
          e : Laur_Sys(1..integer32(Length_Of(c))) := Equations(lmp.all,c);
        begin
          put("The circuits of dimension "); put(dim,1); put_line(" :");
          put(c);
          put_line("The equations :"); put(e);
          Clear(c); Clear(e);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Compute_Circuit_Equations_for_Maps;

  procedure Write_Equations ( map : in Monomial_Map ) is

  -- DESCRIPTION :
  --   Writes the equations defined by the map to screen.
 
     eqs : constant Laur_Sys := Equations(map);

   begin
     put(eqs);
   end Write_Equations;

  procedure Compute_Equations_for_Maps is

    lp : Link_to_Laur_Sys;
    sols,tmp : Monomial_Map_List;

  begin
    Read_System_and_Maps(lp,sols);
    new_line;
    put("Product of the degrees in the system : ");
    put(Product_of_Degrees(lp.all),1); put_line(".");
    new_line;
    put("Read "); put(Length_Of(sols),1);
    put_line(" monomial maps from file.");
    new_line;
    if not Is_Null(sols) then
      Remove_Parameter_Symbols(natural32(Head_Of(sols).n));
      tmp := sols;
      for i in 1..Length_Of(sols) loop
        put("equations for map "); put(i,1); put_line(" :");
        Write_Equations(Head_Of(tmp).all);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Compute_Equations_for_Maps;

  procedure Call_Black_Box_Solver is

    file : file_type;
    lp : Link_to_Laur_Sys;
    name : Link_to_String;
    found,append_sols,fail : boolean;
    ans : character;
    sols : Link_to_Array_of_Monomial_Map_Lists;

  begin
    put_line("Reading a binomial system ...");
    Read_Name_and_Open_File(file,name);
    get(file,lp);
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    append_sols := not found;
    close(file);
    new_line;
    put("Do you want intermediate output ? ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Black_Box_Binomial_Solver(standard_output,lp.all,sols,fail);
     else Black_Box_Binomial_Solver(lp.all,false,sols,fail);
    end if;
    if not fail and sols /= null then
      if not append_sols then
        Write_to_File(sols.all);
      else
        put_line("Appending solutions to the file " & name.all);
        Append(name.all,sols.all);
      end if;
    end if;
  end Call_Black_Box_Solver;

  procedure Main is

    lp : Link_to_Laur_Sys;
    ans : character;

  begin
    new_line;
    put_line("Monomial maps represent solution sets of binomial systems.");
    new_line;
    put_line("MENU to test the operations on monomial maps :");
    put_line("  1. solve binomial system and write monomial maps;");
    put_line("  2. read monomial maps and process the maps;");
    put_line("  3. read binomial system and its monomial maps;");
    put_line("  4. solve with permanent and filter submaps;");
    put_line("  5. read system and its solution maps to compute the degree;");
    put_line("  6. read system and its maps to compute circuit equations;");
    put_line("  7. compute equations for solution maps of a binomial system;");
    put_line("  8. read system and call the black box solver.");
    put("Type 1, 2, 3, 4, 5, 6, 7, or 8 to select : ");
    Ask_Alternative(ans,"1234567");
    new_line;
    case ans is
      when '1' =>
        put_line("Reading a binomial system ...");
        get(lp);
        new_line;
        put_line("solving the system ...");
        new_line;
        Solve(lp.all);
      when '2' => Read_Solution_Maps;
      when '3' => Read_System_and_Maps;
      when '4' => Call_Solve_and_Filter;
      when '5' => Compute_Degrees_of_Maps;
      when '6' => Compute_Circuit_Equations_for_Maps;
      when '7' => Compute_Equations_for_Maps;
      when '8' => Call_Black_Box_Solver;
      when others => null;   
    end case;
  end Main;

begin
  Main;
end ts_monmap;
