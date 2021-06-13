with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Floating_Lifting_Utilities;
with Integer_Lifting_Utilities;
with Arrays_of_Floating_Vector_Lists;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Supports_of_Polynomial_Systems;
with Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Integer_Mixed_Subdivisions;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;
with Main_Lifting_Functions;
with Cells_Container;
with Integer_Cells_Container;

procedure ts_celcon is

-- DESCRIPTION :
--   This procedure is an interactive test on the operations
--   in the cells container.

  function Get_Cell ( mcc : Integer_Mixed_Subdivisions.Mixed_Subdivision;
                      k : natural32 )
                    return Integer_Mixed_Subdivisions.Mixed_Cell is

  -- DESCRIPTION :
  --   Returns the k-th cells of the mixed-subdivision,
  --   assuming k >= 0 and k <= Length_Of(mcc).

    use Integer_Mixed_Subdivisions;

    res : Mixed_Cell;
    tmp : Mixed_Subdivision := mcc;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = k
       then return Head_Of(tmp);
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return res;
  end Get_Cell;

  function Get_Cell ( mcc : Floating_Mixed_Subdivisions.Mixed_Subdivision;
                      k : natural32 )
                    return Floating_Mixed_Subdivisions.Mixed_Cell is

  -- DESCRIPTION :
  --   Returns the k-th cells of the mixed-subdivision,
  --   assuming k >= 0 and k <= Length_Of(mcc).

    use Floating_Mixed_Subdivisions;

    res : Mixed_Cell;
    tmp : Mixed_Subdivision := mcc;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = k
       then return Head_Of(tmp);
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return res;
  end Get_Cell;

  procedure Integer_Retrieve_and_Show_Cell
              ( k,len,n : in natural32;
                mix : in Standard_Integer_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Retrieves the k-th cell from the container for cells induced
  --   by an integer valued lifting function, and if no failure,
  --   the the cell is shown to the user.

  -- ON INPUT :
  --   k        index to desired cell, must be >= 0;
  --   len      len = Cells_Container.Length;
  --   n        dimension of the lifted supports;
  --   mix      type of mixture needed for volume computation.

    use Integer_Mixed_Subdivisions;

    mic : Mixed_Cell;
    fail : boolean;
    mv : natural32;
    cnt,lab,normal : Standard_Integer_Vectors.Link_to_Vector;

  begin
    Integer_Cells_Container.Retrieve(k,mic,fail);
    if fail then
      put("Retrieval of cell "); put(k,1); put(" failed, is ");
      put(k,1); put(" > "); put(len,1); put_line(" ?");
    else
      put("Retrieval of cell "); put(k,1); put_line(" succeeded.");
      put(n-1,mix.all,mic,mv);
      put(" mixed volume of the cell : "); put(mv,1); new_line;
      Integer_Cells_Container.Retrieve_Mixed_Cell(k,fail,cnt,lab,normal);
      put("Number of points in each support of the cell :");
      put(cnt); new_line;
      put("Labels of the points : "); put(lab); new_line;
      put_line("The inner normal to the cell : "); put_line(normal);
    end if;
  end Integer_Retrieve_and_Show_Cell;

  procedure Floating_Retrieve_and_Show_Cell
              ( k,len,n : in natural32;
                mix : in Standard_Integer_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Retrieves the k-th cell from the container for cells induced
  --   by a floating-point lifting function, and if no failure,
  --   the the cell is shown to the user.

  -- ON INPUT :
  --   k        index to desired cell, must be >= 0;
  --   len      len = Cells_Container.Length;
  --   n        dimension of the lifted supports;
  --   mix      type of mixture needed for volume computation.

    use Floating_Mixed_Subdivisions;

    mic : Mixed_Cell;
    fail : boolean;
    mv : natural32;
    cnt,lab : Standard_Integer_Vectors.Link_to_Vector;
    normal : Standard_Floating_Vectors.Link_to_Vector;

  begin
    Cells_Container.Retrieve(k,mic,fail);
    if fail then
      put("Retrieval of cell "); put(k,1); put(" failed, is ");
      put(k,1); put(" > "); put(len,1); put_line(" ?");
    else
      put("Retrieval of cell "); put(k,1); put_line(" succeeded.");
      put(n-1,mix.all,mic,mv);
      put(" mixed volume of the cell : "); put(mv,1); new_line;
      Cells_Container.Retrieve_Mixed_Cell(k,fail,cnt,lab,normal);
      put("Number of points in each support of the cell :");
      put(cnt); new_line;
      put("Labels of the points : "); put(lab); new_line;
      put_line("The inner normal to the cell : "); put_line(normal);
    end if;
  end Floating_Retrieve_and_Show_Cell;

  procedure Integer_Test_Retrievals is

  -- DESRIPTION :
  --   Test the retrieval of the cells in the container,
  --   for cells induced by an integer value lifting function.

    use Standard_Integer_Vectors;
    use Arrays_of_Integer_Vector_Lists;
    use Arrays_of_Integer_Vector_Lists_io;

    mix : Standard_Integer_Vectors.Link_to_Vector;
    len,n,k : natural32 := 0;
    ans : character;
    lif : Link_to_Array_of_Lists;

  begin
    new_line;
    put_line("Testing the retrieval of cells...");
    new_line;
    len := Integer_Cells_Container.Length;
    put("Number of cells in the container : "); put(len,1); new_line;
    n := Integer_Cells_Container.Dimension;
    put("Dimension of the points in the cells : "); put(n,1); new_line;
    mix := Integer_Cells_Container.Type_of_Mixture;
    if mix = null then
      put_line("No type of mixture returned, no retrieving.");
    else
      put("Type of mixture : "); put(mix); new_line;
      put("Do you want to see the lifted supports ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        lif := Integer_Cells_Container.Lifted_Supports;
        if lif /= null
         then put_line("The lifted supports :"); put(lif.all);
         else put_line("The lifted supports are empty...");
        end if;
      end if;
      loop
        put("Give number of a cell to retrieve (0 to exit) : ");
        get(k);
        exit when (k = 0);
        Integer_Retrieve_and_Show_Cell(k,len,n,mix);
      end loop;
    end if;
  end Integer_Test_Retrievals;

  procedure Floating_Test_Retrievals is

  -- DESRIPTION :
  --   Test the retrieval of the cells in the container,
  --   for cells induced by a floating-point lifting function.

    use Standard_Integer_Vectors;
    use Arrays_of_Floating_Vector_Lists;

    mix : Standard_Integer_Vectors.Link_to_Vector;
    len,n,k : natural32 := 0;
    ans : character;
    lif : Link_to_Array_of_Lists;

  begin
    new_line;
    put_line("Testing the retrieval of cells...");
    new_line;
    len := Cells_Container.Length;
    put("Number of cells in the container : "); put(len,1); new_line;
    n := Cells_Container.Dimension;
    put("Dimension of the points in the cells : "); put(n,1); new_line;
    mix := Cells_Container.Type_of_Mixture;
    if mix = null then
      put_line("No type of mixture returned, no retrieving.");
    else
      put("Type of mixture : "); put(mix); new_line;
      put("Do you want to see the lifted supports ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        lif := Cells_Container.Lifted_Supports;
        if lif /= null
         then put_line("The lifted supports :"); put(lif.all);
         else put_line("The lifted supports are empty...");
        end if;
      end if;
      loop
        put("Give number of a cell to retrieve (0 to exit) : ");
        get(k);
        exit when (k = 0);
        Floating_Retrieve_and_Show_Cell(k,len,n,mix);
      end loop;
    end if;
  end Floating_Test_Retrievals;

  procedure Floating_Read_and_Initialize is

  -- DESCRIPTION :
  --   Prompts the user to supply a file name for the cells,
  --   and initializes the cells container, for cells induced
  --   by a floating-point lifting function.

    use Arrays_of_Floating_Vector_Lists;
    use Floating_Lifting_Utilities;
    use Floating_Mixed_Subdivisions;

    file : file_type;
    n,m : natural32 := 0;
    r : integer32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    put_line("Reading the file name for a mixed-cell configuration.");
    Read_Name_and_Open_File(file);
    get(file,n,m,mix,mcc);
    r := mix'last;
    lif := new Array_of_Lists'(Lifted_Supports(r,mcc));
    Cells_Container.Initialize(mix,lif,mcc);
  end Floating_Read_and_Initialize;

  procedure Integer_Read_and_Initialize is

  -- DESCRIPTION :
  --   Prompts the user to supply a file name for the cells,
  --   and initializes the cells container, for cells induced
  --   by an integer valued lifting function.

    use Arrays_of_Integer_Vector_Lists;
    use Integer_Lifting_Utilities;
    use Integer_Mixed_Subdivisions;

    file : file_type;
    n,m : natural32 := 0;
    r : integer32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    put_line("Reading the file name for a mixed-cell configuration.");
    Read_Name_and_Open_File(file);
    get(file,n,m,mix,mcc);
    r := mix'last;
    lif := new Array_of_Lists'(Lifted_Supports(r,mcc));
    Integer_Cells_Container.Initialize(mix,lif,mcc);
  end Integer_Read_and_Initialize;

  procedure Floating_Read_and_Construct is

  -- DESCRIPTION :
  --   Prompts the user to supply a file name for the cells,
  --   for a regular mixed cell configuration induced by a floating lifting,
  --   and lets the user pick the cells to include in the container.

    use Arrays_of_Floating_Vector_Lists;
    use Floating_Lifting_Utilities;
    use Floating_Mixed_Subdivisions;

    file : file_type;
    len,k,n,m : natural32 := 0;
    r : integer32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    put_line("Reading the file name for a mixed-cell configuration.");
    Read_Name_and_Open_File(file);
    get(file,n,m,mix,mcc);
    r := mix'last;
    lif := new Array_of_Lists'(Lifted_Supports(r,mcc));
    Cells_Container.Initialize(mix);
    Cells_Container.Initialize(lif);
    len := Length_Of(mcc);
    loop
      put("Give number of cell to add to container (0 to exit) : ");
      get(k);
      exit when (k = 0);
      if k > len then 
        put(k,1); put(" > "); put(len,1); put_line(", the number of cells");
        put_line("Please try again...");
      else
        Cells_Container.Append(Get_Cell(mcc,k));
      end if;
    end loop;
  end Floating_Read_and_Construct;

  procedure Integer_Read_and_Construct is

  -- DESCRIPTION :
  --   Prompts the user to supply a file name for the cells,
  --   for a regular mixed cell configuration induced by a integer lifting,
  --   and lets the user pick the cells to include in the container.

    use Arrays_of_Integer_Vector_Lists;
    use Integer_Lifting_Utilities;
    use Integer_Mixed_Subdivisions;

    file : file_type;
    len,k,n,m : natural32 := 0;
    r : integer32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    put_line("Reading the file name for a mixed-cell configuration.");
    Read_Name_and_Open_File(file);
    get(file,n,m,mix,mcc);
    r := mix'last;
    lif := new Array_of_Lists'(Lifted_Supports(r,mcc));
    Integer_Cells_Container.Initialize(mix);
    Integer_Cells_Container.Initialize(lif);
    len := Length_Of(mcc);
    loop
      put("Give number of cell to add to container (0 to exit) : ");
      get(k);
      exit when (k = 0);
      if k > len then 
        put(k,1); put(" > "); put(len,1); put_line(", the number of cells");
        put_line("Please try again...");
      else
        Integer_Cells_Container.Append(Get_Cell(mcc,k));
      end if;
    end loop;
  end Integer_Read_and_Construct;

  procedure Lift_and_Prune
              ( p : in Poly_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
                mcc : out Integer_Mixed_Subdivisions.Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Extracts the supports and computes the mixture.

    use Arrays_of_Integer_Vector_Lists;

    sup : Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);
    perms : Standard_Integer_Vectors.Link_to_Vector;
    r : integer32;

    use Main_Lifting_Functions;

  begin
    Mixed_Volume_Computation.Compute_Mixture(sup,mix,perms);
    r := mix'last;
    put("Number of distinct supports : "); put(r,1); new_line;
    declare
      mixsup : constant Array_of_Lists(mix'range)
             := Mixed_Volume_Computation.Typed_Lists(mix.all,sup);
      lifsup : Array_of_Lists(mix'range);
      lifted : Link_to_Array_of_Lists;
    begin
      for k in mixsup'range loop
        lifsup(k) := Read_Integer_Lifting(mixsup(k));
      end loop;
      lifted := new Array_of_Lists'(lifsup);
      put_line("The lifted supports : ");
      Arrays_of_Integer_Vector_Lists_io.put(lifsup);
      Integer_Cells_Container.Initialize(mix);
      Integer_Cells_Container.Initialize(lifted);
      Integer_Cells_Container.Make_Subdivision;
      mcc := Integer_Cells_Container.Retrieve;
    end;
  end Lift_and_Prune;

  procedure Integer_Make_Subdivision is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   extracts the supports, computes the type of mixture,
  --   and asks the user for a lifting for each point.

    lp : Link_to_Poly_Sys;
    mcc : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    n,mv : natural32;

  begin
    put_line("Reading a polynomial system ...");
    get(lp);
    Lift_and_Prune(lp.all,mix,mcc);
    put_line("The mixed cell configuration :");
    n := natural32(lp'last);
    Integer_Mixed_Subdivisions_io.put(n,mix.all,mcc,mv);
    put("The mixed volume : "); put(mv,1); new_line;
  end Integer_Make_Subdivision;

  procedure Integer_Test is

  -- DESCRIPTION :
  --   Test on cell configurations defined by integer lifting functions.

    ans : character;

  begin
    new_line;
    put_line("Choose one of the following : ");
    put_line("  1. initialize cells container and test retrievals;");
    put_line("  2. selectively append cells to the cells container;");
    put_line("  3. give polynomial system and make a subdivision.");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Integer_Read_and_Initialize;
      when '2' => Integer_Read_and_Construct;
      when '3' => Integer_Make_Subdivision;
      when others => null;
    end case;
    Integer_Test_Retrievals;
    Integer_Cells_Container.Clear;
  end Integer_Test;

  procedure Floating_Test is

  -- DESCRIPTION :
  --   Test on cell configurations defined by floating lifting functions.

    ans : character;
  
  begin
    new_line;
    put_line("Choose one of the following : ");
    put_line("  1. initialize cells container and test retrievals;");
    put_line("  2. selectively append cells to the cells container.");
    put("Type 1 or 2 to make your choice : "); Ask_Alternative(ans,"12");
    new_line;
    if ans = '1'
     then Floating_Read_and_Initialize;
     else Floating_Read_and_Construct;
    end if;
    Floating_Test_Retrievals;
    Cells_Container.Clear;
  end Floating_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing the operations in the cells container...");
    new_line;
    put("Integer valued lifting function ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Integer_Test;
     else Floating_Test;
    end if;
  end Main;

begin
  Main;
end ts_celcon;
