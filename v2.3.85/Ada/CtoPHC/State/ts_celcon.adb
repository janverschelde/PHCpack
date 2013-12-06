with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Floating_Lifting_Utilities;         use Floating_Lifting_Utilities;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Cells_Container;

procedure ts_celcon is

-- DESCRIPTION :
--   This procedure is an interactive test on the operations
--   in the cells container.

  function Get_Cell ( mcc : Mixed_Subdivision; k : natural32 )
                    return Mixed_Cell is

  -- DESCRIPTION :
  --   Returns the k-th cells of the mixed-subdivision,
  --   assuming k >= 0 and k <= Length_Of(mcc).

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

  procedure Retrieve_and_Show_Cell
              ( k,len,n : in natural32;
                mix : in Standard_Integer_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Retrieves the k-th cell from the container, and if no failure,
  --   the the cell is shown to the user.

  -- ON INPUT :
  --   k        index to desired cell, must be >= 0;
  --   len      len = Cells_Container.Length;
  --   n        dimension of the lifted supports;
  --   mix      type of mixture needed for volume computation.

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
  end Retrieve_and_Show_Cell;

  procedure Test_Retrievals is

  -- DESRIPTION :
  --   Test the retrieval of the cells in the container.

    use Standard_Integer_Vectors;

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
        Retrieve_and_Show_Cell(k,len,n,mix);
      end loop;
    end if;
  end Test_Retrievals;

  procedure Read_and_Initialize is

  -- DESCRIPTION :
  --   Prompts the user to supply a file name for the cells,
  --   and initializes the cells container.

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
  end Read_and_Initialize;

  procedure Read_and_Construct is

  -- DESCRIPTION :
  --   Prompts the user to supply a file name for the cells,
  --   and lets the user pick the cells to include in the container.

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
  end Read_and_Construct;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing the operations in the cells container...");
    new_line;
    put_line("Choose one of the following : ");
    put_line("  1. initialize cells container and test retrievals;");
    put_line("  2. selectively append cells to the cells container.");
    put("Type 1 or 2 to make your choice : "); Ask_Alternative(ans,"12");
    new_line;
    if ans = '1'
     then Read_and_Initialize;
     else Read_and_Construct;
    end if;
    Test_Retrievals;
    Cells_Container.Clear;
  end Main;

begin
  Main;
end ts_celcon;
