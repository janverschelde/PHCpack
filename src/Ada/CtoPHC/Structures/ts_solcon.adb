with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;     use Multprec_Complex_Solutions_io;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Multprec_Solutions_Container;

procedure ts_solcon is

-- DESCRIPTION :
--   This procedure is an interactive test on the operations
--   in the solution container.

  procedure Standard_Test_Retrievals is

  -- DESRIPTION :
  --   Test the retrieval of the solutions in the container.

    use Standard_Complex_Solutions;

    fail : boolean;
    ind : natural32 := 0;
    ls : Link_to_Solution;

  begin
    new_line;
    put_line("Testing the retrieval of solutions...");
    new_line;
    put("Number of solutions in the container : ");
    put(Standard_Solutions_Container.Length,1); new_line;
    put("Dimension of solution vectors : ");
    put(Standard_Solutions_Container.Dimension,1); new_line;
    loop
      ind := ind + 1;
      Standard_Solutions_Container.Retrieve(ind,ls,fail);
      exit when fail;
      put("Solution "); put(ind,1); put_line(" :");
      put(ls.all); new_line;
    end loop;
  end Standard_Test_Retrievals;

  procedure DoblDobl_Test_Retrievals is

  -- DESRIPTION :
  --   Test the retrieval of the solutions in the container.

    use DoblDobl_Complex_Solutions;

    fail : boolean;
    ind : natural32 := 0;
    ls : Link_to_Solution;

  begin
    new_line;
    put_line("Testing the retrieval of solutions...");
    new_line;
    put("Number of solutions in the container : ");
    put(DoblDobl_Solutions_Container.Length,1); new_line;
    put("Dimension of solution vectors : ");
    put(DoblDobl_Solutions_Container.Dimension,1); new_line;
    loop
      ind := ind + 1;
      DoblDobl_Solutions_Container.Retrieve(ind,ls,fail);
      exit when fail;
      put("Solution "); put(ind,1); put_line(" :");
      put(ls.all); new_line;
    end loop;
  end DoblDobl_Test_Retrievals;

  procedure QuadDobl_Test_Retrievals is

  -- DESRIPTION :
  --   Test the retrieval of the solutions in the container.

    use QuadDobl_Complex_Solutions;

    fail : boolean;
    ind : natural32 := 0;
    ls : Link_to_Solution;

  begin
    new_line;
    put_line("Testing the retrieval of solutions...");
    new_line;
    put("Number of solutions in the container : ");
    put(QuadDobl_Solutions_Container.Length,1); new_line;
    put("Dimension of solution vectors : ");
    put(QuadDobl_Solutions_Container.Dimension,1); new_line;
    loop
      ind := ind + 1;
      QuadDobl_Solutions_Container.Retrieve(ind,ls,fail);
      exit when fail;
      put("Solution "); put(ind,1); put_line(" :");
      put(ls.all); new_line;
    end loop;
  end QuadDobl_Test_Retrievals;

  procedure Multprec_Test_Retrievals is

  -- DESRIPTION :
  --   Test the retrieval of the solutions in the container.

    use Multprec_Complex_Solutions;

    fail : boolean;
    ind : natural32 := 0;
    ls : Link_to_Solution;

  begin
    new_line;
    put_line("Testing the retrieval of solutions...");
    new_line;
    put("Number of solutions in the container : ");
    put(Multprec_Solutions_Container.Length,1); new_line;
    put("Dimension of solution vectors : ");
    put(Multprec_Solutions_Container.Dimension,1); new_line;
    loop
      ind := ind + 1;
      Multprec_Solutions_Container.Retrieve(ind,ls,fail);
      exit when fail;
      put("Solution "); put(ind,1); put_line(" :");
      put(ls.all); new_line;
    end loop;
  end Multprec_Test_Retrievals;

  procedure Standard_Test_Additions
              ( sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Tests the constructors of the container.

    use Standard_Complex_Solutions;

    tmp,retrieved : Solution_List;

  begin
    new_line;
    put_line("Testing the addition of solutions...");
    new_line;
    tmp := sols;
    while not Is_Null(tmp) loop
      Standard_Solutions_Container.Append(Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    retrieved := Standard_Solutions_Container.Retrieve;
    put_line("The retrieved solution list :");
    if not Is_Null(retrieved)
     then put(Standard_Output,
              Length_Of(retrieved),natural32(Head_Of(retrieved).n),retrieved);
    end if;
  end Standard_Test_Additions;

  procedure DoblDobl_Test_Additions
              ( sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Tests the constructors of the container.

    use DoblDobl_Complex_Solutions;

    tmp,retrieved : Solution_List;

  begin
    new_line;
    put_line("Testing the addition of solutions...");
    new_line;
    tmp := sols;
    while not Is_Null(tmp) loop
      DoblDobl_Solutions_Container.Append(Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    retrieved := DoblDobl_Solutions_Container.Retrieve;
    put_line("The retrieved solution list :");
    if not Is_Null(retrieved)
     then put(Standard_Output,
              Length_Of(retrieved),natural32(Head_Of(retrieved).n),retrieved);
    end if;
  end DoblDobl_Test_Additions;

  procedure QuadDobl_Test_Additions
              ( sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Tests the constructors of the container.

    use QuadDobl_Complex_Solutions;

    tmp,retrieved : Solution_List;

  begin
    new_line;
    put_line("Testing the addition of solutions...");
    new_line;
    tmp := sols;
    while not Is_Null(tmp) loop
      QuadDobl_Solutions_Container.Append(Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    retrieved := QuadDobl_Solutions_Container.Retrieve;
    put_line("The retrieved solution list :");
    if not Is_Null(retrieved)
     then put(Standard_Output,
              Length_Of(retrieved),natural32(Head_Of(retrieved).n),retrieved);
    end if;
  end QuadDobl_Test_Additions;

  procedure Multprec_Test_Additions
              ( sols : in Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Tests the constructors of the container.

    use Multprec_Complex_Solutions;

    tmp,retrieved : Solution_List;

  begin
    new_line;
    put_line("Testing the addition of solutions...");
    new_line;
    tmp := sols;
    while not Is_Null(tmp) loop
      Multprec_Solutions_Container.Append(Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    retrieved := Multprec_Solutions_Container.Retrieve;
    put_line("The retrieved solution list :");
    if not Is_Null(retrieved)
     then put(Standard_Output,
              Length_Of(retrieved),natural32(Head_Of(retrieved).n),retrieved);
    end if;
  end Multprec_Test_Additions;

  procedure Main is

    st_sols : Standard_Complex_Solutions.Solution_List;
    dd_sols : DoblDobl_Complex_Solutions.Solution_List;
    qd_sols : QuadDobl_Complex_Solutions.Solution_List;
    mp_sols : Multprec_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put_line("MENU to test the operations in the solutions container :");
    put_line("  1. test with standard double complex numbers;");
    put_line("  2. test with double double complex numbers;");
    put_line("  3. test with quad double complex numbers;");
    put_line("  4. test with multiprecision complex numbers.");
    put("Type 1, 2, 3, or 4 to select the precision : ");
    Ask_Alternative(ans,"1234");
    new_line;
    case ans is
      when '1' =>
        Read(st_sols);
        Standard_Solutions_Container.Initialize(st_sols);
        Standard_Test_Retrievals;
        Standard_Solutions_Container.Clear;
        Standard_Test_Additions(st_sols);
      when '2' =>
        Read(dd_sols);
        DoblDobl_Solutions_Container.Initialize(dd_sols);
        DoblDobl_Test_Retrievals;
        DoblDobl_Solutions_Container.Clear;
        DoblDobl_Test_Additions(dd_sols);
      when '3' =>
        Read(qd_sols);
        QuadDobl_Solutions_Container.Initialize(qd_sols);
        QuadDobl_Test_Retrievals;
        QuadDobl_Solutions_Container.Clear;
        QuadDobl_Test_Additions(qd_sols);
      when '4' =>
        Read(mp_sols);
        Multprec_Solutions_Container.Initialize(mp_sols);
        Multprec_Test_Retrievals;
        Multprec_Solutions_Container.Clear;
        Multprec_Test_Additions(mp_sols);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_solcon;
