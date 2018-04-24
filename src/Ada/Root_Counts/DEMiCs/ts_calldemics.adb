with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;
with Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with DEMiCs_Command_Line;                use DEMiCs_Command_Line;

procedure ts_calldemics is

-- DESCRIPTION :
--   Tests the basic file based command line interface to DEMiCs.

  function Apply_Lifting
             ( sup : Lists_of_Integer_Vectors.List;
               lif : Standard_Floating_Vectors.Vector )
             return Lists_of_Floating_Vectors.List is

  -- DESCRIPTION :
  --   Given in sup is the support of a polynomial
  --   and in lif are the corresponding lifting values.
  --   Returns the support in sup, with the lifting applied
  --   as defined by the lifting values in lif.

  -- REQUIRED : Length_Of(sup) = lif'last.

    res,res_last : Lists_of_Floating_Vectors.List;
    tmp : Lists_of_Integer_Vectors.List := sup;
    ilv : Standard_Integer_Vectors.Link_to_Vector;
    idx : integer32 := 0;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      idx := idx + 1;
      ilv := Lists_of_Integer_Vectors.Head_Of(tmp);
      declare
        flv : Standard_Floating_Vectors.Vector(ilv'first..ilv'last+1);
      begin
        for k in ilv'range loop
          flv(k) := double_float(ilv(k));
        end loop;
        flv(flv'last) := lif(idx);
        Lists_of_Floating_Vectors.Append(res,res_last,flv);
      end;
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Apply_Lifting;

  function Apply_Lifting
             ( sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lif : Standard_Floating_VecVecs.VecVec )
             return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

  -- DESCRIPTION :
  --   Given in sup are the supports of a polynomial system
  --   and in lif are the corresponding lifting values.
  --   Returns the supports in sup, with the lifting applied
  --   as defined by the lifting values in lif.

  -- REQUIRED :
  --   for k in lif'range: length_of(sup(k)) = lif(k)'last.

    res : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range);

  begin
    for k in sup'range loop
      res(k) := Apply_Lifting(sup(k),lif(k).all);
    end loop;
    return res;
  end Apply_Lifting;

  function Make_Mixed_Cell
             ( dim : integer32;
               lbl : Standard_Integer_Vectors.Vector;
               lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
             return Mixed_Cell is

  -- DESCRIPTION :
  --   Returns the mixed cell corresponding to the labeled points.

  -- ON ENTRY :
  --   dim     dimension before lifting;
  --   lbl     labels to the points that span the mixed cell;
  --   lifsup  lifted supports.

    res : Mixed_Cell;
    idxfirst,idxsecond : integer32;
    tmp,last : Lists_of_Floating_Vectors.List;
    idx,done : integer32;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    res.nor := new Standard_Floating_Vectors.Vector'(1..dim+1 => 0.0);
    res.pts := new Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..dim);
    res.sub := null;
    for k in 1..dim loop
      last := res.pts(k);
      idxfirst := lbl(2*(k-1)+1);
      idxsecond := lbl(2*(k-1)+2);
      idx := 0;
      done := 2;
      tmp := lifsup(k);
      while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
        idx := idx + 1;
        if idx = idxfirst then
          lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
          Lists_of_Floating_Vectors.Append(res.pts(k),last,lpt.all);
          done := done - 1;
        elsif idx = idxsecond then
          lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
          Lists_of_Floating_Vectors.Append(res.pts(k),last,lpt.all);
          done := done - 1;
        end if;
        exit when (done = 0);
        tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Make_Mixed_Cell;

  function Make_Mixed_Cells
             ( dim : integer32;
               labels : Lists_of_Integer_Vectors.List;
               lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
             return Mixed_Subdivision is

  -- DESCRIPTION :
  --   Returns the mixed cells corresponding to the labeled points.

  -- ON ENTRY :
  --   dim     dimension before lifting;
  --   labels  labels to the points that span the mixed cell;
  --   lifsup  lifted supports.

    res,res_last : Mixed_Subdivision;
    tmp : Lists_of_Integer_Vectors.List := labels;
    lbl : Standard_Integer_Vectors.Link_to_Vector;
    mic : Mixed_Cell;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      lbl := Lists_of_Integer_Vectors.Head_Of(tmp);
      mic := Make_Mixed_Cell(dim,lbl.all,lifsup);
      Append(res,res_last,mic);
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Make_Mixed_Cells;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   prepares the input for demics, calls demics,
  --   and then parses the output file.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    randnbr : constant string := Random_Name("");
    iptname : constant string := "/tmp/demics_input" & randnbr;
    optname : constant string := "/tmp/demics_output" & randnbr;
    ans : character;
    verbose : boolean;
    mv : natural32;

  begin
    new_line;
    put_line("Calling DEMiCs for mixed volume computation ...");
    new_line;
    get(lp);
    new_line;
    put_line("Writing input to " & iptname & ".");
    put_line("Writing output to " & optname & ".");
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    declare
      sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(lp'range);
      lif : Standard_Floating_VecVecs.VecVec(lp'range);
      cells : Lists_of_Integer_Vectors.List;
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range);
      mix : Standard_Integer_Vectors.Vector(sup'range) := (sup'range => 1);
      mcc : Mixed_Subdivision;
    begin
      Prepare_Input(lp.all,iptname,sup,verbose);
      Call_DEMiCs(iptname,optname,verbose=>verbose);
      Process_Output(lp'last,optname,mv,lif,cells,verbose);
      new_line;
      put_line("The lifting values for the supports : "); put(lif);
      new_line;
      put_line("The indices to the mixed cells :"); put(cells);
      new_line;
      put("The mixed volume : "); put(mv,1); new_line;
      lifsup := Apply_Lifting(sup,lif);
      new_line;
      put_line("The lifting applied to the supports :");
      Floating_Mixed_Subdivisions_io.put(lifsup);
      mcc := Make_Mixed_Cells(sup'last,cells,lifsup);
      new_line;
      put_line("The mixed-cell configuration :");
      Floating_Mixed_Subdivisions_io.put(natural32(sup'last),mix,mcc,mv);
      put("The mixed volume : "); put(mv,1); new_line;
      
    end;
  end Main;

begin
  Main;
end ts_calldemics;
