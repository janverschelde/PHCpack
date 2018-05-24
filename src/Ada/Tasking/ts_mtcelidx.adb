with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with Lists_of_Strings;
with DEMiCs_Command_Line;
with DEMiCs_Output_Convertors;
with DEMiCs_Algorithm;                   use DEMiCs_Algorithm;
with DEMiCs_Output_Data;
with Pipelined_Cell_Indices;
with Semaphore;
with Multitasking;

procedure ts_mtcelidx is

-- DESCRIPTION :
--   Development of multitasked pipelined processing of the cell indices
--   produced by DEMiCs.

  procedure Write_Cell_Indices is

  -- DESCRIPTION :
  --   Writes the cell indices stored in DEMiCs_Output_Data.

    cell : String_Splitters.Link_to_String;

    use String_Splitters;

  begin
    DEMiCs_Output_Data.Initialize_Cell_Indices_Pointer;
    loop
      cell := DEMiCs_Output_Data.Get_Next_Cell_Indices;
      exit when (cell = null);
      put_line(cell.all);
    end loop;
  end Write_Cell_Indices;

  procedure Write_DEMiCs_Output
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

  -- DESCRIPTIN :
  --   Writes the cells computed by DEMiCs.

  -- ON ENTRY :
  --   dim      dimension of the points in each support;
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   verbose  flag for more information.

    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    mcc : Mixed_Subdivision;
    mv : natural32;

  begin
    Process_Output(dim,mix,sup,lifsup,mcc,verbose);
    put_line("The lifted supports :");
    Floating_Mixed_Subdivisions_io.put(lifsup);
    put_line("The mixed-cell configuration :");
    Floating_Mixed_Subdivisions_io.put(natural32(dim),mix.all,mcc,mv);
    put("The mixed volume : "); put(mv,1); new_line;
  end Write_DEMiCs_Output;

  function Random_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return Standard_Floating_VecVecs.Link_to_VecVec is

    res : Standard_Floating_VecVecs.Link_to_VecVec;
    resrep : Standard_Floating_VecVecs.VecVec(mix'range);
    idx : integer32 := 1;
    len : integer32;

  begin
    for i in resrep'range loop
      len := integer32(Lists_of_Integer_Vectors.Length_Of(sup(idx)));
      declare
        vals : Standard_Floating_Vectors.Vector(1..len);
      begin
        for j in 1..len loop
          vals(j) := Standard_Random_Numbers.Random;
        end loop;
        resrep(i) := new Standard_Floating_Vectors.Vector'(vals);
      end;
      idx := idx + mix(i);
    end loop;
    res := new Standard_Floating_VecVecs.VecVec'(resrep);
    put_line("The random lifting : ");
    Standard_Floating_VecVecs_io.put(res.all);
    return res;
  end Random_Lifting;

  procedure Compute_Mixed_Volume ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Computes the mixed volume of the Newton polytopes
  --   spanned by the supports of p.

    dim : constant integer32 := p'last;
    ans : character;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    verbose : boolean;
    nbtasks : integer32 := 0;
    cellcnt : natural32 := 0;
    sem_write : Semaphore.Lock;

    procedure Count_Cell
                ( i : in integer32;
                  m : in Standard_Integer_Vectors.Link_to_Vector;
                  indices : in Standard_Integer_Vectors.Vector ) is

    -- DESCRIPTION :
    --   Adds one to cellcnt each time when called.

    begin
      cellcnt := cellcnt + 1;
    end Count_Cell;

    procedure Write_Cell
                ( i : in integer32;
                  m : in Standard_Integer_Vectors.Link_to_Vector;
                  indices : in Standard_Integer_Vectors.Vector ) is

    -- DESCRIPTION :
    --   Task i writes the cell indices to screen.

    begin
      Semaphore.Request(sem_write);
      put("Task "); put(i,1);
      put(" writes cell indices : "); put(indices); new_line;
      cellcnt := cellcnt + 1;
      Semaphore.Release(sem_write);
    end Write_Cell;

    procedure Two_Stage_Test is

    -- DESCRIPTION :
    --   Runs the production before the processing
    --   as a 2-stage process without pipelining.

    begin
      put_line("The cell indices : ");
      if nbtasks <= 0 then
        Write_Cell_Indices;
      else
        Pipelined_Cell_Indices.Consume_Cells(nbtasks,mix);
        Pipelined_Cell_Indices.Consume_Cells
          (nbtasks,mix,Count_Cell'access,false);
        put("The number of cells : "); put(cellcnt,1); new_line;
      end if;
    end Two_Stage_Test;

    use Lists_of_Strings;

  begin
    new_line;
    put("Verbose ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Monitor the adding of cell indices ? (y/n) ");
    Ask_Yes_or_No(ans);
    DEMiCs_Output_Data.monitor := (ans = 'y');
    new_line;
    put("Give the number of tasks (0 for no tasks) : ");
    get(nbtasks);
    new_line;
    Extract_Supports(p,mix,sup,verbose);
    if nbtasks < 2 then
      Pipelined_Cell_Indices.Produce_Cells(mix,sup,verbose);
      if verbose
       then Write_DEMiCs_Output(dim,mix,sup);
      end if;
      Two_Stage_Test;
    else
      Pipelined_Cell_Indices.Pipelined_Mixed_Indices
        (nbtasks,mix,sup,Write_Cell'access,false);
      put("Number of cells written : "); put(cellcnt,1); new_line;
      put("Number of cells : ");
      put(Length_Of(DEMiCs_Output_Data.Retrieve_Cell_Indices),1); new_line;
    end if;
  end Compute_Mixed_Volume;

  procedure Run_Pipeline
              ( dim,nbtasks : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Computes the lifted supports and runs the pipeline.

    lif : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Random_Lifting(mix,sup);
    lsp : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
        := DEMiCs_Output_Convertors.Apply_Lifting(mix.all,sup,lif.all);
    cellcnt : natural32 := 0;
    sem_write : Semaphore.Lock;

    procedure Write_Cell
                ( i : in integer32;
                  m : in Standard_Integer_Vectors.Link_to_Vector;
                  mic : in Mixed_Cell ) is

    -- DESCRIPTION :
    --   Task i writes the cell to screen.

      sub : Mixed_Subdivision;

      use Standard_Floating_VecVecs;

    begin
      Semaphore.Request(sem_write);
      put("Task "); put(i,1); put_line(" writes cell : ");
      Floating_Mixed_Subdivisions_io.put(natural32(dim),mix.all,mic);
      cellcnt := cellcnt + 1;
      Semaphore.Release(sem_write);
    end Write_Cell;

    use Lists_of_Strings;

  begin
    DEMiCs_Output_Data.allocate := true;
    DEMiCs_Output_Data.Store_Dimension_and_Mixture(dim,mix);
    DEMiCs_Output_Data.Initialize_Allocated_Cell_Pointer;
    if nbtasks > 1 then
      Pipelined_Cell_Indices.Pipelined_Mixed_Cells
        (nbtasks,dim,mix,sup,lif,lsp,Write_Cell'access,verbose);
      put("Number of cells processed : "); put(cellcnt,1); new_line;
      put("Number of cells computed  : ");
      put(Length_Of(DEMiCs_Output_Data.Retrieve_Cell_Indices),1);
      new_line;
    end if;
  end Run_Pipeline;

  procedure Construct_Mixed_Cells ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Constructs the mixed cells in a regular subdivision of 
  --   the Newton polytopes spanned by the supports of p.

    dim : constant integer32 := p'last;
    ans : character;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    verbose : boolean;
    nbtasks : integer32 := 0;

  begin
    new_line;
    put("Verbose ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Monitor the adding of cell indices ? (y/n) ");
    Ask_Yes_or_No(ans);
    DEMiCs_Output_Data.monitor := (ans = 'y');
    new_line;
    put("Give the number of tasks (at least 2) : ");
    get(nbtasks);
    new_line;
    Extract_Supports(p,mix,sup,verbose);
    Run_Pipeline(dim,nbtasks,mix,sup,verbose);
  end Construct_Mixed_Cells;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and then prepares the input for DEMiCs.

    lp : Link_to_Poly_Sys;
    ans : character;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put_line("MENU to test pipelined computations :");
    put_line("  1. write cell indices with pipelining");
    put_line("  2. pipelined mixed cell construction");
    put("Type 1 or 2 to make a choice : ");
    Ask_Alternative(ans,"12");
    case ans is
      when '1' => Compute_Mixed_Volume(lp.all);
      when '2' => Construct_Mixed_Cells(lp.all);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtcelidx;
