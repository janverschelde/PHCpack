with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
-- with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Lists_of_Integer_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Floating_Lifting_Functions;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with Drivers_for_Static_Lifting;
-- with Lists_of_Strings;
with DEMiCs_Output_Convertors;
-- with DEMiCs_Algorithm;
-- with DEMiCs_Output_Data;
with DEMiCs_Translated;
with DEMiCs_Output_Cells;
with Pipelined_Cell_Indices;
with Semaphore;
with Pipelined_Polyhedral_Homotopies;   use Pipelined_Polyhedral_Homotopies;
with use_c2phc; -- must be added for the linking to succeed ...

package body Test_Multitasked_Pipelining is

-- DESCRIPTION :
--   Development of multitasked pipelined processing of the cell indices
--   produced by DEMiCs.

  procedure Write_Cell_Indices is

  -- DESCRIPTION :
  --   Writes the cell indices stored in DEMiCs_Output_Data,
  --   or with DEMiCs_Translated, in DEMiCs_Output_Cells.

   -- cell : String_Splitters.Link_to_String;
   -- use String_Splitters;

    cell : Standard_Integer_Vectors.Link_to_Vector;
    use Standard_Integer_Vectors;

  begin
   -- DEMiCs_Output_Data.Initialize_Cell_Indices_Pointer;
    DEMiCs_Output_Cells.Initialize_Cell_Indices_Pointer;
    loop
     -- cell := DEMiCs_Output_Data.Get_Next_Cell_Indices;
      cell := DEMiCs_Output_Cells.Get_Next_Cell_Indices;
      exit when (cell = null);
     -- put_line(cell.all);
      put(cell); new_line;
    end loop;
  end Write_Cell_Indices;

  procedure Write_DEMiCs_Output
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stable : in boolean; stlb : in double_float;
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
    mv,smv,tmv : natural32;

  begin
   -- DEMiCs_Algorithm.Process_Output(dim,mix,sup,lifsup,mcc,verbose);
    if verbose
     then DEMiCs_Translated.Process_Output(dim,mix,sup,lifsup,mcc,1);
     else DEMiCs_Translated.Process_Output(dim,mix,sup,lifsup,mcc);
    end if;
    put_line("The lifted supports :");
    Floating_Mixed_Subdivisions_io.put(lifsup);
    if not stable then
      put_line("The mixed-cell configuration :");
      Floating_Mixed_Subdivisions_io.put(natural32(dim),mix.all,mcc,mv);
      put("The mixed volume : "); put(mv,1); new_line;
    else
      Drivers_for_Static_Lifting.Floating_Volume_Computation
        (standard_output,dim,stlb,mix.all,mcc,mv,smv,tmv);
    end if;
  end Write_DEMiCs_Output;

  procedure Compute_Mixed_Volume ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Computes the mixed volume of the Newton polytopes
  --   spanned by the supports of p.

    dim : constant integer32 := p'last;
    ans : character;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    verbose,stable : boolean;
    genuine : constant boolean
            := Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(p);
    stlb : double_float;
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
      cellcnt := cellcnt + 1;
      put("Task "); put(i,1);
      put(" writes cell "); put(cellcnt,1);
      put(" : "); put(indices); new_line;
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

   -- use Lists_of_Strings;
    use Lists_of_Integer_Vectors;

  begin
    new_line;
    put("Verbose ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Monitor the adding of cell indices ? (y/n) ");
    Ask_Yes_or_No(ans);
   -- DEMiCs_Output_Data.monitor := (ans = 'y');
    DEMiCs_Output_Cells.monitor := (ans = 'y');
    if genuine then
      stable := false;
    else
      new_line;
      put("Compute the stable mixed volume ? (y/n) ");
      Ask_Yes_or_No(ans);
      stable := (ans = 'y');
    end if;
   -- DEMiCs_Output_Data.stable := stable;
    DEMiCs_Output_Cells.stable := stable;
    new_line;
    put("Give the number of tasks (0 for no tasks) : ");
    get(nbtasks);
    new_line;
   -- DEMiCs_Algorithm.Extract_Supports(p,mix,sup,verbose);
    DEMiCs_Translated.Extract_Supports(p,mix,sup,verbose);
    if stable
     then stlb := Floating_Lifting_Functions.Lifting_Bound(p);
     else stlb := 0.0;
    end if;
    if nbtasks < 2 then
      if stable
       then Pipelined_Cell_Indices.Produce_Cells(mix,sup,true,stlb,verbose);
       else Pipelined_Cell_Indices.Produce_Cells(mix,sup,verbose);
      end if;
      if verbose
       then Write_DEMiCs_Output(dim,mix,sup,stable,stlb);
      end if;
      Two_Stage_Test;
    else
      Pipelined_Cell_Indices.Pipelined_Mixed_Indices
        (nbtasks,mix,sup,Write_Cell'access,false);
      put("Number of cells written : "); put(cellcnt,1); new_line;
      put("Number of cells : ");
     -- put(Length_Of(DEMiCs_Output_Data.Retrieve_Cell_Indices),1); new_line;
      put(Length_Of(DEMiCs_Output_Cells.Retrieve_Cell_Indices),1); new_line;
    end if;
  end Compute_Mixed_Volume;

  procedure Test_Pipeline
              ( dim,nbtasks : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Computes the lifted supports and runs the pipeline as a test.

   -- lif : constant Standard_Floating_VecVecs.Link_to_VecVec
   --     := DEMiCs_Algorithm.Random_Lifting(mix,sup);
    lif : constant Standard_Floating_VecVecs.Link_to_VecVec
        := DEMiCs_Translated.Random_Lifting(mix,sup);
    lsp : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
        := DEMiCs_Output_Convertors.Apply_Lifting(mix.all,sup,lif.all);
    cellcnt : natural32 := 0;
    sem_write : Semaphore.Lock;

    procedure Write_Cell
                ( i : in integer32;
                  m : in Standard_Integer_Vectors.Link_to_Vector;
                  mic : in Mixed_Cell ) is

    -- DESCRIPTION :
    --   Task i writes the cell to screen.

    begin
      Semaphore.Request(sem_write);
      put("Task "); put(i,1); put_line(" writes cell : ");
      Floating_Mixed_Subdivisions_io.put(natural32(dim),mix.all,mic);
      cellcnt := cellcnt + 1;
      Semaphore.Release(sem_write);
    end Write_Cell;

   -- use Lists_of_Strings;
    use Lists_of_Integer_Vectors;

  begin
    if vrblvl > 0
     then put_line("-> in Test_Multitasked_Pipelining.test_pipeline ...");
    end if;
   -- DEMiCs_Output_Data.allocate := true;
   -- DEMiCs_Output_Data.Store_Dimension_and_Mixture(dim,mix);
   -- DEMiCs_Output_Data.Initialize_Allocated_Cell_Pointer;
    DEMiCs_Output_Cells.allocate := true;
    DEMiCs_Output_Cells.Store_Dimension_and_Mixture(dim,mix);
    DEMiCs_Output_Cells.Initialize_Allocated_Cell_Pointer;
    if nbtasks > 1 then
      Pipelined_Cell_Indices.Pipelined_Mixed_Cells
        (nbtasks,dim,mix,sup,lif,lsp,Write_Cell'access,verbose,vrblvl-1);
      put("Number of cells processed : "); put(cellcnt,1); new_line;
      put("Number of cells computed  : ");
     -- put(Length_Of(DEMiCs_Output_Data.Retrieve_Cell_Indices),1);
      put(Length_Of(DEMiCs_Output_Cells.Retrieve_Cell_Indices),1);
      new_line;
    end if;
  end Test_Pipeline;

  procedure Random_Coefficient_System
              ( dim,nbtasks : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Generates random lifting values, prompts the user for a file name,
  --   and then runs a 2-stage pipeline with number of tasks nbtasks
  --   to solve a random coefficient system with polyhedral homotopies.

   -- lif : constant Standard_Floating_VecVecs.Link_to_VecVec
   --     := DEMiCs_Algorithm.Random_Lifting(mix,sup);
    lif : constant Standard_Floating_VecVecs.Link_to_VecVec
        := DEMiCs_Translated.Random_Lifting(mix,sup);
    q : Standard_Complex_Laur_Systems.Laur_Sys(sup'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    file : file_type;

  begin
    put("Reading a file name ");
    put_line("to write the random coefficient system ...");
    skip_line; -- capture new line symbol
    Read_Name_and_Create_File(file);
    new_line;
    Pipeline_Cells_to_Paths(dim,nbtasks,mix,sup,lif,q,qsols,verbose);
    put_line(file,q);
    put_line(file,"THE SOLUTIONS :");
    put(file,Standard_Complex_Solutions.Length_Of(qsols),
        natural32(dim),qsols);
  end Random_Coefficient_System;

  procedure Construct_Mixed_Cells
              ( p : in Laur_Sys; randstart : in boolean;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Constructs the mixed cells in a regular subdivision of 
  --   the Newton polytopes spanned by the supports of p.
  --   If randstart, then a random coefficient system will be
  --   constructed and solved by pipelined polyhedral homotopies.

    dim : constant integer32 := p'last;
    ans : character;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    verbose : boolean;
    nbtasks : integer32 := 0;

  begin
    if vrblvl > 0 then
      put_line("-> in Test_Multitasked_Pipelining.construct_mixed_cells ...");
    end if;
    new_line;
    put("Verbose ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Monitor the adding of cell indices ? (y/n) ");
    Ask_Yes_or_No(ans);
   -- DEMiCs_Output_Data.monitor := (ans = 'y');
    DEMiCs_Output_Cells.monitor := (ans = 'y');
    new_line;
    put("Give the number of tasks (at least 2) : ");
    get(nbtasks);
    new_line;
   -- DEMiCs_Algorithm.Extract_Supports(p,mix,sup,verbose);
    DEMiCs_Translated.Extract_Supports(p,mix,sup,verbose);
    if not randstart
     then Test_Pipeline(dim,nbtasks,mix,sup,verbose,vrblvl-1);
     else Random_Coefficient_System(dim,nbtasks,mix,sup,verbose);
    end if;
  end Construct_Mixed_Cells;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for a polynomial system
  --   and then prepares the input for DEMiCs.

    lp : Link_to_Laur_Sys;
    ans,vrb : character;
    vrblvl : integer32 := 99;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put_line("MENU to test pipelined computations :");
    put_line("  1. write cell indices with pipelining");
    put_line("  2. pipelined mixed cell construction");
    put_line("  3. apply pipelined polyhedral homotopies");
    put("Type 1, 2, or 3 to make a choice : ");
    Ask_Alternative(ans,"123");
    new_line;
    put("Run with positive verbose level ? (y/n) ");
    Ask_Yes_or_No(vrb);
    if vrb = 'n'
     then vrblvl := 0;
    end if;
    case ans is
      when '1' => Compute_Mixed_Volume(lp.all);
      when '2' => Construct_Mixed_Cells(lp.all,false,vrblvl);
      when '3' => Construct_Mixed_Cells(lp.all,true,vrblvl);
      when others => null;
    end case;
  end Main;

end Test_Multitasked_Pipelining;
