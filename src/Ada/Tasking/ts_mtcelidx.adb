with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with DEMiCs_Command_Line;
with DEMiCs_Algorithm;                   use DEMiCs_Algorithm;
with DEMiCs_Output_Data;
with Semaphore;
with Multitasking;

procedure ts_mtcelidx is

-- DESCRIPTION :
--   Development of multitasked pipelined processing of the cell indices
--   produced by DEMiCs.

  procedure Produce_Cells
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Calls DEMiCs to produce the cells.

  -- ON ENTRY :
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   verbose  flag for more information.

  begin
    Call_DEMiCs(mix,sup,verbose);
    Show_Output;
  end Produce_Cells;

  procedure Write_Cell_Indices is

  -- DESCRIPTION :
  --   Writes the cell indices stored in DEMiCs_Output_Data.

    cell : String_Splitters.Link_to_String;

    use String_Splitters;

  begin
    DEMiCs_Output_Data.Initialize_Cell_Pointer;
    loop
      cell := DEMiCs_Output_Data.Get_Next_Cell;
      exit when (cell = null);
      put_line(cell.all);
    end loop;
  end Write_Cell_Indices;

  procedure Multitasked_Write_Cells
               ( nt : in integer32;
                 mix : in Standard_Integer_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Writes the cell indices stored in DEMiCs_Output_Data,
  --   using nt tasks.  In mix is the type of mixture,
  --   needed to compute the cell indices.

    nbr : constant integer32
        := DEMiCs_Command_Line.Number_of_Points_in_Cell(mix.all);
    done : boolean := false;
    sem_data : Semaphore.Lock;
    sem_write : Semaphore.Lock;

    procedure Write_Cell ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Gets the next cell indices and write the indices to screen.

      inds : Standard_Integer_Vectors.Vector(1..nbr); 
      cell : String_Splitters.Link_to_String;
      use String_Splitters;

    begin
      loop
        Semaphore.Request(sem_data);
        cell := DEMiCs_Output_Data.Get_Next_Cell;
        done := (cell = null);
        Semaphore.Release(sem_data);
        exit when done;
        DEMiCs_Command_Line.Line2Cell_Indices(cell.all,nbr,mix,inds,false);
        Semaphore.Request(sem_write);
        put("Task "); put(i,1); put_line(" writes :");
        put_line(cell.all);
        put("Cell indices : "); put(inds); new_line;
        Semaphore.Release(sem_write);
      end loop;
    end Write_Cell;

    procedure Write_Cells is new Multitasking.Silent_Workers(Write_Cell);

  begin
    DEMiCs_Output_Data.Initialize_Cell_Pointer;
    Write_Cells(nt);
  end Multitasked_Write_Cells;

  procedure Process_Cells
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

  -- DESCRIPTIN :
  --   Processing of the cells computed by DEMiCs.

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
  end Process_Cells;

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
    Extract_Supports(p,mix,sup,verbose);
    Produce_Cells(mix,sup,verbose);
    put_line("The cell indices : ");
    if nbtasks <= 0
     then Write_Cell_Indices;
     else Multitasked_Write_Cells(nbtasks,mix);
    end if;
   -- Process_Cells(dim,mix,sup,verbose);
  end Compute_Mixed_Volume;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and then prepares the input for DEMiCs.

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    Compute_Mixed_Volume(lp.all);
  end Main;

begin
  Main;
end ts_mtcelidx;
