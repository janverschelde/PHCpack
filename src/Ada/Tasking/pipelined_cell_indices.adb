with Ada.Text_IO;                        use Ada.Text_IO;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
-- with String_Splitters;
with Floating_Mixed_Subdivisions_io;
-- with DEMiCs_Command_Line;
with DEMiCs_Output_Convertors;
-- with DEMiCs_Algorithm;
-- with DEMiCs_Output_Data;
with DEMiCs_Translated;
with DEMiCs_Output_Cells;
with Semaphore;
with Multitasking;

package body Pipelined_Cell_Indices is

  procedure Produce_Cells
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Pipelined_Cell_Indices.produce_cells 1 ...");
    end if;
   -- DEMiCs_Output_Data.done := false;
    DEMiCs_Output_Cells.done := false;
    if verbose
     then put_line("Calling DEMiCs ...");
    end if;
   -- DEMiCs_Algorithm.Call_DEMiCs(mix,sup,verbose);
   -- DEMiCs_Output_Data.done := true;
    if verbose 
     then DEMiCs_Translated.Call_DEMiCs(mix,sup,vrblvl=>99);
     else DEMiCs_Translated.Call_DEMiCs(mix,sup);
    end if;
    DEMiCs_Output_Cells.done := true;
    if verbose then
     -- DEMiCs_Algorithm.Show_Output;
      DEMiCs_Translated.Show_Output;
    end if;
  end Produce_Cells;

  procedure Produce_Cells
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stable : in boolean; stlb : in double_float;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Pipelined_Cell_Indices.produce_cells 2 ...");
    end if;
   -- DEMiCs_Output_Data.done := false;
    DEMiCs_Output_Cells.done := false;
    if verbose
     then put_line("Calling DEMiCs ...");
    end if;
   -- DEMiCs_Algorithm.Call_DEMiCs(mix,sup,stable,stlb,verbose);
   -- DEMiCs_Output_Data.done := true;
    if verbose
     then DEMiCs_Translated.Call_DEMiCs(mix,sup,stlb,vrblvl=>99);
     else DEMiCs_Translated.Call_DEMiCs(mix,sup,stlb);
    end if;
    DEMiCs_Output_Cells.done := true;
    if verbose then
     -- DEMiCs_Algorithm.Show_Output;
      DEMiCs_Translated.Show_Output;
    end if;
  end Produce_Cells;

  procedure Produce_Cells
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : in Standard_Floating_VecVecs.Link_to_VecVec;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is

   -- lifting : constant Standard_Floating_Vectors.Vector
   --         := DEMiCs_Algorithm.Flatten(lif);
    lifting : constant Standard_Floating_Vectors.Link_to_Vector
            := DEMiCs_Translated.Flatten(lif,vrblvl-1);

  begin
    if vrblvl > 0
     then put_line("-> in Pipelined_Cell_Indices.produce_cells 3 ...");
    end if;
   -- DEMiCs_Output_Data.done := false;
    DEMiCs_Output_Cells.done := false;
    if verbose
     then put_line("Calling DEMiCs ...");
    end if;
   -- DEMiCs_Algorithm.Call_DEMiCs(mix,sup,lifting,verbose);
    if verbose
     then DEMiCs_Translated.Call_DEMiCs(mix,sup,0.0,lifting,vrblvl=>99);
     else DEMiCs_Translated.Call_DEMiCs(mix,sup,0.0,lifting);
    end if;
   -- DEMiCs_Output_Data.done := true;
    DEMiCs_Output_Cells.done := true;
    if verbose then
     -- DEMiCs_Algorithm.Show_Output;
      DEMiCs_Translated.Show_Output;
    end if;
  end Produce_Cells;

  procedure Consume_Cells
              ( nt : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                process : access procedure
                  ( idtask : in integer32;
                    mix : in Standard_Integer_Vectors.Link_to_Vector;
                    idx : in Standard_Integer_Vectors.Vector ) := null;
                verbose : in boolean := true ) is

   -- nbr : constant integer32
   --     := DEMiCs_Command_Line.Number_of_Points_in_Cell(mix.all);
    done : boolean := false;
    sem_data,sem_write : Semaphore.Lock;

    procedure Process_Cell ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Gets the string representation for the next cell indices
    --   and computes the integer indices to pass to the process.

     -- inds : Standard_Integer_Vectors.Vector(1..nbr); 
     -- cell : String_Splitters.Link_to_String;
     -- use String_Splitters;
      cell : Standard_Integer_Vectors.Link_to_Vector;

      use Standard_Integer_Vectors;

    begin
      loop
        Semaphore.Request(sem_data);
       -- cell := DEMiCs_Output_Data.Get_Next_Cell_Indices;
        cell := DEMiCs_Output_Cells.Get_Next_Cell_Indices;
        done := (cell = null);
        Semaphore.Release(sem_data);
        exit when done;
        if verbose then
          Semaphore.Request(sem_write);
         -- put("Task "); put(i,1); put_line(" extract cell indices :");
         -- DEMiCs_Command_Line.Line2Cell_Indices(cell.all,nbr,mix,inds,true);
          put("Task "); put(i,1); put(" has cell "); put(cell); new_line;
          Semaphore.Release(sem_write);
       -- else
       --   DEMiCs_Command_Line.Line2Cell_Indices(cell.all,nbr,mix,inds,false);
       -- DEMiCs_Output_Cells has already the indices
        end if;
        if verbose or process = null then
          Semaphore.Request(sem_write);
         -- put("Task "); put(i,1); put_line(" writes :");
         -- put_line(cell.all);
         -- put("Cell indices : "); put(inds); new_line;
          put("Task "); put(i,1); put(" has cell "); put(cell); new_line;
          Semaphore.Release(sem_write);
        end if;
        if process /= null then
         -- process(i,mix,inds);
          process(i,mix,cell.all);
        end if;
      end loop;
    end Process_Cell;

    procedure Run_Tasks is new Multitasking.Silent_Workers(Process_Cell);

  begin
   -- DEMiCs_Output_Data.Initialize_Cell_Indices_Pointer;
    DEMiCs_Output_Cells.Initialize_Cell_Indices_Pointer;
    Run_Tasks(nt);
  end Consume_Cells;

  procedure Pipelined_Mixed_Indices
              ( nt : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                process : access procedure
                  ( idtask : in integer32;
                    mix : in Standard_Integer_Vectors.Link_to_Vector;
                    idx : in Standard_Integer_Vectors.Vector ) := null;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is

   -- nbr : constant integer32
   --     := DEMiCs_Command_Line.Number_of_Points_in_Cell(mix.all);
    sem_data,sem_write : Semaphore.Lock;

    procedure do_job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   The first task is the producer, the others are consumers.

     -- inds : Standard_Integer_Vectors.Vector(1..nbr); 
     -- cell : String_Splitters.Link_to_String;
     -- use String_Splitters;
      cell : Standard_Integer_Vectors.Link_to_Vector;

      use Standard_Integer_Vectors;

    begin
      if i = 1 then
        Produce_Cells(mix,sup,verbose,vrblvl-1);
      else
        loop
         -- put_line("Task " & Multitasking.to_string(i) & " entered loop.");
          Semaphore.Request(sem_data);
         -- put_line("Task " & Multitasking.to_string(i) & " inside lock.");
         -- cell := DEMiCs_Output_Data.Get_Next_Cell_Indices;
          cell := DEMiCs_Output_Cells.Get_Next_Cell_Indices;
          Semaphore.Release(sem_data);
         -- only done if empty cell and production stopped
         -- exit when (cell = null) and DEMiCs_Output_Data.done;
          exit when (cell = null) and DEMiCs_Output_Cells.done;
          if cell /= null then
           -- DEMiCs_Command_Line.Line2Cell_Indices
           --   (cell.all,nbr,mix,inds,false);
            if verbose or process = null then
              Semaphore.Request(sem_write);
             -- put("Task "); put(i,1); put_line(" writes :");
             -- put_line(cell.all);
             -- put("Cell indices : "); put(inds); new_line;
              put("Task "); put(i,1); put_line(" has cell ");
              put(cell); new_line;
              Semaphore.Release(sem_write);
            end if;
            if process /= null then
             -- process(i,mix,inds);
              process(i,mix,cell.all);
            end if;
          end if;
        end loop;
      end if;
    end do_job;

    procedure Run_Tasks is new Multitasking.Silent_Workers(do_job);
   -- procedure Run_Tasks is new Multitasking.Reporting_Workers(do_job);

  begin
   -- DEMiCs_Output_Data.Initialize_Cell_Indices_Pointer;
    DEMiCs_Output_Cells.Initialize_Cell_Indices_Pointer;
    Run_Tasks(nt);
  end Pipelined_Mixed_Indices;

  procedure Pipelined_Mixed_Cells
              ( nt,dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : in Standard_Floating_VecVecs.Link_to_VecVec;
                lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                process : access procedure
                  ( idtask : in integer32;
                    mix : in Standard_Integer_Vectors.Link_to_Vector;
                    mic : in Mixed_Cell ) := null;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is

   -- nbr : constant integer32
   --     := DEMiCs_Command_Line.Number_of_Points_in_Cell(mix.all);
    sem_data,sem_write : Semaphore.Lock;
    cells_processed : integer32 := 0;
    cells_computed : integer32 := -1;

    procedure do_job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   The first task is the producer, the others are consumers.

     -- inds : Standard_Integer_Vectors.Vector(1..nbr); 
     -- cell : String_Splitters.Link_to_String;
      cell : Standard_Integer_Vectors.Link_to_Vector;
      sub : Mixed_Subdivision;
      mic : Mixed_Cell;
      wrk : Standard_Floating_Vectors.Vector(1..dim+1);
      mat : Standard_Floating_Matrices.Matrix(1..dim,1..dim);
      rhs : Standard_Floating_Vectors.Vector(1..dim);
      ipvt : Standard_Integer_Vectors.Vector(1..dim);

     -- use String_Splitters;
      use Standard_Integer_Vectors;

    begin
      if vrblvl > 0 then
        put("-> in Pipelined_Cell_Indices.pipelined_mixed_cells.do_job,");
        put(" i : "); put(i,1); put_line(" ...");
      end if;
      if i = 1 then
        Produce_Cells(mix,sup,lif,verbose,vrblvl-1);
     -- cells_computed := integer32(DEMiCs_Output_Data.Number_of_Cell_Indices);
        cells_computed := integer32(DEMiCs_Output_Cells.Number_of_Cell_Indices);
      else
        loop
          Semaphore.Request(sem_data);
         -- cell := DEMiCs_Output_Data.Get_Next_Cell_Indices;
          cell := DEMiCs_Output_Cells.Get_Next_Cell_Indices;
          if cell /= null then
            cells_processed := cells_processed + 1;
           -- sub := DEMiCs_Output_Data.Get_Next_Allocated_Cell;
            sub := DEMiCs_Output_Cells.Get_Next_Allocated_Cell;
          end if;
          Semaphore.Release(sem_data);
          if cell /= null then
            mic := Head_Of(sub);
            if verbose then
              Semaphore.Request(sem_write);
             -- put("Task "); put(i,1); put_line(" extract cell indices :");
             -- DEMiCs_Command_Line.Line2Cell_Indices
             --   (cell.all,nbr,mix,inds,true);
              put("Task "); put(i,1); put(" has cell ");
              put(cell); new_line;
              Semaphore.Release(sem_write);
            -- else
            --   DEMiCs_Command_Line.Line2Cell_Indices
            --     (cell.all,nbr,mix,inds,false);
            end if;
           -- DEMiCs_Output_Convertors.Make_Mixed_Cell
           --   (mic,dim,mix.all,inds,lifsup,mat,rhs,ipvt,wrk,verbose);
            DEMiCs_Output_Convertors.Make_Mixed_Cell
              (mic,dim,mix.all,cell.all,lifsup,mat,rhs,ipvt,wrk,verbose);
            if verbose or process = null then
              Semaphore.Request(sem_write);
             -- put("Task "); put(i,1); put_line(" writes :");
             -- put_line(cell.all);
             -- put("Cell indices : "); put(inds); new_line;
              put("Task "); put(i,1); put(" has cell ");
              put(cell); new_line;
              Floating_Mixed_Subdivisions_io.put(natural32(dim),mix.all,mic);
              Semaphore.Release(sem_write);
            end if;
            if process /= null then
              process(i,mix,mic);
            end if;
          end if;
          exit when (cells_computed = cells_processed);
        end loop;
      end if;
    end do_job;

    procedure Run_Tasks is new Multitasking.Silent_Workers(do_job);
   -- procedure Run_Tasks is new Multitasking.Reporting_Workers(do_job);

  begin
    if vrblvl > 0 then
      put_line("-> in Pipelined_Cell_Indices.pipelined_mixed_cells ...");
    end if;
   -- DEMiCs_Output_Data.Initialize_Cell_Indices_Pointer;
    DEMiCs_Output_Cells.Initialize_Cell_Indices_Pointer;
    Run_Tasks(nt);
  end Pipelined_Mixed_Cells;

end Pipelined_Cell_Indices;
