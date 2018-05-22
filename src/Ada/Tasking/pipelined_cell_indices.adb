with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with String_Splitters;
with DEMiCs_Command_Line;
with DEMiCs_Algorithm;
with DEMiCs_Output_Data;
with Semaphore;
with Multitasking;

package body Pipelined_Cell_Indices is

  procedure Produce_Cells
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is
  begin
    DEMiCs_Output_Data.done := false;
    if verbose
     then put_line("Calling DEMiCs ...");
    end if;
    DEMiCs_Algorithm.Call_DEMiCs(mix,sup,verbose);
    DEMiCs_Output_Data.done := true;
    if verbose
     then DEMiCs_Algorithm.Show_Output;
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

    nbr : constant integer32
        := DEMiCs_Command_Line.Number_of_Points_in_Cell(mix.all);
    done : boolean := false;
    sem_data,sem_write : Semaphore.Lock;

    procedure Process_Cell ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Gets the string representation for the next cell indices
    --   and computes the integer indices to pass to the process.

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
        if verbose or process = null then
          Semaphore.Request(sem_write);
          put("Task "); put(i,1); put_line(" writes :");
          put_line(cell.all);
          put("Cell indices : "); put(inds); new_line;
          Semaphore.Release(sem_write);
        end if;
        if process /= null then
          process(i,mix,inds);
        end if;
      end loop;
    end Process_Cell;

    procedure Run_Tasks is new Multitasking.Silent_Workers(Process_Cell);

  begin
    DEMiCs_Output_Data.Initialize_Cell_Pointer;
    Run_Tasks(nt);
  end Consume_Cells;

end Pipelined_Cell_Indices;
