with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;  use Standard_System_and_Solutions_io;
with Multitasking,Semaphore;

procedure ts_mtsols is

-- DESCRIPTION :
--   Test on multitasking a solution list of a polynomial system.

  procedure Report_Solution
               ( id,nb : in integer32; ls : Link_to_Solution ) is

  -- DESCRIPTION :
  --   Task with identification number id reports the receipt of
  --   solution with number nb, with data in ls.

  begin
    put_line("Task " & Multitasking.to_String(id)
                     & " received solution " 
                     & Multitasking.to_String(nb) & ".");
    delay 1.0;
  end Report_Solution;

  procedure Multitask_on_Solutions
               ( p : in Poly_Sys; sols : in Solution_List;
                 n : in integer32 ) is

  -- DESCRIPTION :
  --   Given a polynomial system p with solutions in sols,
  --   n threads will be launched to multitask on the solution list.

    ptr : Solution_List;
    cnt : integer32 := 0;

    procedure Next_Solution ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   The n threads will run through the solution list,
    --   advancing the pointer ptr in a critical section,
    --   simultanuously with the counter.

    -- CORRESPONDENCE BETWEEN cnt AND ptr :
    --   cnt = 0                   <=> ptr not set to sols yet
    --   cnt in 1..Length_Of(sols) <=> ptr points to current solution
    --   cnt = Length_Of(sols) + 1 <=> Is_Null(ptr)

      s : Semaphore.Lock;
      myjob : integer32;
      myptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        Semaphore.Request(s);
        if cnt = 0 then
          cnt := 1;
          ptr := sols;
        else
          cnt := cnt + 1;
          if not Is_Null(ptr)
           then ptr := Tail_Of(ptr);
          end if;
        end if;
        myjob := cnt;
        myptr := ptr;
        Semaphore.Release(s);
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        Report_Solution(i,myjob,ls);
      end loop;
    end Next_Solution;
    procedure do_jobs is new Multitasking.Reporting_Workers(Next_Solution);

  begin
    do_jobs(n);
  end Multitask_on_Solutions;

  procedure Main is

    p : Link_to_Poly_Sys;
    s : Solution_List;
    n : integer32 := 0;

  begin
    new_line;
    put_line("Test on multitasking a solution list of a polynomial system...");
    get(p,s);
    new_line;
    put("Read "); put(Length_Of(s),1); put_line(" solutions.");
    new_line;
    put("Give the number of threads : "); get(n);
    Multitask_on_Solutions(p.all,s,n);
  end Main;

begin
  Main;
end ts_mtsols;
