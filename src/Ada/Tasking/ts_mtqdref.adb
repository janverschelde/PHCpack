with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Polynomial_Convertors;    use QuadDobl_Polynomial_Convertors;
with QuadDobl_Complex_Poly_SysFun;      use QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;    use QuadDobl_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with QuadDobl_Root_Refiners;            use QuadDobl_Root_Refiners;
with Multitasking,Semaphore;

procedure ts_mtqdref is

-- DESCRIPTION :
--   Test on multitasking a solution list of a polynomial system
--   with double double complex arithmetic.

  procedure Refine_Solution
               ( id,nb : in integer32; output : in boolean;
                 ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat ) is

  -- DESCRIPTION :
  --   Task with identification number id reports the receipt of
  --   solution with number nb, with data in ls.

    epsxa : constant double_float := 1.0E-60;
    epsfa : constant double_float := 1.0E-60;
    numit : natural32 := 0;
    maxit : constant natural32 := 5;
    fail : boolean;

  begin
    if output then
      put_line("Task " & Multitasking.to_String(id)
                       & " received solution " 
                       & Multitasking.to_String(nb) & ".");
    end if;
    Silent_Newton(f,jf,ls.all,epsxa,epsfa,numit,maxit,fail);
  end Refine_Solution;

  procedure Multitasking_Refinement
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List; 
                 n : in integer32; output : in boolean ) is

  -- DESCRIPTION :
  --   Given a polynomial system p with solutions in sols,
  --   n threads will be launched to multitask on the solution list.

    use QuadDobl_Complex_Solutions;

    f : Eval_Poly_Sys(p'range) := Create(p);
    jf : Jaco_Mat(p'range,p'range) := Create(p);
    ejf : Eval_Jaco_Mat(jf'range(1),jf'range(2)) := Create(jf);
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
        Refine_Solution(i,myjob,output,ls,f,ejf);
      end loop;
    end Next_Solution;
    procedure do_jobs is new Multitasking.Reporting_Workers(Next_Solution);

  begin
    do_jobs(n);
    Clear(f); Clear(jf); Clear(ejf);
  end Multitasking_Refinement;

  procedure Refine ( file : in file_type;
                     n : in integer32; output : in boolean;
                     p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Refines the solutions s of p with complex double double arithmetic,
  --   using n threads.

    dd_p : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range)
         := Standard_Poly_Sys_to_QuadDobl_Complex(p);
    dd_s : QuadDobl_Complex_Solutions.Solution_List
         := QuadDobl_Complex_Solutions.Create(s);
    timer : Timing_Widget;

  begin
    tstart(timer);
    Multitasking_Refinement(dd_p,dd_s,n,output);
    tstop(timer);
    QuadDobl_System_and_Solutions_io.put(file,dd_p,dd_s);
    new_line(file);
    print_times(file,timer,"multitasking refinement with double doubles");
    QuadDobl_Complex_Poly_Systems.Clear(dd_p);
  end Refine;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, solutions,
  --   and the number of threads.

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    s : Standard_Complex_Solutions.Solution_List;
    n : integer32 := 0;
    file : file_type;
    ans : character;

  begin
    new_line;
    put_line("Test on multitasking refinement ...");
    new_line;
    put_line("Reading a system and its solutions ...");
    Standard_System_and_Solutions_io.get(p,s);
    new_line;
    put("Read "); put(Standard_Complex_Solutions.Length_Of(s),1);
    put_line(" solutions.");
    new_line;
    put("Give the number of threads : "); get(n);
    new_line;
    skip_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Do you want to monitor the progress of the refinement ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Refine(file,n,true,p.all,s);
     else Refine(file,n,false,p.all,s);
    end if;
  end Main;

begin
  Main;
end ts_mtqdref;
