with text_io;                           use text_io;
with DoblDobl_Root_Refiners;            use DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;            use QuadDobl_Root_Refiners;
with DoblDobl_Solutions_Queue;
with QuadDobl_Solutions_Queue;
with Multitasking;

package body Multitasked_DD_QD_Refiners is

  procedure Refine_Solution
               ( id,nbt,solno : in integer32; output : in boolean;
                 ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 epsxa,epsfa : in double_float; maxit : natural32 ) is

    numit : natural32 := 0;
    fail : boolean;

  begin
    if output then
      put_line("Task " & Multitasking.to_String(id) & " out of "
                       & Multitasking.to_String(nbt)
                       & " received solution " 
                       & Multitasking.to_String(solno) & ".");
    end if;
    Silent_Newton(f,jf,ls.all,epsxa,epsfa,numit,maxit,fail);
  end Refine_Solution;

  procedure Refine_Solution
               ( id,nbt,solno : in integer32; output : in boolean;
                 ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 epsxa,epsfa : in double_float; maxit : natural32 ) is

    numit : natural32 := 0;
    fail : boolean;

  begin
    if output then
      put_line("Task " & Multitasking.to_String(id) & " out of "
                       & Multitasking.to_String(nbt)
                       & " received solution " 
                       & Multitasking.to_String(solno) & ".");
    end if;
    Silent_Newton(f,jf,ls.all,epsxa,epsfa,numit,maxit,fail);
  end Refine_Solution;

  procedure Multitasking_Refinement
               ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List; 
                 n : in integer32; output : in boolean;
                 epsxa,epsfa : in double_float; maxit : in natural32 ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;

    f : Eval_Poly_Sys(p'range) := Create(p);
    jf : Jaco_Mat(p'range,p'range) := Create(p);
    ejf : Eval_Jaco_Mat(jf'range(1),jf'range(2)) := Create(jf);

    procedure Silent_Refiner ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   The n threads will run through the solution list,
    --   and apply the root refiner without reporting of output.

      myjob : integer32;
      myptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        myptr := DoblDobl_Solutions_Queue.Next;
        myjob := DoblDobl_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        Refine_Solution(i,n,myjob,output,ls,f,ejf,epsxa,epsfa,maxit);
      end loop;
    end Silent_Refiner;

    procedure Reporting_Refiner ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   The n threads will run through the solution list,
    --   and apply the root refiner with reporting of output.

      myjob : integer32;
      myptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        myptr := DoblDobl_Solutions_Queue.Next;
        myjob := DoblDobl_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        Refine_Solution(i,n,myjob,output,ls,f,ejf,epsxa,epsfa,maxit);
      end loop;
    end Reporting_Refiner;

    procedure Multitasked_Silent_Refiner is
      new Multitasking.Silent_Workers(Silent_Refiner);

    procedure Multitasked_Reporting_Refiner is
      new Multitasking.Reporting_Workers(Reporting_Refiner);

  begin
    DoblDobl_Solutions_Queue.Initialize(sols);
    if output
     then Multitasked_Reporting_Refiner(n);
     else Multitasked_Silent_Refiner(n);
    end if;
    Clear(f); Clear(jf); Clear(ejf);
  end Multitasking_Refinement;

  procedure Multitasking_Refinement
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List; 
                 n : in integer32; output : in boolean;
                 epsxa,epsfa : in double_float; maxit : in natural32 ) is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    f : Eval_Poly_Sys(p'range) := Create(p);
    jf : Jaco_Mat(p'range,p'range) := Create(p);
    ejf : Eval_Jaco_Mat(jf'range(1),jf'range(2)) := Create(jf);

    procedure Silent_Refiner ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   The n threads will run through the solution list,
    --   and apply the root refiner without reporting of output.

      myjob : integer32;
      myptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        myptr := QuadDobl_Solutions_Queue.Next;
        myjob := QuadDobl_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        Refine_Solution(i,n,myjob,output,ls,f,ejf,epsxa,epsfa,maxit);
      end loop;
    end Silent_Refiner;

    procedure Reporting_Refiner ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   The n threads will run through the solution list,
    --   and apply the root refiner with reporting of output.

      myjob : integer32;
      myptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        myptr := QuadDobl_Solutions_Queue.Next;
        myjob := QuadDobl_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        Refine_Solution(i,n,myjob,output,ls,f,ejf,epsxa,epsfa,maxit);
      end loop;
    end Reporting_Refiner;

    procedure Multitasked_Silent_Refiner is
      new Multitasking.Silent_Workers(Silent_Refiner);

    procedure Multitasked_Reporting_Refiner is
      new Multitasking.Reporting_Workers(Reporting_Refiner);

  begin
    QuadDobl_Solutions_Queue.Initialize(sols);
    if output
     then Multitasked_Reporting_Refiner(n);
     else Multitasked_Silent_Refiner(n);
    end if;
    Clear(f); Clear(jf); Clear(ejf);
  end Multitasking_Refinement;

end Multitasked_DD_QD_Refiners;
