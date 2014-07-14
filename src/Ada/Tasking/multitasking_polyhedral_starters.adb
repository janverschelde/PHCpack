with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with Standard_Integer_Linear_Solvers;
with Standard_Integer64_Linear_Solvers;
with Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Standard_Tableau_Formats;
with Standard_Radial_Solvers;
with Standard_Binomial_Systems;
with Standard_Binomial_Solvers;
with DoblDobl_Tableau_Formats;
with DoblDobl_Radial_Solvers;
with DoblDobl_Binomial_Systems;
with DoblDobl_Binomial_Solvers;
with QuadDobl_Tableau_Formats;
with QuadDobl_Radial_Solvers;
with QuadDobl_Binomial_Systems;
with QuadDobl_Binomial_Solvers;
with Multitasking;
with Multitasking_Volume_Computation;    use Multitasking_Volume_Computation;
with Polyhedral_Start_Systems;           use Polyhedral_Start_Systems;

package body Multitasking_Polyhedral_Starters is

  procedure Silent_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector ) is

  -- NOTE :
  --   In an attempt to avoid threads racing for common data,
  --   three modfications (compared to the reporting versions):
  --   (1) every task has its own copy of the tableau data structure;
  --   (2) every task has a deep copy of the tableau data structure;
  --   (3) the mixed cells are distributed in advance of task launch.

    use Standard_Complex_Solutions;

    n : constant integer32 := q'last;
    cells,cells_last : array(1..nt) of Mixed_Subdivision;

    procedure Job ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   This code gets executed by task i (= task_id) out of a number of
    --   tasks equal to nb_tasks.  The workload assignment is static:
    --   as number k runs from 1 till the total number of cells,
    --   task i takes those cells with k mod nb_tasks = i - 1.

      cff : Standard_Complex_VecVecs.VecVec(q'range);
      exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
      mic : Mixed_Cell;
      s_c : Standard_Complex_Vectors.Vector(1..2*n); -- for fully mixed
      C : Standard_Complex_Matrices.Matrix(q'range,q'range); -- semi mixed
      piv : Standard_Integer_Vectors.Vector(q'range);
      info : integer32;
      A,M,U : Standard_Integer_Matrices.Matrix(q'range,q'range);
      b,wrk : Standard_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Standard_Floating_Vectors.Vector(b'range);
      bsc : Standard_Complex_Vectors.Vector(b'range);
      pdetU : natural32 := 0;
      mysols : Solution_List := sols(task_id);
      ptr : Solution_List;
      ls : Link_to_Solution;
      mycells : Mixed_Subdivision := cells(task_id);

    begin
      Standard_Tableau_Formats.Extract_Coefficients_and_Exponents_Copies
        (q,cff,exp);
      while not Is_Null(mycells) loop
        mic := Head_Of(mycells); 
        if r = n then
          Select_Coefficients(cff,exp,mic.pts.all,s_c);
          Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
        else
          Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
          Standard_Complex_Linear_Solvers.lufac(C,n,piv,info);
          Standard_Complex_Linear_Solvers.lusolve(C,n,piv,b);
        end if;
        U := A;
        Standard_Integer_Linear_Solvers.Upper_Triangulate(M,U);
        pdetU := Volume_of_Diagonal(U);
        brd := Standard_Radial_Solvers.Radii(b);
        bsc := Standard_Radial_Solvers.Scale(b,brd);
        Standard_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
        logbrd := Standard_Radial_Solvers.Log10(brd);
        logx := Standard_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
        logx := Standard_Radial_Solvers.Multiply(M,logx);
        e10x := Standard_Radial_Solvers.Exp10(logx);
        ptr := mysols;
        for i in 1..pdetU loop
          ls := Head_Of(ptr);
          Standard_Binomial_Systems.Eval(M,ls.v,wrk);
          ls.v := wrk;
          Standard_Radial_Solvers.Multiply(ls.v,e10x);
          ptr := Tail_Of(ptr);
        end loop;
        for j in 1..pdetU loop
          mysols := Tail_Of(mysols);
        end loop;
        mycells := Tail_Of(mycells);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   First the distribution of the mixed volume along the static 
    --   workload assignment is computed.
    --   Then the lists of start solutions can be allocated in advance,
    --   before the launching of the tasks.  
    --   Each task modifies its own list of solutions.
    --   At the end, the checking of start solution follows the same
    --   static workload assignment.

      vol : Standard_Integer_Vectors.Vector(1..nt);
      tmp : Mixed_Subdivision := mcc;
      mic : Mixed_Cell;
      ind : integer32;

    begin
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,mcc,vol,mixvol);
      for i in 1..nt loop
        sols(i) := Create(n,vol(i));
      end loop;
      for k in 1..Length_Of(mcc) loop
        mic := Head_Of(tmp);
        ind := 1 + (integer32(k) mod nt);
        Append(cells(ind),cells_last(ind),mic);
        tmp := Tail_Of(tmp);
      end loop;
      do_jobs(nt);
      if r = n
       then Check_Solutions(q,mcc,sols,res);
       else Check_Solutions(q,mix,mcc,sols,res);
      end if;
    end Main;

  begin
    Main;
  end Silent_Multithreaded_Solve_Start_Systems;

  procedure Reporting_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector ) is

    use Standard_Complex_Solutions;

    n : constant integer32 := q'last;
    vol : Standard_Integer_Vectors.Vector(1..nt);
    cff : Standard_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);

    procedure Job ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   This code gets executed by task i (= task_id) out of a number of 
    --   tasks equal to nb_tasks.  The workload assignment is static:
    --   as number k runs from 1 till the total number of cells,
    --   task i takes those cells with k mod nb_tasks = i - 1.

      tmp : Mixed_Subdivision := mcc;
      mic : Mixed_Cell;
      s_c : Standard_Complex_Vectors.Vector(1..2*n); -- for fully mixed
      C : Standard_Complex_Matrices.Matrix(q'range,q'range); -- semi mixed
      piv : Standard_Integer_Vectors.Vector(q'range);
      info : integer32;
      A,M,U : Standard_Integer_Matrices.Matrix(q'range,q'range);
      b,wrk : Standard_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Standard_Floating_Vectors.Vector(b'range);
      bsc : Standard_Complex_Vectors.Vector(b'range);
      pdetU,sum : natural32 := 0;
      mysols : Solution_List := sols(task_id);
      ptr : Solution_List;
      ls : Link_to_Solution;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      for k in 1..Length_Of(mcc) loop
        if integer32(k) mod nb_tasks = task_id-1 then
          put_line("task " & Multitasking.to_string(natural32(task_id))
                           & " is at cell "
                           & Multitasking.to_string(k));
          mic := Head_Of(tmp); 
          if r = n then
            Select_Coefficients(cff,exp,mic.pts.all,s_c);
            Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
          else
            Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
            Standard_Complex_Linear_Solvers.lufac(C,n,piv,info);
            Standard_Complex_Linear_Solvers.lusolve(C,n,piv,b);
          end if;
          U := A;
          Standard_Integer_Linear_Solvers.Upper_Triangulate(M,U);
          pdetU := Volume_of_Diagonal(U);
          brd := Standard_Radial_Solvers.Radii(b);
          bsc := Standard_Radial_Solvers.Scale(b,brd);
          Standard_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
          logbrd := Standard_Radial_Solvers.Log10(brd);
          logx := Standard_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
          logx := Standard_Radial_Solvers.Multiply(M,logx);
          e10x := Standard_Radial_Solvers.Exp10(logx);
          ptr := mysols;
          for i in 1..pdetU loop
            ls := Head_Of(ptr);
            Standard_Binomial_Systems.Eval(M,ls.v,wrk);
            ls.v := wrk;
            Standard_Radial_Solvers.Multiply(ls.v,e10x);
            ptr := Tail_Of(ptr);
          end loop;
          for j in 1..pdetU loop
            mysols := Tail_Of(mysols);
          end loop;
          sum := sum + pdetU;
          put_line("Task " & Multitasking.to_string(natural32(task_id))
                           & " computed "
                           & Multitasking.to_string(pdetU)
                           & " start solutions.");
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      put_line("Task " & Multitasking.to_string(natural32(task_id))
                       & " computed in total "
                       & Multitasking.to_string(sum)
                       & " out of "
                       & Multitasking.to_string(natural32(vol(task_id)))
                       & " start solutions.");

    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   First the distribution of the mixed volume along the static 
    --   workload assignment is computed.
    --   Then the lists of start solutions can be allocated in advance,
    --   before the launching of the tasks.  
    --   Each task modifies its own list of solutions.
    --   At the end, the checking of start solution follows the same
    --   static workload assignment.

    begin
      new_line;
      put_line("computing volumes ...");
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,mcc,vol,mixvol);
      new_line;
      put("volume distribution : "); put(vol); new_line;
      put("the mixed volume : "); put(mixvol,1); new_line;
      new_line;
      put_line("allocating solution lists ...");
      for i in 1..nt loop
        sols(i) := Create(n,vol(i));
      end loop;
      Standard_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
      new_line;
      put_line("launching tasks ...");
      new_line;
      do_jobs(nt);
      new_line;
      put_line("checking the solutions ...");
      new_line;
      if r = q'last
       then Check_Solutions(cff,exp,mcc,sols,res);
       else Check_Solutions(q,mix,mcc,sols,res);
      end if;
    end Main;

  begin
    Main;
  end Reporting_Multithreaded_Solve_Start_Systems;

  procedure Silent_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector ) is

  -- NOTE :
  --   In an attempt to avoid threads racing for common data,
  --   three modfications (compared to the reporting versions):
  --   (1) every task has its own copy of the tableau data structure;
  --   (2) every task has a deep copy of the tableau data structure;
  --   (3) the mixed cells are distributed in advance of task launch.

    use DoblDobl_Complex_Solutions;

    n : constant integer32 := q'last;
    cells,cells_last : array(1..nt) of Mixed_Subdivision;

    procedure Job ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   This code gets executed by task i (= task_id) out of a number of
    --   tasks equal to nb_tasks.  The workload assignment is static:
    --   as number k runs from 1 till the total number of cells,
    --   task i takes those cells with k mod nb_tasks = i - 1.

      cff : DoblDobl_Complex_VecVecs.VecVec(q'range);
      exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
      mic : Mixed_Cell;
      s_c : DoblDobl_Complex_Vectors.Vector(1..2*n); -- for fully mixed
      C : DoblDobl_Complex_Matrices.Matrix(q'range,q'range); -- semi mixed
      piv : Standard_Integer_Vectors.Vector(q'range);
      info : integer32;
      A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
      b,wrk : DoblDobl_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Double_Double_Vectors.Vector(b'range);
      bsc : DoblDobl_Complex_Vectors.Vector(b'range);
      pdetU : natural64 := 0;
      mysols : Solution_List := sols(task_id);
      ptr : Solution_List;
      ls : Link_to_Solution;
      mycells : Mixed_Subdivision := cells(task_id);

    begin
      DoblDobl_Tableau_Formats.Extract_Coefficients_and_Exponents_Copies
        (q,cff,exp);
      while not Is_Null(mycells) loop
        mic := Head_Of(mycells); 
        if r = n then
          Select_Coefficients(cff,exp,mic.pts.all,s_c);
          Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
        else
          Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
          DoblDobl_Complex_Linear_Solvers.lufac(C,n,piv,info);
          DoblDobl_Complex_Linear_Solvers.lusolve(C,n,piv,b);
        end if;
        U := A;
        Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
        pdetU := Volume_of_Diagonal(U);
        brd := DoblDobl_Radial_Solvers.Radii(b);
        bsc := DoblDobl_Radial_Solvers.Scale(b,brd);
        DoblDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
        logbrd := DoblDobl_Radial_Solvers.Log10(brd);
        logx := DoblDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
        logx := DoblDobl_Radial_Solvers.Multiply(M,logx);
        e10x := DoblDobl_Radial_Solvers.Exp10(logx);
        ptr := mysols;
        for i in 1..pdetU loop
          ls := Head_Of(ptr);
          DoblDobl_Binomial_Systems.Eval(M,ls.v,wrk);
          ls.v := wrk;
          DoblDobl_Radial_Solvers.Multiply(ls.v,e10x);
          ptr := Tail_Of(ptr);
        end loop;
        for j in 1..pdetU loop
          mysols := Tail_Of(mysols);
        end loop;
        mycells := Tail_Of(mycells);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   First the distribution of the mixed volume along the static 
    --   workload assignment is computed.
    --   Then the lists of start solutions can be allocated in advance,
    --   before the launching of the tasks.  
    --   Each task modifies its own list of solutions.
    --   At the end, the checking of start solution follows the same
    --   static workload assignment.

      vol : Standard_Integer_Vectors.Vector(1..nt);
      tmp : Mixed_Subdivision := mcc;
      mic : Mixed_Cell;
      ind : integer32;

    begin
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,mcc,vol,mixvol);
      for i in 1..nt loop
        sols(i) := Create(n,vol(i));
      end loop;
      for k in 1..Length_Of(mcc) loop
        mic := Head_Of(tmp);
        ind := 1 + (integer32(k) mod nt);
        Append(cells(ind),cells_last(ind),mic);
        tmp := Tail_Of(tmp);
      end loop;
      do_jobs(nt);
      if r = n
       then Check_Solutions(q,mcc,sols,res);
       else Check_Solutions(q,mix,mcc,sols,res);
      end if;
    end Main;

  begin
    Main;
  end Silent_Multithreaded_Solve_Start_Systems;

  procedure Reporting_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector ) is

    use DoblDobl_Complex_Solutions;

    n : constant integer32 := q'last;
    vol : Standard_Integer_Vectors.Vector(1..nt);
    cff : DoblDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);

    procedure Job ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   This code gets executed by task i (= task_id) out of a number of 
    --   tasks equal to nb_tasks.  The workload assignment is static:
    --   as number k runs from 1 till the total number of cells,
    --   task i takes those cells with k mod nb_tasks = i - 1.

      tmp : Mixed_Subdivision := mcc;
      mic : Mixed_Cell;
      s_c : DoblDobl_Complex_Vectors.Vector(1..2*n); -- for fully mixed
      C : DoblDobl_Complex_Matrices.Matrix(q'range,q'range); -- semi mixed
      piv : Standard_Integer_Vectors.Vector(q'range);
      info : integer32;
      A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
      b,wrk : DoblDobl_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Double_Double_Vectors.Vector(b'range);
      bsc : DoblDobl_Complex_Vectors.Vector(b'range);
      pdetU,sum : natural64 := 0;
      mysols : Solution_List := sols(task_id);
      ptr : Solution_List;
      ls : Link_to_Solution;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      for k in 1..Length_Of(mcc) loop
        if integer32(k) mod nb_tasks = task_id-1 then
          put_line("task " & Multitasking.to_string(natural32(task_id))
                           & " is at cell "
                           & Multitasking.to_string(k));
          mic := Head_Of(tmp); 
          if r = n then
            Select_Coefficients(cff,exp,mic.pts.all,s_c);
            Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
          else
            Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
            DoblDobl_Complex_Linear_Solvers.lufac(C,n,piv,info);
            DoblDobl_Complex_Linear_Solvers.lusolve(C,n,piv,b);
          end if;
          U := A;
          Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
          pdetU := Volume_of_Diagonal(U);
          brd := DoblDobl_Radial_Solvers.Radii(b);
          bsc := DoblDobl_Radial_Solvers.Scale(b,brd);
          DoblDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
          logbrd := DoblDobl_Radial_Solvers.Log10(brd);
          logx := DoblDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
          logx := DoblDobl_Radial_Solvers.Multiply(M,logx);
          e10x := DoblDobl_Radial_Solvers.Exp10(logx);
          ptr := mysols;
          for i in 1..pdetU loop
            ls := Head_Of(ptr);
            DoblDobl_Binomial_Systems.Eval(M,ls.v,wrk);
            ls.v := wrk;
            DoblDobl_Radial_Solvers.Multiply(ls.v,e10x);
            ptr := Tail_Of(ptr);
          end loop;
          for j in 1..pdetU loop
            mysols := Tail_Of(mysols);
          end loop;
          sum := sum + pdetU;
          put_line("Task " & Multitasking.to_string(natural32(task_id))
                           & " computed "
                           & Multitasking.to_string(natural32(pdetU))
                           & " start solutions.");
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      put_line("Task " & Multitasking.to_string(natural32(task_id))
                       & " computed in total "
                       & Multitasking.to_string(natural32(sum))
                       & " out of "
                       & Multitasking.to_string(natural32(vol(task_id)))
                       & " start solutions.");

    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   First the distribution of the mixed volume along the static 
    --   workload assignment is computed.
    --   Then the lists of start solutions can be allocated in advance,
    --   before the launching of the tasks.  
    --   Each task modifies its own list of solutions.
    --   At the end, the checking of start solution follows the same
    --   static workload assignment.

    begin
      new_line;
      put_line("computing volumes ...");
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,mcc,vol,mixvol);
      new_line;
      put("volume distribution : "); put(vol); new_line;
      put("the mixed volume : "); put(mixvol,1); new_line;
      new_line;
      put_line("allocating solution lists ...");
      for i in 1..nt loop
        sols(i) := Create(n,vol(i));
      end loop;
      DoblDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
      new_line;
      put_line("launching tasks ...");
      new_line;
      do_jobs(nt);
      new_line;
      put_line("checking the solutions ...");
      new_line;
      if r = q'last
       then Check_Solutions(cff,exp,mcc,sols,res);
       else Check_Solutions(q,mix,mcc,sols,res);
      end if;
    end Main;

  begin
    Main;
  end Reporting_Multithreaded_Solve_Start_Systems;

  procedure Silent_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector ) is

  -- NOTE :
  --   In an attempt to avoid threads racing for common data,
  --   three modfications (compared to the reporting versions):
  --   (1) every task has its own copy of the tableau data structure;
  --   (2) every task has a deep copy of the tableau data structure;
  --   (3) the mixed cells are distributed in advance of task launch.

    use QuadDobl_Complex_Solutions;

    n : constant integer32 := q'last;
    cells,cells_last : array(1..nt) of Mixed_Subdivision;

    procedure Job ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   This code gets executed by task i (= task_id) out of a number of
    --   tasks equal to nb_tasks.  The workload assignment is static:
    --   as number k runs from 1 till the total number of cells,
    --   task i takes those cells with k mod nb_tasks = i - 1.

      cff : QuadDobl_Complex_VecVecs.VecVec(q'range);
      exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
      mic : Mixed_Cell;
      s_c : QuadDobl_Complex_Vectors.Vector(1..2*n); -- for fully mixed
      C : QuadDobl_Complex_Matrices.Matrix(q'range,q'range); -- semi mixed
      piv : Standard_Integer_Vectors.Vector(q'range);
      info : integer32;
      A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
      b,wrk : QuadDobl_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Quad_Double_Vectors.Vector(b'range);
      bsc : QuadDobl_Complex_Vectors.Vector(b'range);
      pdetU : natural64 := 0;
      mysols : Solution_List := sols(task_id);
      ptr : Solution_List;
      ls : Link_to_Solution;
      mycells : Mixed_Subdivision := cells(task_id);

    begin
      QuadDobl_Tableau_Formats.Extract_Coefficients_and_Exponents_Copies
        (q,cff,exp);
      while not Is_Null(mycells) loop
        mic := Head_Of(mycells); 
        if r = n then
          Select_Coefficients(cff,exp,mic.pts.all,s_c);
          Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
        else
          Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
          QuadDobl_Complex_Linear_Solvers.lufac(C,n,piv,info);
          QuadDobl_Complex_Linear_Solvers.lusolve(C,n,piv,b);
        end if;
        U := A;
        Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
        pdetU := Volume_of_Diagonal(U);
        brd := QuadDobl_Radial_Solvers.Radii(b);
        bsc := QuadDobl_Radial_Solvers.Scale(b,brd);
        QuadDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
        logbrd := QuadDobl_Radial_Solvers.Log10(brd);
        logx := QuadDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
        logx := QuadDobl_Radial_Solvers.Multiply(M,logx);
        e10x := QuadDobl_Radial_Solvers.Exp10(logx);
        ptr := mysols;
        for i in 1..pdetU loop
          ls := Head_Of(ptr);
          QuadDobl_Binomial_Systems.Eval(M,ls.v,wrk);
          ls.v := wrk;
          QuadDobl_Radial_Solvers.Multiply(ls.v,e10x);
          ptr := Tail_Of(ptr);
        end loop;
        for j in 1..pdetU loop
          mysols := Tail_Of(mysols);
        end loop;
        mycells := Tail_Of(mycells);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   First the distribution of the mixed volume along the static 
    --   workload assignment is computed.
    --   Then the lists of start solutions can be allocated in advance,
    --   before the launching of the tasks.  
    --   Each task modifies its own list of solutions.
    --   At the end, the checking of start solution follows the same
    --   static workload assignment.

      vol : Standard_Integer_Vectors.Vector(1..nt);
      tmp : Mixed_Subdivision := mcc;
      mic : Mixed_Cell;
      ind : integer32;

    begin
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,mcc,vol,mixvol);
      for i in 1..nt loop
        sols(i) := Create(n,vol(i));
      end loop;
      for k in 1..Length_Of(mcc) loop
        mic := Head_Of(tmp);
        ind := 1 + (integer32(k) mod nt);
        Append(cells(ind),cells_last(ind),mic);
        tmp := Tail_Of(tmp);
      end loop;
      do_jobs(nt);
      if r = n
       then Check_Solutions(q,mcc,sols,res);
       else Check_Solutions(q,mix,mcc,sols,res);
      end if;
    end Main;

  begin
    Main;
  end Silent_Multithreaded_Solve_Start_Systems;

  procedure Reporting_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector ) is

    use QuadDobl_Complex_Solutions;

    n : constant integer32 := q'last;
    vol : Standard_Integer_Vectors.Vector(1..nt);
    cff : QuadDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);

    procedure Job ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   This code gets executed by task i (= task_id) out of a number of 
    --   tasks equal to nb_tasks.  The workload assignment is static:
    --   as number k runs from 1 till the total number of cells,
    --   task i takes those cells with k mod nb_tasks = i - 1.

      tmp : Mixed_Subdivision := mcc;
      mic : Mixed_Cell;
      s_c : QuadDobl_Complex_Vectors.Vector(1..2*n); -- for fully mixed
      C : QuadDobl_Complex_Matrices.Matrix(q'range,q'range); -- semi mixed
      piv : Standard_Integer_Vectors.Vector(q'range);
      info : integer32;
      A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
      b,wrk : QuadDobl_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Quad_Double_Vectors.Vector(b'range);
      bsc : QuadDobl_Complex_Vectors.Vector(b'range);
      pdetU,sum : natural64 := 0;
      mysols : Solution_List := sols(task_id);
      ptr : Solution_List;
      ls : Link_to_Solution;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      for k in 1..Length_Of(mcc) loop
        if integer32(k) mod nb_tasks = task_id-1 then
          put_line("task " & Multitasking.to_string(natural32(task_id))
                           & " is at cell "
                           & Multitasking.to_string(k));
          mic := Head_Of(tmp); 
          if r = n then
            Select_Coefficients(cff,exp,mic.pts.all,s_c);
            Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
          else
            Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
            QuadDobl_Complex_Linear_Solvers.lufac(C,n,piv,info);
            QuadDobl_Complex_Linear_Solvers.lusolve(C,n,piv,b);
          end if;
          U := A;
          Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
          pdetU := Volume_of_Diagonal(U);
          brd := QuadDobl_Radial_Solvers.Radii(b);
          bsc := QuadDobl_Radial_Solvers.Scale(b,brd);
          QuadDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
          logbrd := QuadDobl_Radial_Solvers.Log10(brd);
          logx := QuadDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
          logx := QuadDobl_Radial_Solvers.Multiply(M,logx);
          e10x := QuadDobl_Radial_Solvers.Exp10(logx);
          ptr := mysols;
          for i in 1..pdetU loop
            ls := Head_Of(ptr);
            QuadDobl_Binomial_Systems.Eval(M,ls.v,wrk);
            ls.v := wrk;
            QuadDobl_Radial_Solvers.Multiply(ls.v,e10x);
            ptr := Tail_Of(ptr);
          end loop;
          for j in 1..pdetU loop
            mysols := Tail_Of(mysols);
          end loop;
          sum := sum + pdetU;
          put_line("Task " & Multitasking.to_string(natural32(task_id))
                           & " computed "
                           & Multitasking.to_string(natural32(pdetU))
                           & " start solutions.");
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      put_line("Task " & Multitasking.to_string(natural32(task_id))
                       & " computed in total "
                       & Multitasking.to_string(natural32(sum))
                       & " out of "
                       & Multitasking.to_string(natural32(vol(task_id)))
                       & " start solutions.");

    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   First the distribution of the mixed volume along the static 
    --   workload assignment is computed.
    --   Then the lists of start solutions can be allocated in advance,
    --   before the launching of the tasks.  
    --   Each task modifies its own list of solutions.
    --   At the end, the checking of start solution follows the same
    --   static workload assignment.

    begin
      new_line;
      put_line("computing volumes ...");
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,mcc,vol,mixvol);
      new_line;
      put("volume distribution : "); put(vol); new_line;
      put("the mixed volume : "); put(mixvol,1); new_line;
      new_line;
      put_line("allocating solution lists ...");
      for i in 1..nt loop
        sols(i) := Create(n,vol(i));
      end loop;
      QuadDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
      new_line;
      put_line("launching tasks ...");
      new_line;
      do_jobs(nt);
      new_line;
      put_line("checking the solutions ...");
      new_line;
      if r = q'last
       then Check_Solutions(cff,exp,mcc,sols,res);
       else Check_Solutions(q,mix,mcc,sols,res);
      end if;
    end Main;

  begin
    Main;
  end Reporting_Multithreaded_Solve_Start_Systems;

end Multitasking_Polyhedral_Starters;
