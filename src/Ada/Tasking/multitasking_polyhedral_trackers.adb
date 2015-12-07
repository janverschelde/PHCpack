with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
--with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
--with Double_Double_Numbers;              use Double_Double_Numbers;
--with Quad_Double_Numbers;                use Quad_Double_Numbers;
--with Standard_Complex_Numbers;
--with DoblDobl_Complex_Numbers;
--with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
--with Standard_Complex_Norms_Equals;
--with DoblDobl_Complex_Vector_Norms;
--with QuadDobl_Complex_Vector_Norms;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Integer_Linear_Solvers;
with Standard_Integer64_Linear_Solvers;
with Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_Functions;
with Supports_of_Polynomial_Systems;
with Floating_Mixed_Subdivisions_io;
--with Mixed_Volume_Computation;
with Standard_Radial_Solvers;
with Standard_Binomial_Systems;
with Standard_Binomial_Solvers;
with DoblDobl_Radial_Solvers;
with DoblDobl_Binomial_Systems;
with DoblDobl_Binomial_Solvers;
with QuadDobl_Radial_Solvers;
with QuadDobl_Binomial_Systems;
with QuadDobl_Binomial_Solvers;
with Floating_Integer_Convertors;
with Floating_Lifting_Utilities;
--with Polyhedral_Coefficient_Homotopies;
with Multitasking,Semaphore;
with Mixed_Cells_Queue;
with Multitasking_Volume_Computation;    use Multitasking_Volume_Computation;
with Polyhedral_Start_Systems;           use Polyhedral_Start_Systems;
with Single_Polyhedral_Trackers;         use Single_Polyhedral_Trackers;
-- added for Check_Solution :
--with Strings_and_Numbers;                use Strings_and_Numbers;
--with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
--with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
--with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;

package body Multitasking_Polyhedral_Trackers is

  procedure Silent_Multitasking_Path_Tracker
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                h : in Standard_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                c : in Standard_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in Standard_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                mf : in Standard_Complex_Laur_JacoMats.Mult_Factors;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    cell_ptr : Mixed_Subdivision := mcc;
    s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : Standard_Complex_VecVecs.Array_of_VecVecs(1..nt);

    procedure Next_Cell ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Task i processes the next cell.

      mycell_ptr : Mixed_Subdivision;
      mic : Mixed_Cell;
      s_c : Standard_Complex_Vectors.Vector(1..2*n);    -- for fully mixed
      CC : Standard_Complex_Matrices.Matrix(1..n,1..n); -- for semi mixed
      piv : Standard_Integer_Vectors.Vector(1..n);
      info : integer32;
      A,M,U : Standard_Integer_Matrices.Matrix(q'range,q'range);
      b,wrk : Standard_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Standard_Floating_Vectors.Vector(b'range);
      bsc : Standard_Complex_Vectors.Vector(b'range);
      pdetU : natural32 := 0;
      mysols,mysols_ptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        mycell_ptr := Mixed_Cells_Queue.Next;
        exit when Is_Null(mycell_ptr);
        mic := Head_Of(mycell_ptr);
        if r = n then
          Select_Coefficients(c,e,mic.pts.all,s_c);
          Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
        else
          Select_Subsystem_to_Matrix_Format(c,e,mix,mic.pts.all,A,CC,b);
          Standard_Complex_Linear_Solvers.lufac(CC,n,piv,info);
          Standard_Complex_Linear_Solvers.lusolve(CC,n,piv,b);
        end if;
        U := A;
        Standard_Integer_Linear_Solvers.Upper_Triangulate(M,U);
        pdetU := Volume_of_Diagonal(U);
        Semaphore.Request(s_sols);    -- occupy solution space
        mysols := sols_ptr;
        for i in 1..pdetU loop
          sols_ptr := Tail_Of(sols_ptr);
        end loop;
       -- Semaphore.Release(s_sols);    -- end of second critical section
        brd := Standard_Radial_Solvers.Radii(b);
        bsc := Standard_Radial_Solvers.Scale(b,brd);
        Standard_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
        Semaphore.Release(s_sols);    -- end of second critical section
       -- make this end consistent with dd and qd ...
        logbrd := Standard_Radial_Solvers.Log10(brd);
        logx := Standard_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
        logx := Standard_Radial_Solvers.Multiply(M,logx);
        e10x := Standard_Radial_Solvers.Exp10(logx);
        mysols_ptr := mysols;
        for i in 1..pdetU loop
          ls := Head_Of(mysols_ptr);
          Standard_Binomial_Systems.Eval(M,ls.v,wrk);
          ls.v := wrk;
          Standard_Radial_Solvers.Multiply(ls.v,e10x);
          Track_Path(mix,lif,mic.nor,c,dpow(task_id),dctm(task_id),e,h,j,mf,ls);
          mysols_ptr := Tail_Of(mysols_ptr);
        end loop;
      end loop;
    end Next_Cell;
    procedure do_jobs is new Multitasking.Silent_Workers(Next_Cell);

    procedure Main is

      vol : Standard_Integer_Vectors.Vector(1..nt);
      mixvol : natural32;

    begin
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,cell_ptr,vol,mixvol);
      sols := Create(n,integer32(mixvol));
      sols_ptr := sols;
      Allocate_Workspace_for_Exponents(e,dpow);
      Allocate_Workspace_for_Coefficients(c,dctm);
      Mixed_Cells_Queue.Initialize(mcc);
      do_jobs(nt);
      for t in 1..nt loop
        Standard_Floating_VecVecs.Deep_Clear(dpow(t));
        Standard_Complex_VecVecs.Deep_Clear(dctm(t));
      end loop;
    end Main;

  begin
    Main;
  end Silent_Multitasking_Path_Tracker;

  procedure Reporting_Multitasking_Path_Tracker
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                h : in Standard_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                c : in Standard_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in Standard_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                mf : in Standard_Complex_Laur_JacoMats.Mult_Factors;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    cell_ptr : Mixed_Subdivision := mcc;
    s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : Standard_Complex_VecVecs.Array_of_VecVecs(1..nt);

    procedure Next_Cell ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   Task i processes the next cell. 

      myjob : natural32;
      mycell_ptr : Mixed_Subdivision;
      mic : Mixed_Cell;
      s_c : Standard_Complex_Vectors.Vector(1..2*n); -- fully mixed
      CC : Standard_Complex_Matrices.Matrix(1..n,1..n); -- for semi mixed
      piv : Standard_Integer_Vectors.Vector(1..n);
      info : integer32;
      A,M,U : Standard_Integer_Matrices.Matrix(q'range,q'range);
      b,wrk : Standard_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Standard_Floating_Vectors.Vector(b'range);
      bsc : Standard_Complex_Vectors.Vector(b'range);
      pdetU : natural32 := 0;
      mysols,mysols_ptr : Solution_List;
      ls : Link_to_Solution;

    --  procedure Check_Solution is
    --
    --    y : constant Standard_Complex_Vectors.Vector
    --      := Standard_Binomial_Systems.Eval(A,b,ls.v);
    --    r : constant double_float
    --      := Standard_Complex_Norms_Equals.Max_Norm(y);
    --    use Polyhedral_Coefficient_Homotopies;
    --    fp : Standard_Floating_VecVecs.VecVec(c'range)
    --       := Power_Transform(e,lif,mix,mic.nor.all);
    --    sq : Laur_Sys(q'range)
    --       := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
    --    yz : constant Standard_Complex_Vectors.Vector
    --       := Standard_Complex_Laur_SysFun.Eval(sq,ls.v);
    -- 
    --   begin
    --     put_line("The start solution evaluated : "); put_line(yz);
    --     put_line("The binomial start system :"); put(sq);
    --     put_line("The points in the cell :");
    --     Floating_Mixed_Subdivisions_io.put(mic.pts.all);
    --     put_line("The matrix A "); put(A);
    --     put_line("The coefficients :"); put_line(s_c);
    --   end Check_Solution;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      loop
        mycell_ptr := Mixed_Cells_Queue.Next;
        exit when Is_Null(mycell_ptr);
        myjob := natural32(Mixed_Cells_Queue.Next_Counter);
        put_line("task " & Multitasking.to_string(natural32(task_id))
                         & " computes cell "
                         & Multitasking.to_string(myjob));
        mic := Head_Of(mycell_ptr);
        if r = n then
          Select_Coefficients(c,e,mic.pts.all,s_c);
          Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
        else
          Select_Subsystem_to_Matrix_Format(c,e,mix,mic.pts.all,A,CC,b);
          Standard_Complex_Linear_Solvers.lufac(CC,n,piv,info);
          Standard_Complex_Linear_Solvers.lusolve(CC,n,piv,b);
        end if;
        U := A;
        Standard_Integer_Linear_Solvers.Upper_Triangulate(M,U);
        pdetU := Volume_of_Diagonal(U);
        Semaphore.Request(s_sols);    -- occupy solution space
        mysols := sols_ptr;
        for k in 1..pdetU loop
          sols_ptr := Tail_Of(sols_ptr);
        end loop;
       -- Semaphore.Release(s_sols);    -- end of second critical section
        brd := Standard_Radial_Solvers.Radii(b);
        bsc := Standard_Radial_Solvers.Scale(b,brd);
        Standard_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
        Semaphore.Release(s_sols);    -- end of second critical section
       -- no problems observed, but for dd/qd, the critical section
       -- had to be extended, so do so too here for consistency
        logbrd := Standard_Radial_Solvers.Log10(brd);
        logx := Standard_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
        logx := Standard_Radial_Solvers.Multiply(M,logx);
        e10x := Standard_Radial_Solvers.Exp10(logx);
        mysols_ptr := mysols;
        for k in 1..pdetU loop
          ls := Head_Of(mysols_ptr);
          Standard_Binomial_Systems.Eval(M,ls.v,wrk);
          ls.v := wrk;
          Standard_Radial_Solvers.Multiply(ls.v,e10x);
         -- Check_Solution;
          Track_Path(mix,lif,mic.nor,c,dpow(task_id),dctm(task_id),e,h,j,mf,ls);
          mysols_ptr := Tail_Of(mysols_ptr);
        end loop;
        put_line("task " & Multitasking.to_string(natural32(task_id))
                         & " computed "
                         & Multitasking.to_string(pdetU)
                         & " solutions");
      end loop;
    end Next_Cell;
    procedure do_jobs is new Multitasking.Silent_Workers(Next_Cell);

    procedure Main is

      vol : Standard_Integer_Vectors.Vector(1..nt);
      mixvol : natural32;

    begin
      new_line;
      put_line("computing volumes ...");
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,cell_ptr,vol,mixvol);
      new_line;
      put("volume distribution : "); put(vol); new_line;
      put("the mixed volume : "); put(mixvol,1); new_line;
      new_line;
      put_line("allocating solution lists and other data...");
      sols := Create(n,integer32(mixvol));
      sols_ptr := sols;
      Allocate_Workspace_for_Exponents(e,dpow);
      Allocate_Workspace_for_Coefficients(c,dctm);
      new_line;
      put_line("launching tasks ...");
      new_line;
      Mixed_Cells_Queue.Initialize(mcc);
      do_jobs(nt);
      for t in 1..nt loop
        Standard_Floating_VecVecs.Deep_Clear(dpow(t));
        Standard_Complex_VecVecs.Deep_Clear(dctm(t));
      end loop;
    end Main;

  begin
    Main;
  end Reporting_Multitasking_Path_Tracker;

  procedure Silent_Multitasking_Path_Tracker
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                h : in DoblDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                c : in DoblDobl_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in DoblDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                mf : in DoblDobl_Complex_Laur_JacoMats.Mult_Factors;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    cell_ptr : Mixed_Subdivision := mcc;
    s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : DoblDobl_Complex_VecVecs.Array_of_VecVecs(1..nt);

    procedure Next_Cell ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Task i processes the next cell.

      mycell_ptr : Mixed_Subdivision;
      mic : Mixed_Cell;
      s_c : DoblDobl_Complex_Vectors.Vector(1..2*n);    -- for fully mixed
      CC : DoblDobl_Complex_Matrices.Matrix(1..n,1..n); -- for semi mixed
      piv : Standard_Integer_Vectors.Vector(1..n);
      info : integer32;
      A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
      b,wrk : DoblDobl_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Double_Double_Vectors.Vector(b'range);
      bsc : DoblDobl_Complex_Vectors.Vector(b'range);
      pdetU : natural64 := 0;
      mysols,mysols_ptr : Solution_List;
      ls : Link_to_Solution;

     -- procedure Check_Solution is
    
     --   y : constant DoblDobl_Complex_Vectors.Vector
     --     := DoblDobl_Binomial_Systems.Eval(A,b,ls.v);
     --   r : constant double_double
     --     := DoblDobl_Complex_Vector_Norms.Max_Norm(y);
     --   residual : constant double_float := hi_part(r);
     --   use Polyhedral_Coefficient_Homotopies;
     --   fp : Standard_Floating_VecVecs.VecVec(c'range)
     --      := Power_Transform(e,lif,mix,mic.nor.all);
     --   sq : DoblDobl_Complex_Laur_Systems.Laur_Sys(q'range)
     --      := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
     --   yz : constant DoblDobl_Complex_Vectors.Vector
     --      := DoblDobl_Complex_Laur_SysFun.Eval(sq,ls.v);
     
     --  begin
     --    put_line("*** residual : " & convert(residual)
     --             & " by task " & Multitasking.to_string(task_id) & "***");
     --    if residual > 1.0E-8 then
     --      put_line("The binomial start system :"); put(sq);
     --      put_line("The start solution evaluated : "); put_line(yz);
     --      put_line("The points in the cell :");
     --      Floating_Mixed_Subdivisions_io.put(mic.pts.all);
     --      put_line("The matrix A "); put(A);
     --      put_line("The coefficients :"); put_line(s_c);
     --    end if;
     --  end Check_Solution;

    begin
      loop
        mycell_ptr := Mixed_Cells_Queue.Next;
        exit when Is_Null(mycell_ptr);
        mic := Head_Of(mycell_ptr);
        if r = n then
          Select_Coefficients(c,e,mic.pts.all,s_c);
          Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
        else
          Select_Subsystem_to_Matrix_Format(c,e,mix,mic.pts.all,A,CC,b);
          DoblDobl_Complex_Linear_Solvers.lufac(CC,n,piv,info);
          DoblDobl_Complex_Linear_Solvers.lusolve(CC,n,piv,b);
        end if;
        U := A;
        Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
        pdetU := Volume_of_Diagonal(U);
        Semaphore.Request(s_sols);    -- occupy solution space
        mysols := sols_ptr;
        for i in 1..pdetU loop
          sols_ptr := Tail_Of(sols_ptr);
        end loop;
       -- Semaphore.Release(s_sols);    -- end of second critical section
        brd := DoblDobl_Radial_Solvers.Radii(b);
        bsc := DoblDobl_Radial_Solvers.Scale(b,brd);
        DoblDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
        Semaphore.Release(s_sols); -- critical section must end here !?
        logbrd := DoblDobl_Radial_Solvers.Log10(brd);
        logx := DoblDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
        logx := DoblDobl_Radial_Solvers.Multiply(M,logx);
        e10x := DoblDobl_Radial_Solvers.Exp10(logx);
        mysols_ptr := mysols;
        for i in 1..pdetU loop
          ls := Head_Of(mysols_ptr);
          DoblDobl_Binomial_Systems.Eval(M,ls.v,wrk);
          ls.v := wrk;
          DoblDobl_Radial_Solvers.Multiply(ls.v,e10x);
         -- Check_Solution;
          Track_Path(mix,lif,mic.nor,c,dpow(task_id),dctm(task_id),e,h,j,mf,ls);
          mysols_ptr := Tail_Of(mysols_ptr);
        end loop;
      end loop;
    end Next_Cell;
    procedure do_jobs is new Multitasking.Silent_Workers(Next_Cell);

    procedure Main is

      vol : Standard_Integer_Vectors.Vector(1..nt);
      mixvol : natural32;

    begin
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,cell_ptr,vol,mixvol);
      sols := Create(n,integer32(mixvol));
      sols_ptr := sols;
      Allocate_Workspace_for_Exponents(e,dpow);
      Allocate_Workspace_for_Coefficients(c,dctm);
      Mixed_Cells_Queue.Initialize(mcc);
      do_jobs(nt);
      for t in 1..nt loop
        Standard_Floating_VecVecs.Deep_Clear(dpow(t));
        DoblDobl_Complex_VecVecs.Deep_Clear(dctm(t));
      end loop;
    end Main;

  begin
    Main;
  end Silent_Multitasking_Path_Tracker;

  procedure Reporting_Multitasking_Path_Tracker
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                h : in DoblDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                c : in DoblDobl_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in DoblDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                mf : in DoblDobl_Complex_Laur_JacoMats.Mult_Factors;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    cell_ptr : Mixed_Subdivision := mcc;
    s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : DoblDobl_Complex_VecVecs.Array_of_VecVecs(1..nt);

    procedure Next_Cell ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   Task i processes the next cell. 

      myjob : natural32;
      mycell_ptr : Mixed_Subdivision;
      mic : Mixed_Cell;
      s_c : DoblDobl_Complex_Vectors.Vector(1..2*n); -- fully mixed
      CC : DoblDobl_Complex_Matrices.Matrix(1..n,1..n); -- for semi mixed
      piv : Standard_Integer_Vectors.Vector(1..n);
      info : integer32;
      A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
      b,wrk : DoblDobl_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Double_Double_Vectors.Vector(b'range);
      bsc : DoblDobl_Complex_Vectors.Vector(b'range);
      pdetU : natural64 := 0;
      mysols,mysols_ptr : Solution_List;
      ls : Link_to_Solution;

     -- procedure Check_Solution is
    
     --   y : constant DoblDobl_Complex_Vectors.Vector
     --     := DoblDobl_Binomial_Systems.Eval(A,b,ls.v);
     --   r : constant double_double
     --     := DoblDobl_Complex_Vector_Norms.Max_Norm(y);
     --   use Polyhedral_Coefficient_Homotopies;
     --   fp : Standard_Floating_VecVecs.VecVec(c'range)
     --      := Power_Transform(e,lif,mix,mic.nor.all);
     --   sq : DoblDobl_Complex_Laur_Systems.Laur_Sys(q'range)
     --      := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
     --   yz : constant DoblDobl_Complex_Vectors.Vector
     --      := DoblDobl_Complex_Laur_SysFun.Eval(sq,ls.v);
     
     --  begin
     --    put_line("The start solution evaluated : "); put_line(yz);
     --    put_line("The binomial start system :"); put(sq);
     --    put_line("The points in the cell :");
     --    Floating_Mixed_Subdivisions_io.put(mic.pts.all);
     --    put_line("The matrix A "); put(A);
     --    put_line("The coefficients :"); put_line(s_c);
     --  end Check_Solution;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      loop
        mycell_ptr := Mixed_Cells_Queue.Next;
        exit when Is_Null(mycell_ptr);
        myjob := natural32(Mixed_Cells_Queue.Next_Counter);
        put_line("task " & Multitasking.to_string(natural32(task_id))
                         & " computes cell "
                         & Multitasking.to_string(myjob));
        mic := Head_Of(mycell_ptr);
        if r = n then
          Select_Coefficients(c,e,mic.pts.all,s_c);
          Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
        else
          Select_Subsystem_to_Matrix_Format(c,e,mix,mic.pts.all,A,CC,b);
          DoblDobl_Complex_Linear_Solvers.lufac(CC,n,piv,info);
          DoblDobl_Complex_Linear_Solvers.lusolve(CC,n,piv,b);
        end if;
        U := A;
        Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
        pdetU := Volume_of_Diagonal(U);
        Semaphore.Request(s_sols);    -- occupy solution space
        mysols := sols_ptr;
        for k in 1..pdetU loop
          sols_ptr := Tail_Of(sols_ptr);
        end loop;
       -- Semaphore.Release(s_sols);    -- end of second critical section
        brd := DoblDobl_Radial_Solvers.Radii(b);
        bsc := DoblDobl_Radial_Solvers.Scale(b,brd);
        DoblDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
        Semaphore.Release(s_sols); -- critical section must end here !?
        logbrd := DoblDobl_Radial_Solvers.Log10(brd);
        logx := DoblDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
        logx := DoblDobl_Radial_Solvers.Multiply(M,logx);
        e10x := DoblDobl_Radial_Solvers.Exp10(logx);
        mysols_ptr := mysols;
        for k in 1..pdetU loop
          ls := Head_Of(mysols_ptr);
          DoblDobl_Binomial_Systems.Eval(M,ls.v,wrk);
          ls.v := wrk;
          DoblDobl_Radial_Solvers.Multiply(ls.v,e10x);
         -- Check_Solution;
          Track_Path(mix,lif,mic.nor,c,dpow(task_id),dctm(task_id),e,h,j,mf,ls);
          mysols_ptr := Tail_Of(mysols_ptr);
        end loop;
        put_line("task " & Multitasking.to_string(natural32(task_id))
                         & " computed "
                         & Multitasking.to_string(natural32(pdetU))
                         & " solutions");
      end loop;
    end Next_Cell;
    procedure do_jobs is new Multitasking.Silent_Workers(Next_Cell);

    procedure Main is

      vol : Standard_Integer_Vectors.Vector(1..nt);
      mixvol : natural32;

    begin
      new_line;
      put_line("computing volumes ...");
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,cell_ptr,vol,mixvol);
      new_line;
      put("volume distribution : "); put(vol); new_line;
      put("the mixed volume : "); put(mixvol,1); new_line;
      new_line;
      put_line("allocating solution lists and other data...");
      sols := Create(n,integer32(mixvol));
      sols_ptr := sols;
      Allocate_Workspace_for_Exponents(e,dpow);
      Allocate_Workspace_for_Coefficients(c,dctm);
      new_line;
      put_line("launching tasks ...");
      new_line;
      Mixed_Cells_Queue.Initialize(mcc);
      do_jobs(nt);
      for t in 1..nt loop
        Standard_Floating_VecVecs.Deep_Clear(dpow(t));
        DoblDobl_Complex_VecVecs.Deep_Clear(dctm(t));
      end loop;
    end Main;

  begin
    Main;
  end Reporting_Multitasking_Path_Tracker;

  procedure Silent_Multitasking_Path_Tracker
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                h : in QuadDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                c : in QuadDobl_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in QuadDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                mf : in QuadDobl_Complex_Laur_JacoMats.Mult_Factors;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    cell_ptr : Mixed_Subdivision := mcc;
    s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : QuadDobl_Complex_VecVecs.Array_of_VecVecs(1..nt);

    procedure Next_Cell ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Task i processes the next cell.

      mycell_ptr : Mixed_Subdivision;
      mic : Mixed_Cell;
      s_c : QuadDobl_Complex_Vectors.Vector(1..2*n);    -- for fully mixed
      CC : QuadDobl_Complex_Matrices.Matrix(1..n,1..n); -- for semi mixed
      piv : Standard_Integer_Vectors.Vector(1..n);
      info : integer32;
      A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
      b,wrk : QuadDobl_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Quad_Double_Vectors.Vector(b'range);
      bsc : QuadDobl_Complex_Vectors.Vector(b'range);
      pdetU : natural64 := 0;
      mysols,mysols_ptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        mycell_ptr := Mixed_Cells_Queue.Next;
        exit when Is_Null(mycell_ptr);
        mic := Head_Of(mycell_ptr);
        if r = n then
          Select_Coefficients(c,e,mic.pts.all,s_c);
          Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
        else
          Select_Subsystem_to_Matrix_Format(c,e,mix,mic.pts.all,A,CC,b);
          QuadDobl_Complex_Linear_Solvers.lufac(CC,n,piv,info);
          QuadDobl_Complex_Linear_Solvers.lusolve(CC,n,piv,b);
        end if;
        U := A;
        Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
        pdetU := Volume_of_Diagonal(U);
        Semaphore.Request(s_sols);    -- occupy solution space
        mysols := sols_ptr;
        for i in 1..pdetU loop
          sols_ptr := Tail_Of(sols_ptr);
        end loop;
       -- Semaphore.Release(s_sols);    -- end of second critical section
        brd := QuadDobl_Radial_Solvers.Radii(b);
        bsc := QuadDobl_Radial_Solvers.Scale(b,brd);
        QuadDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
        Semaphore.Release(s_sols); -- critical section must end here !?
        logbrd := QuadDobl_Radial_Solvers.Log10(brd);
        logx := QuadDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
        logx := QuadDobl_Radial_Solvers.Multiply(M,logx);
        e10x := QuadDobl_Radial_Solvers.Exp10(logx);
        mysols_ptr := mysols;
        for i in 1..pdetU loop
          ls := Head_Of(mysols_ptr);
          QuadDobl_Binomial_Systems.Eval(M,ls.v,wrk);
          ls.v := wrk;
          QuadDobl_Radial_Solvers.Multiply(ls.v,e10x);
          Track_Path(mix,lif,mic.nor,c,dpow(task_id),dctm(task_id),e,h,j,mf,ls);
          mysols_ptr := Tail_Of(mysols_ptr);
        end loop;
      end loop;
    end Next_Cell;
    procedure do_jobs is new Multitasking.Silent_Workers(Next_Cell);

    procedure Main is

      vol : Standard_Integer_Vectors.Vector(1..nt);
      mixvol : natural32;

    begin
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,cell_ptr,vol,mixvol);
      sols := Create(n,integer32(mixvol));
      sols_ptr := sols;
      Allocate_Workspace_for_Exponents(e,dpow);
      Allocate_Workspace_for_Coefficients(c,dctm);
      Mixed_Cells_Queue.Initialize(mcc);
      do_jobs(nt);
      for t in 1..nt loop
        Standard_Floating_VecVecs.Deep_Clear(dpow(t));
        QuadDobl_Complex_VecVecs.Deep_Clear(dctm(t));
      end loop;
    end Main;

  begin
    Main;
  end Silent_Multitasking_Path_Tracker;

  procedure Reporting_Multitasking_Path_Tracker
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                h : in QuadDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                c : in QuadDobl_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in QuadDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                mf : in QuadDobl_Complex_Laur_JacoMats.Mult_Factors;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    cell_ptr : Mixed_Subdivision := mcc;
    s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : QuadDobl_Complex_VecVecs.Array_of_VecVecs(1..nt);

    procedure Next_Cell ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   Task i processes the next cell. 

      myjob : natural32;
      mycell_ptr : Mixed_Subdivision;
      mic : Mixed_Cell;
      s_c : QuadDobl_Complex_Vectors.Vector(1..2*n); -- fully mixed
      CC : QuadDobl_Complex_Matrices.Matrix(1..n,1..n); -- for semi mixed
      piv : Standard_Integer_Vectors.Vector(1..n);
      info : integer32;
      A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
      b,wrk : QuadDobl_Complex_Vectors.Vector(q'range);
      brd,logbrd,logx,e10x : Quad_Double_Vectors.Vector(b'range);
      bsc : QuadDobl_Complex_Vectors.Vector(b'range);
      pdetU : natural64 := 0;
      mysols,mysols_ptr : Solution_List;
      ls : Link_to_Solution;

    --  procedure Check_Solution is
    --
    --    y : constant Standard_Complex_Vectors.Vector
    --      := Standard_Binomial_Systems.Eval(A,b,ls.v);
    --    r : constant double_float
    --      := Standard_Complex_Norms_Equals.Max_Norm(y);
    --    use Polyhedral_Coefficient_Homotopies;
    --    fp : Standard_Floating_VecVecs.VecVec(c'range)
    --       := Power_Transform(e,lif,mix,mic.nor.all);
    --    sq : Laur_Sys(q'range)
    --       := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
    --    yz : constant Standard_Complex_Vectors.Vector
    --       := Standard_Complex_Laur_SysFun.Eval(sq,ls.v);
    -- 
    --   begin
    --     put_line("The start solution evaluated : "); put_line(yz);
    --     put_line("The binomial start system :"); put(sq);
    --     put_line("The points in the cell :");
    --     Floating_Mixed_Subdivisions_io.put(mic.pts.all);
    --     put_line("The matrix A "); put(A);
    --     put_line("The coefficients :"); put_line(s_c);
    --   end Check_Solution;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      loop
        mycell_ptr := Mixed_Cells_Queue.Next;
        exit when Is_Null(mycell_ptr);
        myjob := natural32(Mixed_Cells_Queue.Next_Counter);
        put_line("task " & Multitasking.to_string(natural32(task_id))
                         & " computes cell "
                         & Multitasking.to_string(myjob));
        mic := Head_Of(mycell_ptr);
        if r = n then
          Select_Coefficients(c,e,mic.pts.all,s_c);
          Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
        else
          Select_Subsystem_to_Matrix_Format(c,e,mix,mic.pts.all,A,CC,b);
          QuadDobl_Complex_Linear_Solvers.lufac(CC,n,piv,info);
          QuadDobl_Complex_Linear_Solvers.lusolve(CC,n,piv,b);
        end if;
        U := A;
        Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
        pdetU := Volume_of_Diagonal(U);
        Semaphore.Request(s_sols);    -- occupy solution space
        mysols := sols_ptr;
        for k in 1..pdetU loop
          sols_ptr := Tail_Of(sols_ptr);
        end loop;
       -- Semaphore.Release(s_sols);    -- end of second critical section
        brd := QuadDobl_Radial_Solvers.Radii(b);
        bsc := QuadDobl_Radial_Solvers.Scale(b,brd);
        QuadDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,mysols);
        Semaphore.Release(s_sols); -- critical section must end here !?
        logbrd := QuadDobl_Radial_Solvers.Log10(brd);
        logx := QuadDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
        logx := QuadDobl_Radial_Solvers.Multiply(M,logx);
        e10x := QuadDobl_Radial_Solvers.Exp10(logx);
        mysols_ptr := mysols;
        for k in 1..pdetU loop
          ls := Head_Of(mysols_ptr);
          QuadDobl_Binomial_Systems.Eval(M,ls.v,wrk);
          ls.v := wrk;
          QuadDobl_Radial_Solvers.Multiply(ls.v,e10x);
         -- Check_Solution;
          Track_Path(mix,lif,mic.nor,c,dpow(task_id),dctm(task_id),e,h,j,mf,ls);
          mysols_ptr := Tail_Of(mysols_ptr);
        end loop;
        put_line("task " & Multitasking.to_string(natural32(task_id))
                         & " computed "
                         & Multitasking.to_string(natural32(pdetU))
                         & " solutions");
      end loop;
    end Next_Cell;
    procedure do_jobs is new Multitasking.Silent_Workers(Next_Cell);

    procedure Main is

      vol : Standard_Integer_Vectors.Vector(1..nt);
      mixvol : natural32;

    begin
      new_line;
      put_line("computing volumes ...");
      Silent_Static_Multithreaded_Mixed_Volume(nt,n,cell_ptr,vol,mixvol);
      new_line;
      put("volume distribution : "); put(vol); new_line;
      put("the mixed volume : "); put(mixvol,1); new_line;
      new_line;
      put_line("allocating solution lists and other data...");
      sols := Create(n,integer32(mixvol));
      sols_ptr := sols;
      Allocate_Workspace_for_Exponents(e,dpow);
      Allocate_Workspace_for_Coefficients(c,dctm);
      new_line;
      put_line("launching tasks ...");
      new_line;
      Mixed_Cells_Queue.Initialize(mcc);
      do_jobs(nt);
      for t in 1..nt loop
        Standard_Floating_VecVecs.Deep_Clear(dpow(t));
        QuadDobl_Complex_VecVecs.Deep_Clear(dctm(t));
      end loop;
    end Main;

  begin
    Main;
  end Reporting_Multitasking_Path_Tracker;

  procedure Silent_Multitasking_Path_Tracker
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Laur_SysFun,Standard_Complex_Laur_JacoMats;

    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(q'range);
    pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists(q'range);
    lif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    h : constant Eval_Coeff_Laur_Sys(q'range) := Create(q);
    c : Standard_Complex_VecVecs.VecVec(q'range) := Coeff(q);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    e := Exponent_Vectors.Create(q);
    Create(q,j,mf);
    Silent_Multitasking_Path_Tracker(q,nt,n,m,mix,lif,mcc,h,c,e,j,mf,sols);
  end Silent_Multitasking_Path_Tracker;

  procedure Reporting_Multitasking_Path_Tracker
              ( file : in file_type;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Laur_SysFun,Standard_Complex_Laur_JacoMats;

    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(q'range);
    pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists(q'range);
    lif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    h : constant Eval_Coeff_Laur_Sys(q'range) := Create(q);
    c : Standard_Complex_VecVecs.VecVec(q'range) := Coeff(q);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    e := Exponent_Vectors.Create(q);
    Create(q,j,mf);
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    Floating_Mixed_Subdivisions_io.put(file,lif);
    Reporting_Multitasking_Path_Tracker(q,nt,n,m,mix,lif,mcc,h,c,e,j,mf,sols);
  end Reporting_Multitasking_Path_Tracker;

  procedure Silent_Multitasking_Path_Tracker
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Laur_SysFun,DoblDobl_Complex_Laur_JacoMats;

    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(q'range);
    pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists(q'range);
    lif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    h : constant Eval_Coeff_Laur_Sys(q'range) := Create(q);
    c : DoblDobl_Complex_VecVecs.VecVec(q'range) := Coeff(q);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    e := Exponent_Vectors.Create(q);
    Create(q,j,mf);
    Silent_Multitasking_Path_Tracker(q,nt,n,m,mix,lif,mcc,h,c,e,j,mf,sols);
  end Silent_Multitasking_Path_Tracker;

  procedure Reporting_Multitasking_Path_Tracker
              ( file : in file_type;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Laur_SysFun,DoblDobl_Complex_Laur_JacoMats;

    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(q'range);
    pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists(q'range);
    lif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    h : constant Eval_Coeff_Laur_Sys(q'range) := Create(q);
    c : DoblDobl_Complex_VecVecs.VecVec(q'range) := Coeff(q);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    e := Exponent_Vectors.Create(q);
    Create(q,j,mf);
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    Floating_Mixed_Subdivisions_io.put(file,lif);
    Reporting_Multitasking_Path_Tracker(q,nt,n,m,mix,lif,mcc,h,c,e,j,mf,sols);
  end Reporting_Multitasking_Path_Tracker;

  procedure Silent_Multitasking_Path_Tracker
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Laur_SysFun,QuadDobl_Complex_Laur_JacoMats;

    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(q'range);
    pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists(q'range);
    lif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    h : constant Eval_Coeff_Laur_Sys(q'range) := Create(q);
    c : QuadDobl_Complex_VecVecs.VecVec(q'range) := Coeff(q);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    e := Exponent_Vectors.Create(q);
    Create(q,j,mf);
    Silent_Multitasking_Path_Tracker(q,nt,n,m,mix,lif,mcc,h,c,e,j,mf,sols);
  end Silent_Multitasking_Path_Tracker;

  procedure Reporting_Multitasking_Path_Tracker
              ( file : in file_type;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Laur_SysFun,QuadDobl_Complex_Laur_JacoMats;

    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(q'range);
    pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists(q'range);
    lif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    h : constant Eval_Coeff_Laur_Sys(q'range) := Create(q);
    c : QuadDobl_Complex_VecVecs.VecVec(q'range) := Coeff(q);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    e := Exponent_Vectors.Create(q);
    Create(q,j,mf);
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    Floating_Mixed_Subdivisions_io.put(file,lif);
    Reporting_Multitasking_Path_Tracker(q,nt,n,m,mix,lif,mcc,h,c,e,j,mf,sols);
  end Reporting_Multitasking_Path_Tracker;

end Multitasking_Polyhedral_Trackers;
