with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Norms_Equals;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
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
with Standard_Floating_VecVecs;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_Functions;
with Continuation_Parameters;
with Standard_Continuation_Data;
with Standard_Path_Trackers;
with DoblDobl_Continuation_Data;
with DoblDobl_Path_Trackers;
with QuadDobl_Continuation_Data;
with QuadDobl_Path_Trackers;
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
with Polyhedral_Coefficient_Homotopies;
with Multitasking,Semaphore;
with Multitasking_Volume_Computation;    use Multitasking_Volume_Computation;
with Polyhedral_Start_Systems;           use Polyhedral_Start_Systems;
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

    use Standard_Complex_Numbers,Standard_Complex_Solutions;
    use Standard_Complex_Laur_SysFun,Standard_Complex_Laur_JacoMats;

    cell_ptr : Mixed_Subdivision := mcc;
    first : boolean := true;
    s_cell,s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : Standard_Complex_VecVecs.Array_of_VecVecs(1..nt);
    t1 : constant Complex_Number := Create(1.0);
    tol : constant double_float := 1.0E-12;
   -- tol_zero : constant double_float := 1.0E-8;
    pp1 : Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_for_Path;
    pp2 : constant Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_End_Game;
    cp1 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_for_Path;
    cp2 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_End_Game;

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

      procedure Call_Path_Tracker is

      -- DESCRIPTION :
      --   Calls the path tracker on the solution ls points to.

        use Standard_Complex_Norms_Equals;
        use Standard_Continuation_Data;
        use Standard_Path_Trackers;
        use Polyhedral_Coefficient_Homotopies;

        function Eval ( x : Standard_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return Standard_Complex_Vectors.Vector is
        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          return Eval(h,dctm(task_id).all,x);
        end Eval;

        function dHt ( x : Standard_Complex_Vectors.Vector;
                       t : Standard_Complex_Numbers.Complex_Number )
                     return Standard_Complex_Vectors.Vector is

          res : Standard_Complex_Vectors.Vector(h'range);
          xtl : constant integer32 := x'last+1;

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in wrk'range loop
            res(k) := Eval(j(k,xtl),mf(k,xtl).all,dctm(task_id)(k).all,x);
          end loop;
          return res;
        end dHt;

        function dHx ( x : Standard_Complex_Vectors.Vector;
                       t : Complex_Number )
                     return Standard_Complex_Matrices.Matrix is

          mt : Standard_Complex_Matrices.Matrix(x'range,x'range);

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in mt'range(1) loop
            for L in mt'range(2) loop
              mt(k,L) := Eval(j(k,L),mf(k,L).all,dctm(task_id)(k).all,x);
            end loop;
          end loop;
          return mt;
        end dHx;

        procedure Track_Path_along_Path is
          new Linear_Single_Normal_Silent_Continue(Max_Norm,Eval,dHt,dHx);

        procedure Track_Path_at_End is
           new Linear_Single_Conditioned_Silent_Continue
                 (Max_Norm,Eval,dHt,dHx);

      begin
        Power_Transform(e,lif,mix,mic.nor.all,dpow(task_id).all);
        Polyhedral_Coefficient_Homotopies.Scale(dpow(task_id).all);
        ls.t := Create(0.0);
        declare
          s : Solu_Info := Shallow_Create(ls);
          v : Standard_Floating_Vectors.Link_to_Vector;
          e : double_float := 0.0;
        begin
          pp1.dist_target := 0.0;
          Track_Path_along_Path(s,t1,tol,false,pp1,cp1);
          Track_Path_at_End(s,t1,tol,false,0,v,e,pp2,cp2);
          ls.err := s.cora; ls.rco := s.rcond; ls.res := s.resa;
        end;
      end Call_Path_Tracker;

    begin
      loop
        Semaphore.Request(s_cell);  -- request new cell
        if first then
          first := false;
        else
          if not Is_Null(cell_ptr)
           then cell_ptr := Tail_Of(cell_ptr);
          end if;
        end if;
        mycell_ptr := cell_ptr;
        Semaphore.Release(s_cell);  -- end of first critical section
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
          Call_Path_Tracker;
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
      for t in 1..nt loop
        dpow(t) := new Standard_Floating_VecVecs.VecVec(e'range);
        dctm(t) := new Standard_Complex_VecVecs.VecVec(c'range);
        for k in dpow(t)'range loop
          dpow(t)(k) := new Standard_Floating_Vectors.Vector(e(k)'range);
          dctm(t)(k) := new Standard_Complex_Vectors.Vector(c(k)'range);
        end loop;
      end loop;
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

    use Standard_Complex_Numbers,Standard_Complex_Solutions;
    use Standard_Complex_Laur_SysFun,Standard_Complex_Laur_JacoMats;

    cell_ptr : Mixed_Subdivision := mcc;
    cell_ind : integer32 := 0;
    s_cell,s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : Standard_Complex_VecVecs.Array_of_VecVecs(1..nt);
    t1 : constant Complex_Number := Create(1.0);
    tol : constant double_float := 1.0E-12;
   -- tol_zero : constant double_float := 1.0E-8;
    pp1 : Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_for_Path;
    pp2 : constant Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_End_Game;
    cp1 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_for_Path;
    cp2 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_End_Game;

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

      procedure Call_Path_Tracker is

      -- DESCRIPTION :
      --   Calls the path tracker on the solution ls points to.

        use Standard_Complex_Norms_Equals;
        use Standard_Continuation_Data;
        use Standard_Path_Trackers;
        use Polyhedral_Coefficient_Homotopies;

        function Eval ( x : Standard_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return Standard_Complex_Vectors.Vector is
        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          return Eval(h,dctm(task_id).all,x);
        end Eval;

        function dHt ( x : Standard_Complex_Vectors.Vector;
                       t : Complex_Number )
                     return Standard_Complex_Vectors.Vector is

          res : Standard_Complex_Vectors.Vector(h'range);
          xtl : constant integer32 := x'last+1;

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in wrk'range loop
            res(k) := Eval(j(k,xtl),mf(k,xtl).all,dctm(task_id)(k).all,x);
          end loop;
          return res;
        end dHt;

        function dHx ( x : Standard_Complex_Vectors.Vector;
                       t : Complex_Number )
                     return Standard_Complex_Matrices.Matrix is

          mt : Standard_Complex_Matrices.Matrix(x'range,x'range);

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in mt'range(1) loop
            for L in mt'range(2) loop
              mt(k,L) := Eval(j(k,L),mf(k,L).all,dctm(task_id)(k).all,x);
            end loop;
          end loop;
          return mt;
        end dHx;

        procedure Track_Path_along_Path is
          new Linear_Single_Normal_Silent_Continue(Max_Norm,Eval,dHt,dHx);

        procedure Track_Path_at_End is
           new Linear_Single_Conditioned_Silent_Continue
                 (Max_Norm,Eval,dHt,dHx);

      begin
        Power_Transform(e,lif,mix,mic.nor.all,dpow(task_id).all);
        Polyhedral_Coefficient_Homotopies.Scale(dpow(task_id).all);
        ls.t := Create(0.0);
        declare
          s : Solu_Info := Shallow_Create(ls);
          v : Standard_Floating_Vectors.Link_to_Vector;
          e : double_float := 0.0;
        begin
          pp1.dist_target := 0.0;
          Track_Path_along_Path(s,t1,tol,false,pp1,cp1);
          Track_Path_at_End(s,t1,tol,false,0,v,e,pp2,cp2);
          ls.err := s.cora; ls.rco := s.rcond; ls.res := s.resa;
        end;
      end Call_Path_Tracker;

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
        Semaphore.Request(s_cell);   -- take next cell
        if cell_ind = 0 then
          cell_ind := 1;
        else
          cell_ind := cell_ind + 1;
          if not Is_Null(cell_ptr)
           then cell_ptr := Tail_Of(cell_ptr);
          end if;
        end if;
        mycell_ptr := cell_ptr;
        myjob := natural32(cell_ind);
        Semaphore.Release(s_cell);   -- end of first critical section
        exit when Is_Null(mycell_ptr);
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
          Call_Path_Tracker;
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
      for t in 1..nt loop
        dpow(t) := new Standard_Floating_VecVecs.VecVec(e'range);
        dctm(t) := new Standard_Complex_VecVecs.VecVec(c'range);
        for k in dpow(t)'range loop
          dpow(t)(k) := new Standard_Floating_Vectors.Vector(e(k)'range);
          dctm(t)(k) := new Standard_Complex_Vectors.Vector(c(k)'range);
        end loop;
      end loop;
      sols_ptr := sols;
      new_line;
      put_line("launching tasks ...");
      new_line;
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

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Laur_SysFun,DoblDobl_Complex_Laur_JacoMats;

    cell_ptr : Mixed_Subdivision := mcc;
    first : boolean := true;
    s_cell,s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : DoblDobl_Complex_VecVecs.Array_of_VecVecs(1..nt);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    tol : constant double_float := 1.0E-12;
   -- tol_zero : constant double_float := 1.0E-8;
    pp1 : Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_for_Path;
    pp2 : constant Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_End_Game;
    cp1 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_for_Path;
    cp2 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_End_Game;

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

      procedure Call_Path_Tracker is

      -- DESCRIPTION :
      --   Calls the path tracker on the solution ls points to.

        use DoblDobl_Complex_Vector_Norms;
        use DoblDobl_Continuation_Data;
        use DoblDobl_Path_Trackers;
        use Polyhedral_Coefficient_Homotopies;

        function Eval ( x : DoblDobl_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return DoblDobl_Complex_Vectors.Vector is
        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          return Eval(h,dctm(task_id).all,x);
        end Eval;

        function dHt ( x : DoblDobl_Complex_Vectors.Vector;
                       t : DoblDobl_Complex_Numbers.Complex_Number )
                     return DoblDobl_Complex_Vectors.Vector is

          res : DoblDobl_Complex_Vectors.Vector(h'range);
          xtl : constant integer32 := x'last+1;

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in wrk'range loop
            res(k) := Eval(j(k,xtl),mf(k,xtl).all,dctm(task_id)(k).all,x);
          end loop;
          return res;
        end dHt;

        function dHx ( x : DoblDobl_Complex_Vectors.Vector;
                       t : Complex_Number )
                     return DoblDobl_Complex_Matrices.Matrix is

          mt : DoblDobl_Complex_Matrices.Matrix(x'range,x'range);

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in mt'range(1) loop
            for L in mt'range(2) loop
              mt(k,L) := Eval(j(k,L),mf(k,L).all,dctm(task_id)(k).all,x);
            end loop;
          end loop;
          return mt;
        end dHx;

        procedure Track_Path_along_Path is
          new Linear_Single_Normal_Silent_Continue(Max_Norm,Eval,dHt,dHx);

        procedure Track_Path_at_End is
           new Linear_Single_Conditioned_Silent_Continue
                 (Max_Norm,Eval,dHt,dHx);

      begin
        Power_Transform(e,lif,mix,mic.nor.all,dpow(task_id).all);
        Polyhedral_Coefficient_Homotopies.Scale(dpow(task_id).all);
        ls.t := Create(zero);
        declare
          s : Solu_Info := Shallow_Create(ls);
          v : Double_Double_Vectors.Link_to_Vector;
          e : double_double := zero;
        begin
          pp1.dist_target := 0.0;
          Track_Path_along_Path(s,t1,tol,false,pp1,cp1);
          Track_Path_at_End(s,t1,tol,false,0,v,e,pp2,cp2);
          ls.err := create(s.cora);
          ls.rco := create(s.rcond); 
          ls.res := create(s.resa);
        end;
      end Call_Path_Tracker;

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
        Semaphore.Request(s_cell);  -- request new cell
        if first then
          first := false;
        else
          if not Is_Null(cell_ptr)
           then cell_ptr := Tail_Of(cell_ptr);
          end if;
        end if;
        mycell_ptr := cell_ptr;
        Semaphore.Release(s_cell);  -- end of first critical section
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
          Call_Path_Tracker;
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
      for t in 1..nt loop
        dpow(t) := new Standard_Floating_VecVecs.VecVec(e'range);
        dctm(t) := new DoblDobl_Complex_VecVecs.VecVec(c'range);
        for k in dpow(t)'range loop
          dpow(t)(k) := new Standard_Floating_Vectors.Vector(e(k)'range);
          dctm(t)(k) := new DoblDobl_Complex_Vectors.Vector(c(k)'range);
        end loop;
      end loop;
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

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Laur_SysFun,DoblDobl_Complex_Laur_JacoMats;

    cell_ptr : Mixed_Subdivision := mcc;
    cell_ind : integer32 := 0;
    s_cell,s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : DoblDobl_Complex_VecVecs.Array_of_VecVecs(1..nt);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    tol : constant double_float := 1.0E-12;
   -- tol_zero : constant double_float := 1.0E-8;
    pp1 : Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_for_Path;
    pp2 : constant Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_End_Game;
    cp1 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_for_Path;
    cp2 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_End_Game;

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

      procedure Call_Path_Tracker is

      -- DESCRIPTION :
      --   Calls the path tracker on the solution ls points to.

        use DoblDobl_Complex_Vector_Norms;
        use DoblDobl_Continuation_Data;
        use DoblDobl_Path_Trackers;
        use Polyhedral_Coefficient_Homotopies;

        function Eval ( x : DoblDobl_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return DoblDobl_Complex_Vectors.Vector is
        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          return Eval(h,dctm(task_id).all,x);
        end Eval;

        function dHt ( x : DoblDobl_Complex_Vectors.Vector;
                       t : Complex_Number )
                     return DoblDobl_Complex_Vectors.Vector is

          res : DoblDobl_Complex_Vectors.Vector(h'range);
          xtl : constant integer32 := x'last+1;

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in wrk'range loop
            res(k) := Eval(j(k,xtl),mf(k,xtl).all,dctm(task_id)(k).all,x);
          end loop;
          return res;
        end dHt;

        function dHx ( x : DoblDobl_Complex_Vectors.Vector;
                       t : Complex_Number )
                     return DoblDobl_Complex_Matrices.Matrix is

          mt : DoblDobl_Complex_Matrices.Matrix(x'range,x'range);

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in mt'range(1) loop
            for L in mt'range(2) loop
              mt(k,L) := Eval(j(k,L),mf(k,L).all,dctm(task_id)(k).all,x);
            end loop;
          end loop;
          return mt;
        end dHx;

        procedure Track_Path_along_Path is
          new Linear_Single_Normal_Silent_Continue(Max_Norm,Eval,dHt,dHx);

        procedure Track_Path_at_End is
           new Linear_Single_Conditioned_Silent_Continue
                 (Max_Norm,Eval,dHt,dHx);

      begin
        Power_Transform(e,lif,mix,mic.nor.all,dpow(task_id).all);
        Polyhedral_Coefficient_Homotopies.Scale(dpow(task_id).all);
        ls.t := Create(zero);
        declare
          s : Solu_Info := Shallow_Create(ls);
          v : Double_Double_Vectors.Link_to_Vector;
          e : double_double := zero;
        begin
          pp1.dist_target := 0.0;
          Track_Path_along_Path(s,t1,tol,false,pp1,cp1);
          Track_Path_at_End(s,t1,tol,false,0,v,e,pp2,cp2);
          ls.err := create(s.cora);
          ls.rco := create(s.rcond);
          ls.res := create(s.resa);
        end;
      end Call_Path_Tracker;

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
        Semaphore.Request(s_cell);   -- take next cell
        if cell_ind = 0 then
          cell_ind := 1;
        else
          cell_ind := cell_ind + 1;
          if not Is_Null(cell_ptr)
           then cell_ptr := Tail_Of(cell_ptr);
          end if;
        end if;
        mycell_ptr := cell_ptr;
        myjob := natural32(cell_ind);
        Semaphore.Release(s_cell);   -- end of first critical section
        exit when Is_Null(mycell_ptr);
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
          Call_Path_Tracker;
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
      for t in 1..nt loop
        dpow(t) := new Standard_Floating_VecVecs.VecVec(e'range);
        dctm(t) := new DoblDobl_Complex_VecVecs.VecVec(c'range);
        for k in dpow(t)'range loop
          dpow(t)(k) := new Standard_Floating_Vectors.Vector(e(k)'range);
          dctm(t)(k) := new DoblDobl_Complex_Vectors.Vector(c(k)'range);
        end loop;
      end loop;
      sols_ptr := sols;
      new_line;
      put_line("launching tasks ...");
      new_line;
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

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Laur_SysFun,QuadDobl_Complex_Laur_JacoMats;

    cell_ptr : Mixed_Subdivision := mcc;
    first : boolean := true;
    s_cell,s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : QuadDobl_Complex_VecVecs.Array_of_VecVecs(1..nt);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    tol : constant double_float := 1.0E-12;
   -- tol_zero : constant double_float := 1.0E-8;
    pp1 : Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_for_Path;
    pp2 : constant Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_End_Game;
    cp1 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_for_Path;
    cp2 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_End_Game;

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

      procedure Call_Path_Tracker is

      -- DESCRIPTION :
      --   Calls the path tracker on the solution ls points to.

        use QuadDobl_Complex_Vector_Norms;
        use QuadDobl_Continuation_Data;
        use QuadDobl_Path_Trackers;
        use Polyhedral_Coefficient_Homotopies;

        function Eval ( x : QuadDobl_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return QuadDobl_Complex_Vectors.Vector is
        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          return Eval(h,dctm(task_id).all,x);
        end Eval;

        function dHt ( x : QuadDobl_Complex_Vectors.Vector;
                       t : QuadDobl_Complex_Numbers.Complex_Number )
                     return QuadDobl_Complex_Vectors.Vector is

          res : QuadDobl_Complex_Vectors.Vector(h'range);
          xtl : constant integer32 := x'last+1;

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in wrk'range loop
            res(k) := Eval(j(k,xtl),mf(k,xtl).all,dctm(task_id)(k).all,x);
          end loop;
          return res;
        end dHt;

        function dHx ( x : QuadDobl_Complex_Vectors.Vector;
                       t : Complex_Number )
                     return QuadDobl_Complex_Matrices.Matrix is

          mt : QuadDobl_Complex_Matrices.Matrix(x'range,x'range);

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in mt'range(1) loop
            for L in mt'range(2) loop
              mt(k,L) := Eval(j(k,L),mf(k,L).all,dctm(task_id)(k).all,x);
            end loop;
          end loop;
          return mt;
        end dHx;

        procedure Track_Path_along_Path is
          new Linear_Single_Normal_Silent_Continue(Max_Norm,Eval,dHt,dHx);

        procedure Track_Path_at_End is
           new Linear_Single_Conditioned_Silent_Continue
                 (Max_Norm,Eval,dHt,dHx);

      begin
        Power_Transform(e,lif,mix,mic.nor.all,dpow(task_id).all);
        Polyhedral_Coefficient_Homotopies.Scale(dpow(task_id).all);
        ls.t := Create(zero);
        declare
          s : Solu_Info := Shallow_Create(ls);
          v : Quad_Double_Vectors.Link_to_Vector;
          e : quad_double := zero;
        begin
          pp1.dist_target := 0.0;
          Track_Path_along_Path(s,t1,tol,false,pp1,cp1);
          Track_Path_at_End(s,t1,tol,false,0,v,e,pp2,cp2);
          ls.err := create(s.cora);
          ls.rco := create(s.rcond); 
          ls.res := create(s.resa);
        end;
      end Call_Path_Tracker;

    begin
      loop
        Semaphore.Request(s_cell);  -- request new cell
        if first then
          first := false;
        else
          if not Is_Null(cell_ptr)
           then cell_ptr := Tail_Of(cell_ptr);
          end if;
        end if;
        mycell_ptr := cell_ptr;
        Semaphore.Release(s_cell);  -- end of first critical section
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
          Call_Path_Tracker;
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
      for t in 1..nt loop
        dpow(t) := new Standard_Floating_VecVecs.VecVec(e'range);
        dctm(t) := new QuadDobl_Complex_VecVecs.VecVec(c'range);
        for k in dpow(t)'range loop
          dpow(t)(k) := new Standard_Floating_Vectors.Vector(e(k)'range);
          dctm(t)(k) := new QuadDobl_Complex_Vectors.Vector(c(k)'range);
        end loop;
      end loop;
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

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Laur_SysFun,QuadDobl_Complex_Laur_JacoMats;

    cell_ptr : Mixed_Subdivision := mcc;
    cell_ind : integer32 := 0;
    s_cell,s_sols : Semaphore.Lock;
    sols_ptr : Solution_List;
    dpow : Standard_Floating_VecVecs.Array_of_VecVecs(1..nt);
    dctm : QuadDobl_Complex_VecVecs.Array_of_VecVecs(1..nt);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    tol : constant double_float := 1.0E-12;
   -- tol_zero : constant double_float := 1.0E-8;
    pp1 : Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_for_Path;
    pp2 : constant Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_End_Game;
    cp1 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_for_Path;
    cp2 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_End_Game;

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

      procedure Call_Path_Tracker is

      -- DESCRIPTION :
      --   Calls the path tracker on the solution ls points to.

        use QuadDobl_Complex_Vector_Norms;
        use QuadDobl_Continuation_Data;
        use QuadDobl_Path_Trackers;
        use Polyhedral_Coefficient_Homotopies;

        function Eval ( x : QuadDobl_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return QuadDobl_Complex_Vectors.Vector is
        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          return Eval(h,dctm(task_id).all,x);
        end Eval;

        function dHt ( x : QuadDobl_Complex_Vectors.Vector;
                       t : Complex_Number )
                     return QuadDobl_Complex_Vectors.Vector is

          res : QuadDobl_Complex_Vectors.Vector(h'range);
          xtl : constant integer32 := x'last+1;

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in wrk'range loop
            res(k) := Eval(j(k,xtl),mf(k,xtl).all,dctm(task_id)(k).all,x);
          end loop;
          return res;
        end dHt;

        function dHx ( x : QuadDobl_Complex_Vectors.Vector;
                       t : Complex_Number )
                     return QuadDobl_Complex_Matrices.Matrix is

          mt : QuadDobl_Complex_Matrices.Matrix(x'range,x'range);

        begin
          Eval(c,REAL_PART(t),dpow(task_id).all,dctm(task_id).all);
          for k in mt'range(1) loop
            for L in mt'range(2) loop
              mt(k,L) := Eval(j(k,L),mf(k,L).all,dctm(task_id)(k).all,x);
            end loop;
          end loop;
          return mt;
        end dHx;

        procedure Track_Path_along_Path is
          new Linear_Single_Normal_Silent_Continue(Max_Norm,Eval,dHt,dHx);

        procedure Track_Path_at_End is
           new Linear_Single_Conditioned_Silent_Continue
                 (Max_Norm,Eval,dHt,dHx);

      begin
        Power_Transform(e,lif,mix,mic.nor.all,dpow(task_id).all);
        Polyhedral_Coefficient_Homotopies.Scale(dpow(task_id).all);
        ls.t := Create(zero);
        declare
          s : Solu_Info := Shallow_Create(ls);
          v : Quad_Double_Vectors.Link_to_Vector;
          e : quad_double := zero;
        begin
          pp1.dist_target := 0.0;
          Track_Path_along_Path(s,t1,tol,false,pp1,cp1);
          Track_Path_at_End(s,t1,tol,false,0,v,e,pp2,cp2);
          ls.err := create(s.cora);
          ls.rco := create(s.rcond);
          ls.res := create(s.resa);
        end;
      end Call_Path_Tracker;

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
        Semaphore.Request(s_cell);   -- take next cell
        if cell_ind = 0 then
          cell_ind := 1;
        else
          cell_ind := cell_ind + 1;
          if not Is_Null(cell_ptr)
           then cell_ptr := Tail_Of(cell_ptr);
          end if;
        end if;
        mycell_ptr := cell_ptr;
        myjob := natural32(cell_ind);
        Semaphore.Release(s_cell);   -- end of first critical section
        exit when Is_Null(mycell_ptr);
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
          Call_Path_Tracker;
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
      for t in 1..nt loop
        dpow(t) := new Standard_Floating_VecVecs.VecVec(e'range);
        dctm(t) := new QuadDobl_Complex_VecVecs.VecVec(c'range);
        for k in dpow(t)'range loop
          dpow(t)(k) := new Standard_Floating_Vectors.Vector(e(k)'range);
          dctm(t)(k) := new QuadDobl_Complex_Vectors.Vector(c(k)'range);
        end loop;
      end loop;
      sols_ptr := sols;
      new_line;
      put_line("launching tasks ...");
      new_line;
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
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector
                 := Standard_Complex_Laur_Functions.Coeff(q(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
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
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector
                 := Standard_Complex_Laur_Functions.Coeff(q(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
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
    c : DoblDobl_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    for i in c'range loop
      declare
        coeff_lq : constant DoblDobl_Complex_Vectors.Vector
                 := DoblDobl_Complex_Laur_Functions.Coeff(q(i));
      begin
        c(i) := new DoblDobl_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
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
    c : DoblDobl_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    for i in c'range loop
      declare
        coeff_lq : constant DoblDobl_Complex_Vectors.Vector
                 := DoblDobl_Complex_Laur_Functions.Coeff(q(i));
      begin
        c(i) := new DoblDobl_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
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
    c : QuadDobl_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    for i in c'range loop
      declare
        coeff_lq : constant QuadDobl_Complex_Vectors.Vector
                 := QuadDobl_Complex_Laur_Functions.Coeff(q(i));
      begin
        c(i) := new QuadDobl_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
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
    c : QuadDobl_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors.Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    mf : Mult_Factors(j'range(1),j'range(2));

  begin
    sup := Supports_of_Polynomial_Systems.Create(q);
    pts := Floating_Integer_Convertors.Convert(sup);
    lif := Floating_Lifting_Utilities.Occurred_Lifting(n,mix,pts,mcc);
    for i in c'range loop
      declare
        coeff_lq : constant QuadDobl_Complex_Vectors.Vector
                 := QuadDobl_Complex_Laur_Functions.Coeff(q(i));
      begin
        c(i) := new QuadDobl_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Exponent_Vectors.Create(q);
    Create(q,j,mf);
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    Floating_Mixed_Subdivisions_io.put(file,lif);
    Reporting_Multitasking_Path_Tracker(q,nt,n,m,mix,lif,mcc,h,c,e,j,mf,sols);
  end Reporting_Multitasking_Path_Tracker;

end Multitasking_Polyhedral_Trackers;
