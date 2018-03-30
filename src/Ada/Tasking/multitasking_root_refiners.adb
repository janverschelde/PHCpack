with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Linear_Solvers;   use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Linear_Solvers;   use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Linear_Solvers;   use QuadDobl_Complex_Linear_Solvers;
with Generic_Lists;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Laur_JacoMats;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_JacoMats;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;
with Standard_Solution_Diagnostics;
with Standard_Condition_Tables;
with DoblDobl_Solution_Diagnostics;
with DoblDobl_Condition_Tables;
with QuadDobl_Solution_Diagnostics;
with QuadDobl_Condition_Tables;
with Standard_Root_Refiners;
with DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;
with Multitasking;

package body Multitasking_Root_Refiners is

-- INTERNAL DATA STRUCTURE :

  type Solution_Metric is record
    label : integer32;        -- label of solution
    initres : double_float;   -- initial residual
    fail : boolean;           -- failed to converge or not
    numbits : natural32;      -- number of iterations
    infty : boolean;          -- at infinity or not
  end record;
  type Link_to_Solution_Metric is access Solution_Metric;

  package List_of_Solution_Metrics is
    new Generic_Lists(Link_to_Solution_Metric);
  type Solution_Metric_List is new List_of_Solution_Metrics.List;

  function Create ( n : integer32 ) return Solution_Metric_List is

  -- DESCRIPTION :
  --   Returns a list of n solution metrics.

    res,res_last : Solution_Metric_List;

  begin
    for i in 1..n loop
      declare
        lsm : Link_to_Solution_Metric;
      begin
        lsm := new Solution_Metric;
        lsm.label := i;
        lsm.initres := 1.0;
        lsm.fail := true;
        lsm.numbits := 0;
        lsm.infty := false;
        Append(res,res_last,lsm);
      end;
    end loop;
    return res;
  end Create;

-- MUTE ROOT REFINING PROCEDURES ON POLYNOMIAL SYSTEMS :

  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Solution_Diagnostics;
    use Standard_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : Standard_Complex_Vectors.Vector(p'range);
      A : Standard_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := Sum_Norm(y);
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              Standard_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + 1.0 = 1.0 then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              Standard_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution.

    begin
      do_jobs(nt);
      Clear(f); Clear(jm); Clear(jf);
    end Main;

  begin
    Main;
  end Mute_Multitasking_Root_Refiner;

  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Diagnostics;
    use DoblDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant double_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : DoblDobl_Complex_Vectors.Vector(p'range);
      A : DoblDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              DoblDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              DoblDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution.

    begin
      do_jobs(nt);
      Clear(f); Clear(jm); Clear(jf);
    end Main;

  begin
    Main;
  end Mute_Multitasking_Root_Refiner;

  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Diagnostics;
    use QuadDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant quad_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : QuadDobl_Complex_Vectors.Vector(p'range);
      A : QuadDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hihi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              QuadDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              QuadDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution.

    begin
      do_jobs(nt);
      Clear(f); Clear(jm); Clear(jf);
    end Main;

  begin
    Main;
  end Mute_Multitasking_Root_Refiner;

-- SILENT ROOT REFINING PROCEDURES ON POLYNOMIAL SYSTEMS :

  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Solution_Diagnostics;
    use Standard_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : Standard_Complex_Vectors.Vector(p'range);
      A : Standard_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := Sum_Norm(y);
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              Standard_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + 1.0 = 1.0 then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              Standard_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                        := Standard_Condition_Tables.Create(15); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;

    begin
      do_jobs(nt);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        Write_Info(file,ls.all,lsm.initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        Standard_Condition_Tables.Update_Corrector(t_err,ls.err);
        Standard_Condition_Tables.Update_Condition(t_rco,ls.rco);
        Standard_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Silent_Multitasking_Root_Refiner;

  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Diagnostics;
    use DoblDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant double_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : DoblDobl_Complex_Vectors.Vector(p'range);
      A : DoblDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              DoblDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              DoblDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..30)
                        := DoblDobl_Condition_Tables.Create(30); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;
      dd_initres : double_double;

    begin
      do_jobs(nt);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        dd_initres := create(lsm.initres);
        Write_Info(file,ls.all,dd_initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        DoblDobl_Condition_Tables.Update_Corrector(t_err,ls.err);
        DoblDobl_Condition_Tables.Update_Condition(t_rco,ls.rco);
        DoblDobl_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      DoblDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Silent_Multitasking_Root_Refiner;

  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Diagnostics;
    use QuadDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant quad_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : QuadDobl_Complex_Vectors.Vector(p'range);
      A : QuadDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hihi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              QuadDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              QuadDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                        := QuadDobl_Condition_Tables.Create(60); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;
      qd_initres : quad_double;

    begin
      do_jobs(nt);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        qd_initres := create(lsm.initres);
        Write_Info(file,ls.all,qd_initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.err);
        QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.rco);
        QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Silent_Multitasking_Root_Refiner;

  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_Solution_Diagnostics;
    use Standard_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : Standard_Complex_Vectors.Vector(p'range);
      A : Standard_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          put_line("task " & Multitasking.to_string(natural32(task_id))
                           & " refines solution "
                           & Multitasking.to_string(natural32(i)));
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := Sum_Norm(y);
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              Standard_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + 1.0 = 1.0 then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              Standard_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                        := Standard_Condition_Tables.Create(15); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;

    begin
      new_line;
      put("launching "); put(nt,1); put_line(" tasks ...");
      new_line;
      do_jobs(nt);
      new_line;
      put_line("writing results to file ...");
      new_line;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        Write_Info(file,ls.all,lsm.initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        Standard_Condition_Tables.Update_Corrector(t_err,ls.err);
        Standard_Condition_Tables.Update_Condition(t_rco,ls.rco);
        Standard_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Reporting_Multitasking_Root_Refiner;

  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Diagnostics;
    use DoblDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant double_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : DoblDobl_Complex_Vectors.Vector(p'range);
      A : DoblDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          put_line("task " & Multitasking.to_string(natural32(task_id))
                           & " refines solution "
                           & Multitasking.to_string(natural32(i)));
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              DoblDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              DoblDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..30)
                        := Standard_Condition_Tables.Create(30); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;
      dd_initres : double_double;

    begin
      new_line;
      put("launching "); put(nt,1); put_line(" tasks ...");
      new_line;
      do_jobs(nt);
      new_line;
      put_line("writing results to file ...");
      new_line;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        dd_initres := create(lsm.initres);
        Write_Info(file,ls.all,dd_initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        DoblDobl_Condition_Tables.Update_Corrector(t_err,ls.err);
        DoblDobl_Condition_Tables.Update_Condition(t_rco,ls.rco);
        DoblDobl_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      DoblDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Reporting_Multitasking_Root_Refiner;

  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Diagnostics;
    use QuadDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant quad_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : QuadDobl_Complex_Vectors.Vector(p'range);
      A : QuadDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          put_line("task " & Multitasking.to_string(natural32(task_id))
                           & " refines solution "
                           & Multitasking.to_string(natural32(i)));
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hihi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              QuadDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              QuadDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                        := Standard_Condition_Tables.Create(60); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;
      qd_initres : quad_double;

    begin
      new_line;
      put("launching "); put(nt,1); put_line(" tasks ...");
      new_line;
      do_jobs(nt);
      new_line;
      put_line("writing results to file ...");
      new_line;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        qd_initres := create(lsm.initres);
        Write_Info(file,ls.all,qd_initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.err);
        QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.rco);
        QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Reporting_Multitasking_Root_Refiner;

-- MUTE ROOT REFINERS ON LAURENT SYSTEMS :

  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;
    use Standard_Solution_Diagnostics;
    use Standard_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Laur_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : Standard_Complex_Vectors.Vector(p'range);
      A : Standard_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := Sum_Norm(y);
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              Standard_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + 1.0 = 1.0 then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              Standard_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution.

    begin
      do_jobs(nt);
      Clear(f); Clear(jm); Clear(jf);
    end Main;

  begin
    Main;
  end Mute_Multitasking_Root_Refiner;

  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Diagnostics;
    use DoblDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Laur_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant double_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : DoblDobl_Complex_Vectors.Vector(p'range);
      A : DoblDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              DoblDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              DoblDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution.

    begin
      do_jobs(nt);
      Clear(f); Clear(jm); Clear(jf);
    end Main;

  begin
    Main;
  end Mute_Multitasking_Root_Refiner;

  procedure Mute_Multitasking_Root_Refiner
              ( nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Diagnostics;
    use QuadDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Laur_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant quad_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : QuadDobl_Complex_Vectors.Vector(p'range);
      A : QuadDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hihi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              QuadDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              QuadDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution.

    begin
      do_jobs(nt);
      Clear(f); Clear(jm); Clear(jf);
    end Main;

  begin
    Main;
  end Mute_Multitasking_Root_Refiner;

-- SILENT ROOT REFINERS ON LAURENT SYSTEMS :

  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;
    use Standard_Solution_Diagnostics;
    use Standard_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Laur_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : Standard_Complex_Vectors.Vector(p'range);
      A : Standard_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := Sum_Norm(y);
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              Standard_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + 1.0 = 1.0 then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              Standard_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                        := Standard_Condition_Tables.Create(15); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;

    begin
      do_jobs(nt);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        Write_Info(file,ls.all,lsm.initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        Standard_Condition_Tables.Update_Corrector(t_err,ls.err);
        Standard_Condition_Tables.Update_Condition(t_rco,ls.rco);
        Standard_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Silent_Multitasking_Root_Refiner;

  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Diagnostics;
    use DoblDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Laur_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant double_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : DoblDobl_Complex_Vectors.Vector(p'range);
      A : DoblDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              DoblDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              DoblDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..30)
                        := Standard_Condition_Tables.Create(30); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;
      dd_initres : double_double;

    begin
      do_jobs(nt);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        dd_initres := create(lsm.initres);
        Write_Info(file,ls.all,dd_initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        DoblDobl_Condition_Tables.Update_Corrector(t_err,ls.err);
        DoblDobl_Condition_Tables.Update_Condition(t_rco,ls.rco);
        DoblDobl_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      DoblDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Silent_Multitasking_Root_Refiner;

  procedure Silent_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Diagnostics;
    use QuadDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Laur_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : constant Solution_Metric_List := Create(integer32(len));
    one : constant quad_double := create(1.0);

    procedure Job ( task_id,nb_tasks : in integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : QuadDobl_Complex_Vectors.Vector(p'range);
      A : QuadDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hihi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              QuadDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              QuadDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                        := Standard_Condition_Tables.Create(60); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;
      qd_initres : quad_double;

    begin
      do_jobs(nt);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        qd_initres := create(lsm.initres);
        Write_Info(file,ls.all,qd_initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.err);
        QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.rco);
        QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Silent_Multitasking_Root_Refiner;

  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;
    use Standard_Solution_Diagnostics;
    use Standard_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Laur_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : Solution_Metric_List := Create(integer32(len));

    procedure Job ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : Standard_Complex_Vectors.Vector(p'range);
      A : Standard_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          put_line("task " & Multitasking.to_string(natural32(task_id))
                           & " refines solution "
                           & Multitasking.to_string(natural32(i)));
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := Sum_Norm(y);
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              Standard_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + 1.0 = 1.0 then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              Standard_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                        := Standard_Condition_Tables.Create(15); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;

    begin
      new_line;
      put("launching "); put(nt,1); put_line(" tasks ...");
      new_line;
      do_jobs(nt);
      new_line;
      put_line("writing results to file ...");
      new_line;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        Write_Info(file,ls.all,lsm.initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        Standard_Condition_Tables.Update_Corrector(t_err,ls.err);
        Standard_Condition_Tables.Update_Condition(t_rco,ls.rco);
        Standard_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Reporting_Multitasking_Root_Refiner;

  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Diagnostics;
    use DoblDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Laur_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : Solution_Metric_List := Create(integer32(len));
    one : constant double_double := create(1.0);

    procedure Job ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : DoblDobl_Complex_Vectors.Vector(p'range);
      A : DoblDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          put_line("task " & Multitasking.to_string(natural32(task_id))
                           & " refines solution "
                           & Multitasking.to_string(natural32(i)));
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              DoblDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              DoblDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..30)
                        := DoblDobl_Condition_Tables.Create(30); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;
      dd_initres : double_double;

    begin
      new_line;
      put("launching "); put(nt,1); put_line(" tasks ...");
      new_line;
      do_jobs(nt);
      new_line;
      put_line("writing results to file ...");
      new_line;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        dd_initres := create(lsm.initres);
        Write_Info(file,ls.all,dd_initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        DoblDobl_Condition_Tables.Update_Corrector(t_err,ls.err);
        DoblDobl_Condition_Tables.Update_Condition(t_rco,ls.rco);
        DoblDobl_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      DoblDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Reporting_Multitasking_Root_Refiner;

  procedure Reporting_Multitasking_Root_Refiner
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                epsxa,epsfa,tolsing : in double_float;
                numit : in out natural32; max : in natural32;
                deflate : in out boolean ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Diagnostics;
    use QuadDobl_Root_Refiners;

    n : constant integer32 := p'length;
    f : Eval_Laur_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    len : constant natural32 := Length_Of(sols);
    slm : Solution_Metric_List := Create(integer32(len));
    one : constant quad_double := create(1.0);

    procedure Job ( task_id,nb_tasks : integer32 ) is

    -- DESCRIPTION :
    --   Defines the job for the task with identification task_id,
    --   out of a total of tasks equal to the number nb_tasks.
    --   Every task does one Newton step on a solution,
    --   depending whether the index of the solution matches
    --   the job identification number, modulo the number of tasks.

      y : QuadDobl_Complex_Vectors.Vector(p'range);
      A : QuadDobl_Complex_Matrices.Matrix(p'range,p'range);
      ipvt : Standard_Integer_Vectors.Vector(1..n);
      index : integer32;
      ptr_sols : Solution_List := sols;
      ptr_slms : Solution_Metric_List := slm;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(task_id)));
      for i in 1..integer32(len) loop
        index := 1 + i mod nb_tasks;
        if index = task_id then
          put_line("task " & Multitasking.to_string(natural32(task_id))
                           & " refines solution "
                           & Multitasking.to_string(natural32(i)));
          ls := Head_Of(ptr_sols);
          lsm := Head_Of(ptr_slms);
          y := Eval(f,ls.v);
          lsm.initres := hihi_part(Sum_Norm(y));
          lsm.infty := At_Infinity(ls.all,false,1.0E+8);
          if not lsm.infty and ls.err < 0.1 and lsm.initres < 0.1 then
            for k in 1..max loop
              QuadDobl_Complex_Vectors.Min(y);
              A := Eval(jf,ls.v);
              lufco(A,n,ipvt,ls.rco);
              if ls.rco + one = one then
                ls.res := Sum_Norm(y);
                lsm.fail := (ls.res > epsfa); exit;
              end if;
              lusolve(A,n,ipvt,y);
              QuadDobl_Complex_Vectors.Add(ls.v,y);
              ls.err := Sum_Norm(y);
              y := Eval(f,ls.v);           
              ls.res := Sum_Norm(y);
              lsm.numbits := lsm.numbits + 1;
              if ls.err < epsxa then
                lsm.fail := false; exit;
              elsif ls.res < epsfa then
                lsm.fail := false; exit;
              end if;
            end loop;
           -- Multiplicity(ls,i,sols,lsm.fail,lsm.infty,deflate,tolsing,epsxa);
          end if;
        end if;
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

    procedure Main is

    -- DESCRIPTION :
    --   Launches a number of tasks to do one Newton step on each solution
    --   and then writes the solution characteristics and the condition
    --   tables to file.

      ptr_sols : Solution_List;
      ptr_slms : Solution_Metric_List;
      ls : Link_to_Solution;
      lsm : Link_to_Solution_Metric;
      t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                        := QuadDobl_Condition_Tables.Create(60); 
      nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing : natural32 := 0;
      qd_initres : quad_double;

    begin
      new_line;
      put("launching "); put(nt,1); put_line(" tasks ...");
      new_line;
      do_jobs(nt);
      new_line;
      put_line("writing results to file ...");
      new_line;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,len,1); put(file," ");
      put(file,n,1); new_line(file);
      Standard_Complex_Solutions_io.put_bar(file);
      ptr_sols := sols;
      ptr_slms := slm;
      for i in 1..len loop
        ls := Head_Of(ptr_sols);
        lsm := Head_Of(ptr_slms);
        qd_initres := create(lsm.initres);
        Write_Info(file,ls.all,qd_initres,i,lsm.numbits,0,lsm.fail,lsm.infty);
        Write_Type(file,ls,lsm.fail,lsm.infty,tolsing,nbfail,nbinfty,
                   nbreal,nbcomp,nbreg,nbsing);
        QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.err);
        QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.rco);
        QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.res);
        ptr_sols := Tail_Of(ptr_sols);
        ptr_slms := Tail_Of(ptr_slms);
      end loop;
      Clear(f); Clear(jm); Clear(jf);
      Write_Global_Info
        (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,0);
      QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    end Main;

  begin
    Main;
  end Reporting_Multitasking_Root_Refiner;

end Multitasking_Root_Refiners;
