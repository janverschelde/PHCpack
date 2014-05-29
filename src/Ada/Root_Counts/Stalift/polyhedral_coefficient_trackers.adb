with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Natural_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Supports_of_Polynomial_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Continuation_Data;         use Standard_Continuation_Data;
with Standard_Continuation_Data_io;      use Standard_Continuation_Data_io;
with Standard_Simpomial_Solvers;
with Floating_Lifting_Utilities;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Floating_integer_Convertors;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Random_Coefficient_Systems;
with Polyhedral_Coefficient_Parameters;
with Polyhedral_Coefficient_Homotopies;  use Polyhedral_Coefficient_Homotopies;
with Polyhedral_Coefficient_Predictors;  use Polyhedral_Coefficient_Predictors;
with Polyhedral_Coefficient_Correctors;  use Polyhedral_Coefficient_Correctors;

package body Polyhedral_Coefficient_Trackers is

  procedure Silent_Track_One_Path
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                pow : in Standard_Floating_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                sol : in out Solu_Info; fail : out boolean ) is

    y : Standard_Complex_Vectors.Vector(hq'range);
    tol : constant double_float
        := Polyhedral_Coefficient_Parameters.tol_root;
    max_nwtit : constant natural32
              := Polyhedral_Coefficient_Parameters.max_corr_iter;
    max_nsteps : constant natural32
               := Polyhedral_Coefficient_Parameters.max_nb_steps;
    min_infty : constant double_float
              := Polyhedral_Coefficient_Parameters.min_infinity;
    max_pstep : constant double_float
              := Polyhedral_Coefficient_Parameters.max_pred_step;
    scs : natural32 := 0;
    nit : natural32;
    nwt_fail : boolean := true;
    s,ds,t,dt : double_float;
    ls : Link_to_Solution := sol.sol;
    backup,px : Standard_Complex_Vectors.Vector(ls.v'range);

  begin
    s := -min_infty; ds := max_pstep; t := 0.0;
    backup := ls.v; px := ls.v;
    sol.nstep := 0; sol.nfail := 0; sol.niter := 0; sol.nsyst := 0;
    loop
      sol.nstep := sol.nstep + 1;
      exit when (sol.nstep > max_nsteps);
      Predictor(s,ds);
      dt := t; t := exp(s); dt := t - dt;
      Eval(coeffv,t,pow,ctm);
      if not nwt_fail
       then Secant_Predictor(ls.v,px,dt);
      end if;
      y := Eval(hq,ctm,ls.v);
      Silent_Apply_Newton
        (hq,ctm,jacmat,mulfac,tol,max_nwtit,ls.v,y,nit,nwt_fail);
      sol.niter := sol.niter + nit; sol.nsyst := sol.niter;
      if nwt_fail
       then sol.nfail := sol.nfail + 1;
      end if;
      Step_Control(nwt_fail,s,ds,max_pstep,scs,ls.v,backup,px);
      exit when (s = 0.0) and (not nwt_fail);
    end loop;
    fail := (sol.nstep > max_nsteps);
    Silent_Newton_Step
      (hq,ctm,jacmat,mulfac,ls.v,y,sol.cora,sol.resa,sol.rcond);
    ls.t := Create(t);
    ls.err := sol.cora; ls.rco := sol.rcond; ls.res := sol.resa;
  end Silent_Track_One_Path;

  procedure Reporting_Track_One_Path
              ( file : in file_type;
                hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                pow : in Standard_Floating_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                sol : in out Solu_Info; fail : out boolean ) is

    y : Standard_Complex_Vectors.Vector(hq'range);
    tol : constant double_float
        := Polyhedral_Coefficient_Parameters.tol_root;
    max_nwtit : constant natural32
              := Polyhedral_Coefficient_Parameters.max_corr_iter;
    max_nsteps : constant natural32
               := Polyhedral_Coefficient_Parameters.max_nb_steps;
    min_infty : constant double_float
              := Polyhedral_Coefficient_Parameters.min_infinity;
    max_pstep : constant double_float
              := Polyhedral_Coefficient_Parameters.max_pred_step;
    nit : natural32;
    scs : natural32 := 0;
    nwt_fail : boolean := true;
    s,t,ds,dt : double_float;
    ls : Link_to_Solution := sol.sol;
    backup,px : Standard_Complex_Vectors.Vector(ls.v'range);

  begin
    s := min_infty; ds := max_pstep; t := 0.0;
    backup := ls.v; px := ls.v;
    sol.nstep := 0; sol.nfail := 0; sol.niter := 0; sol.nsyst := 0;
    loop
      sol.nstep := sol.nstep + 1;
      exit when (sol.nstep > max_nsteps);
      Predictor(s,ds);
      put(file,"step "); put(file,sol.nstep,1); put(file," at s =");
      put(file,s,3); put(file," with ds ="); put(file,ds,3);
      put_line(file," :");
      dt := t; t := exp(s); dt := t - dt;
      if not nwt_fail
       then Secant_Predictor(ls.v,px,dt);
      end if;
      Eval(coeffv,t,pow,ctm);
      y := Eval(hq,ctm,ls.v);
      Reporting_Apply_Newton
        (file,hq,ctm,jacmat,mulfac,tol,max_nwtit,ls.v,y,nit,nwt_fail);
      sol.niter := sol.niter + nit; sol.nsyst := sol.niter;
      if nwt_fail
       then sol.nfail := sol.nfail + 1;
      end if;
      Step_Control(nwt_fail,s,ds,max_pstep,scs,ls.v,backup,px);
      exit when (s = 0.0) and (not nwt_fail);
    end loop;
    put(file,"end :");
    Reporting_Newton_Step
      (file,hq,ctm,jacmat,mulfac,ls.v,y,sol.cora,sol.resa,sol.rcond);
    new_line(file);
    ls.t := Create(t);
    ls.err := sol.cora; ls.rco := sol.rcond; ls.res := sol.resa;
    fail := (sol.nstep > max_nsteps);
  end Reporting_Track_One_Path;

  procedure Track_Paths_for_Cell
              ( file : in file_type; report,monitor : in boolean;
                lq : in Laur_Sys;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                ind : in integer32; mic: in Mixed_Cell;
                cnt : in out integer32 ) is

    q : Laur_Sys(lq'range)
      := Supports_of_Polynomial_Systems.Select_Terms(lq,mix,mic.pts.all);
    tol_zero : constant double_float := 1.0E-12;
    fail,zero_y : boolean;
    pow : Standard_Floating_VecVecs.VecVec(coeffv'range)
        := Power_Transform(expvec,ls,mix,mic.nor.all);
    scapow : Standard_Floating_VecVecs.VecVec(coeffv'range) := Scale(pow);
    sols,tmp : Solution_list;
    sol : Link_to_Solution;
    npaths : natural32;
    s : Solu_Info;
    nbfail,nbregu,nbsing,kind : natural32 := 0;

  begin
    Standard_Simpomial_Solvers.Solve(q,tol_zero,sols,fail,zero_y);
    if monitor then
      if cnt > 0 then
        put("Tracked "); put(cnt,1); put(" paths.  ");
      end if;
      npaths := Length_Of(sols);
      put("Cell "); put(ind,1); put(" leads to "); put(npaths,1);
      put_line(" start solutions.");
    end if;
    tmp := sols;
    while not Is_Null(tmp) loop
      sol := Head_Of(tmp);
      s := Shallow_Create(sol);
      if report then
        Reporting_Track_One_Path
          (file,hq,ctm,scapow,coeffv,jacmat,mulfac,s,fail);
      else
        Silent_Track_One_Path(hq,ctm,scapow,coeffv,jacmat,mulfac,s,fail);
      end if;
      if monitor then
        put("#steps : ");  put(s.nstep,3);
        put("  |dx| =");   put(sol.err,2);
        put("  rco =");    put(sol.rco,2);
        put("  |f(x)| ="); put(sol.res,2);
        if fail
         then put_line("  diverged ");
         else put_line("  converged");
        end if;
      end if;
      Write_Next_Solution
        (file,natural32(cnt),s,tol_zero,tol_zero,nbfail,nbregu,nbsing,kind);
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Floating_VecVecs.Clear(pow);
    Standard_Floating_VecVecs.Clear(scapow);
  end Track_Paths_for_Cell;

  procedure Track_Paths_for_Subdivision
              ( file : in file_type; report,monitor : in boolean;
                lq : in Laur_Sys;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                sub : in Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := sub;
    mic : Mixed_Cell;
    ctm : Standard_Complex_VecVecs.VecVec(coeffv'range);
    cnt_cell : integer32 := 0;
    cnt_path : integer32 := 0;

  begin
    for i in coeffv'range loop
      ctm(i) := new Standard_Complex_Vectors.Vector'
                      (coeffv(i).all'range => Create(integer(0)));
    end loop;
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      cnt_cell := cnt_cell + 1;
      Track_Paths_for_Cell
        (file,report,monitor,lq,ls,hq,ctm,coeffv,expvec,jacmat,mulfac,
         mix,cnt_cell,mic,cnt_path);
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Complex_VecVecs.Clear(ctm);
  end Track_Paths_for_Subdivision;

  procedure Track_Paths_for_Subdivision
              ( infile,outfile : in file_type;
                report,monitor : in boolean;
                m : in integer32; lq : in Laur_Sys;
                vs : in Standard_Floating_VecVecs.Array_of_VecVecs;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector ) is

    mlb : Mixed_Labels;
    mic : Mixed_Cell;
    ctm : Standard_Complex_VecVecs.VecVec(coeffv'range);
    cnt_cell : integer32 := 0;
    cnt_path : integer32 := 0;
    fail : boolean;

  begin
    for i in coeffv'range loop
      ctm(i) := new Standard_Complex_Vectors.Vector'
                      (coeffv(i).all'range => Create(integer(0)));
    end loop;
    while cnt_cell < m loop
      cnt_cell := cnt_cell + 1;
      Read_Next(infile,natural32(lq'last),natural32(mix'last),
                natural32(cnt_cell),mlb,fail);
      if fail then
        put("Failed to read cell "); put(cnt_cell,1); put_line("...");
      else
        mic := Create_Coordinates(vs,mlb);
        Track_Paths_for_Cell
          (outfile,report,monitor,lq,ls,hq,ctm,coeffv,expvec,jacmat,mulfac,
           mix,cnt_cell,mic,cnt_path);
        Deep_Clear(mic);
      end if;
      Clear(mlb);
    end loop;
    Standard_Complex_VecVecs.Clear(ctm);
  exception
     when others
       => put("Exception raised when processing cell "); put(cnt_cell,1);
          new_line; return;
  end Track_Paths_for_Subdivision;

  procedure Polyhedral_Continuation
              ( file : in file_type; report,monitor : in boolean;
                n : in integer32; mv : in natural32; p : in Poly_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                sub : in Mixed_Subdivision; q : out Poly_Sys ) is

    s : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
      := Supports_of_Polynomial_Systems.Create(p);
    fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(s'range)
       := Floating_Integer_Convertors.Convert(s);
    ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
       := Floating_Lifting_Utilities.Occurred_Lifting(p'last,mix,fs,sub);
    lq : Laur_Sys(1..n);
    hq : Eval_Coeff_Laur_Sys(1..n);
    expvec : Exponent_Vectors_Array(1..n);
    coeffv : Standard_Complex_VecVecs.VecVec(1..n);
    jacmat : Eval_Coeff_Jaco_Mat(q'range,q'first..q'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));
    timer : Timing_Widget;

  begin
    q := Random_Coefficient_Systems.Create(natural32(n),mix,ls);
    put_line(file,q);
    new_line(file);
    put_line(file,"TITLE : a random coefficient system");
    new_line(file);
    Polyhedral_Coefficient_Parameters.Write(file);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,mv,1); put(file," "); put(file,n,1); new_line(file);
    put(file,"===========================================================");
    new_line(file);
    lq := Polynomial_to_Laurent_System(q);
    hq := Create(lq);
    expvec := Create(q);
    for i in q'range loop
      declare
        c : constant Standard_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        coeffv(i) := new Standard_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    Create(lq,jacmat,mulfac);
   -- Mixed_Solve(file,lq,ls,hq,coeffv,expvec,jacmat,mulfac,mix,sub,qsols);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    tstart(timer);
    Track_Paths_for_Subdivision
      (file,report,monitor,lq,ls,
       hq,coeffv,expvec,jacmat,mulfac,mix,sub);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"polyhedral coefficient path tracking");
  end Polyhedral_Continuation;

  procedure Polyhedral_Continuation
              ( infile,outfile : in file_type;
                report,monitor : in boolean;
                n : in integer32; mv : in natural32;
                mix : in Standard_Integer_Vectors.Vector;
                q : out Poly_Sys ) is

    vs : Standard_Floating_VecVecs.Array_of_VecVecs(mix'range);
    ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    m : integer32 := 0;
    lq : Laur_Sys(1..n);
    hq : Eval_Coeff_Laur_Sys(1..n);
    expvec : Exponent_Vectors_Array(1..n);
    coeffv : Standard_Complex_VecVecs.VecVec(1..n);
    jacmat : Eval_Coeff_Jaco_Mat(q'range,q'first..q'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));
    fail : boolean;
    timer : Timing_Widget;

  begin
    Read_Lifted_Supports(infile,natural32(n),natural32(mix'last),vs,fail);
    if fail then
      put_line("Failed to read the lifted supports.");
    else
      ls := Arrays_of_Floating_Vector_Lists.Shallow_Create(vs);
      q := Random_Coefficient_Systems.Create(natural32(n),mix,ls);
      put_line(outfile,q);
      new_line(outfile);
      put_line(outfile,"TITLE : a random coefficient system");
      new_line(outfile);
      Polyhedral_Coefficient_Parameters.Write(outfile);
      new_line(outfile);
      put_line(outfile,"THE SOLUTIONS :");
      put(outfile,mv,1); put(outfile," ");
      put(outfile,n,1); new_line(outfile);
      put_line(outfile,
        "===========================================================");
      lq := Polynomial_to_Laurent_System(q);
      hq := Create(lq);
      expvec := Create(q);
      for i in q'range loop
        declare
          c : constant Standard_Complex_Vectors.Vector := Coeff(lq(i));
        begin
          coeffv(i) := new Standard_Complex_Vectors.Vector(c'range);
          for j in c'range loop
            coeffv(i)(j) := c(j);
          end loop;
        end;
      end loop;
      Create(lq,jacmat,mulfac);
      get(infile,m);
      put("ready to process "); put(m,1); put_line(" mixed cells...");
      new_line;
      put_line("See the output file for results ...");
      new_line;
      tstart(timer);
      Track_Paths_for_Subdivision
        (infile,outfile,report,monitor,m,lq,vs,ls,
         hq,coeffv,expvec,jacmat,mulfac,mix);
      tstop(timer);
      new_line(outfile);
      print_times(outfile,timer,"polyhedral coefficient path tracking");
    end if;
  end Polyhedral_Continuation;

end Polyhedral_Coefficient_Trackers;
