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
with Standard_Sparse_Solvers;
with Floating_Lifting_Utilities;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Floating_integer_Convertors;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Random_Coefficient_Systems;
with Polyhedral_Coefficient_Homotopies;  use Polyhedral_Coefficient_Homotopies;

package body Polyhedral_Coefficient_Tracking is

  procedure Silent_Newton_Step
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y : out double_float ) is

    jm : Matrix(y'range,x'range);
    nrm : double_float;
    ipvt : Standard_Natural_Vectors.Vector(x'range);
    info : natural;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        jm(i,j) := Eval(jacmat(i,j),mulfac(i,j).all,ctm(i).all,x);
      end loop;
    end loop;
    lufac(jm,x'last,ipvt,info);
    Standard_Complex_Vectors.Min(y);
    lusolve(jm,x'last,ipvt,y);
    nrm := Max_Norm(y); norm_dx := nrm;
    Standard_Complex_Vectors.Add(x,y);
    y := Eval(hq,ctm,x);
    nrm := Max_Norm(y); norm_y := nrm;
  end Silent_Newton_Step;

  procedure Silent_Newton_Step
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y,rcond : out double_float ) is

    jm : Matrix(y'range,x'range);
    nrm : double_float;
    ipvt : Standard_Natural_Vectors.Vector(x'range);

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        jm(i,j) := Eval(jacmat(i,j),mulfac(i,j).all,ctm(i).all,x);
      end loop;
    end loop;
    lufco(jm,x'last,ipvt,rcond);
    Standard_Complex_Vectors.Min(y);
    lusolve(jm,x'last,ipvt,y);
    nrm := Max_Norm(y); norm_dx := nrm;
    Standard_Complex_Vectors.Add(x,y);
    y := Eval(hq,ctm,x);
    nrm := Max_Norm(y); norm_y := nrm;
  end Silent_Newton_Step;

  procedure Reporting_Newton_Step
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y : out double_float ) is

    jm : Matrix(y'range,x'range);
    nrm : double_float;
    ipvt : Standard_Natural_Vectors.Vector(x'range);
    info : natural;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        jm(i,j) := Eval(jacmat(i,j),mulfac(i,j).all,ctm(i).all,x);
      end loop;
    end loop;
    lufac(jm,x'last,ipvt,info);
    Standard_Complex_Vectors.Min(y);
    lusolve(jm,x'last,ipvt,y);
    nrm := Max_Norm(y); norm_dx := nrm;
    put(file,"  |dx| ="); put(file,nrm,2); 
    Standard_Complex_Vectors.Add(x,y);
    y := Eval(hq,ctm,x);
    nrm := Max_Norm(y); norm_y := nrm;
    put(file,"  |f(x)| ="); put(file,nrm,2);
  end Reporting_Newton_Step;

  procedure Reporting_Newton_Step
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y,rcond : out double_float ) is

    jm : Matrix(y'range,x'range);
    nrm : double_float;
    ipvt : Standard_Natural_Vectors.Vector(x'range);

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        jm(i,j) := Eval(jacmat(i,j),mulfac(i,j).all,ctm(i).all,x);
      end loop;
    end loop;
    lufco(jm,x'last,ipvt,rcond);
    Standard_Complex_Vectors.Min(y);
    lusolve(jm,x'last,ipvt,y);
    nrm := Max_Norm(y); norm_dx := nrm;
    put(file,"  |dx| ="); put(file,nrm,2); 
    Standard_Complex_Vectors.Add(x,y);
    y := Eval(hq,ctm,x);
    nrm := Max_Norm(y); norm_y := nrm;
    put(file,"  |f(x)| ="); put(file,nrm,2);
  end Reporting_Newton_Step;

  procedure Silent_Apply_Newton
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                tol : in double_float; max : in natural;
                x,y : in out Standard_Complex_Vectors.Vector;
                nit : out natural; fail : out boolean ) is

    prev_ndx,prev_ny,new_ndx,new_ny : double_float;

  begin
    nit := 1;
    Silent_Newton_Step(hq,ctm,jacmat,mulfac,x,y,prev_ndx,prev_ny);
    fail := (prev_ndx > tol) and (prev_ny > tol);
    while nit < max loop
      nit := nit + 1;
      Silent_Newton_Step(hq,ctm,jacmat,mulfac,x,y,new_ndx,new_ny);
      exit when (new_ndx > prev_ndx) or (new_ny > prev_ny);  -- divergence
      fail := (new_ndx > tol) and (new_ny > tol);
      exit when not fail;
      prev_ndx := new_ndx; prev_ny := new_ny;
    end loop;
  end Silent_Apply_Newton;

  procedure Reporting_Apply_Newton
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                tol : in double_float; max : in natural;
                x,y : in out Standard_Complex_Vectors.Vector;
                nit : out natural; fail : out boolean ) is

    prev_ndx,prev_ny,new_ndx,new_ny : double_float;

  begin
    nit := 1;
    put(file,nit,3); put(file," :");
    Reporting_Newton_Step(file,hq,ctm,jacmat,mulfac,x,y,prev_ndx,prev_ny);
    new_line(file);
    fail := (prev_ndx > tol) and (prev_ny > tol);
    while nit < max loop
      nit := nit + 1;
      put(file,nit,3); put(file," :");
      Reporting_Newton_Step(file,hq,ctm,jacmat,mulfac,x,y,new_ndx,new_ny);
      exit when (new_ndx > prev_ndx) or (new_ny > prev_ny);  -- divergence
      fail := (new_ndx > tol) and (new_ny > tol);
      exit when not fail;
      prev_ndx := new_ndx; prev_ny := new_ny;
      exit when (nit = max);
      new_line(file);
    end loop;
    if fail
     then put_line(file,"  failure");
     else put_line(file,"  success");
    end if;
  end Reporting_Apply_Newton;

  procedure Silent_Track_One_Path
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                pow : in Standard_Floating_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                sol : in out Solution;
                nsteps : out natural; fail : out boolean ) is

    y : Standard_Complex_Vectors.Vector(hq'range);
    tol : constant double_float := 1.0E-12;
    max : constant natural := 4;
    nit : natural;
    nwt_fail : boolean;
    s,t,ds,nrm_dx,nrm_y,rcond : double_float;
    backup : Standard_Complex_Vectors.Vector(sol.v'range);
    step : natural := 0;
    max_steps : constant natural := 500;

  begin
    s := -5.0;
    ds := 0.1;
    backup := sol.v;
    loop
      step := step + 1;
      exit when (step > max_steps);
      s := s + ds;
      if s > 0.0
       then ds := s - ds;  -- old value for s
            s := 0.0;
            ds := abs(ds); -- do not overshoot again
      end if;
      t := exp(s);
      Eval(coeffv,t,pow,ctm);
      y := Eval(hq,ctm,sol.v);
      Silent_Apply_Newton
        (hq,ctm,jacmat,mulfac,tol,max,sol.v,y,nit,nwt_fail);
      if not nwt_fail then
        if (ds < 0.1) and (ds < -s)
         then ds := 2.0*ds;
        end if;
        backup := sol.v;
      else
        s := s - ds;
        ds := ds/2.0;
        sol.v := backup;
      end if;
      exit when (s = 0.0) and (not nwt_fail);
    end loop;
    fail := (step > max_steps);
    nsteps := step;
    Silent_Newton_Step(hq,ctm,jacmat,mulfac,sol.v,y,nrm_dx,nrm_y,rcond);
    sol.err := nrm_dx; sol.rco := rcond; sol.res := nrm_y;
  end Silent_Track_One_Path;

  procedure Reporting_Track_One_Path
              ( file : in file_type;
                hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                pow : in Standard_Floating_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                sol : in out Solution;
                nsteps : out natural; fail : out boolean ) is

    y : Standard_Complex_Vectors.Vector(hq'range);
    nrm : double_float;
    tol : constant double_float := 1.0E-12;
    max : constant natural := 4;
    nit : natural;
    nwt_fail : boolean;
    s,t,ds,nrm_dx,nrm_y,rcond : double_float;
    backup : Standard_Complex_Vectors.Vector(sol.v'range);
    step : natural := 0;
    max_steps : constant natural := 500;

  begin
    Eval(coeffv,EXP(-10.0),pow,ctm);
    y := Eval(hq,ctm,sol.v);
    nrm := Max_Norm(y);
    put(file,"Residual at the start for t = exp(-10) : ");
    put(file,nrm); new_line(file);
    s := -5.0;
    ds := 0.1;
    backup := sol.v;
    loop
      step := step + 1;
      exit when (step > max_steps);
      s := s + ds;
      if s > 0.0
       then ds := s - ds;  -- old value for s
            s := 0.0;
            ds := abs(ds); -- do not overshoot again
      end if;
      put(file,"step "); put(file,step,1); put(file," at s =");
      put(file,s,3); put(file," with ds ="); put(file,ds,3);
      put_line(file," :");
      t := exp(s);
      Eval(coeffv,t,pow,ctm);
      y := Eval(hq,ctm,sol.v);
      Reporting_Apply_Newton
        (file,hq,ctm,jacmat,mulfac,tol,max,sol.v,y,nit,nwt_fail);
      if not nwt_fail then
        if (ds < 0.1) and (ds < -s)
         then ds := 2.0*ds;
        end if;
        backup := sol.v;
      else
        s := s - ds;
        ds := ds/2.0;
        sol.v := backup;
      end if;
      exit when (s = 0.0) and (not nwt_fail);
    end loop;
    put(file,"end :");
    Reporting_Newton_Step
      (file,hq,ctm,jacmat,mulfac,sol.v,y,nrm_dx,nrm_y,rcond);
    new_line(file);
    sol.err := nrm_dx; sol.rco := rcond; sol.res := nrm_y;
    if step > max_steps
     then fail := true;
     else fail := false;
    end if;
    nsteps := step;
  end Reporting_Track_One_Path;

  procedure Track_Paths_for_Cell
              ( file : in file_type; lq : in Laur_Sys;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                ind : in natural; mic: in Mixed_Cell;
                qsols,last : in out Solution_List ) is

    q : Laur_Sys(lq'range)
      := Supports_of_Polynomial_Systems.Select_Terms(lq,mix,mic.pts.all);
    tol_zero : constant double_float := 1.0E-12;
    fail,zero_y : boolean;
    pow : Standard_Floating_VecVecs.VecVec(coeffv'range)
        := Power_Transform(expvec,ls,mix,mic.nor.all);
    scapow : Standard_Floating_VecVecs.VecVec(coeffv'range) := Scale(pow);
    sols,tmp : Solution_list;
    sol : Link_to_Solution;
    nsteps : natural;

  begin
    Standard_Sparse_Solvers.Solve(q,tol_zero,sols,fail,zero_y);
    put(file,"Mixed cell "); put(file,ind,1); put(file," leads to ");
    put(file,Length_Of(sols),1);
    put_line(file," start solutions.");
    tmp := sols;
    while not Is_Null(tmp) loop
      sol := Head_Of(tmp);
     -- Reporting_Track_One_Path
     --   (file,hq,ctm,scapow,coeffv,jacmat,mulfac,sol.all,nsteps,fail);
      Silent_Track_One_Path
        (hq,ctm,scapow,coeffv,jacmat,mulfac,sol.all,nsteps,fail);
      put(file,"#steps : "); put(file,nsteps,3);
      put(file,"  |dx| ="); put(file,sol.err,2);
      put(file,"  rco ="); put(file,sol.rco,3);
      put(file,"  |f(x)| ="); put(file,sol.res,2);
      if fail
       then put_line(file,"  diverged ");
       else put_line(file,"  converged");
      end if;
      Append(qsols,last,sol.all);
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Floating_VecVecs.Clear(pow);
    Standard_Floating_VecVecs.Clear(scapow);
  end Track_Paths_for_Cell;

  procedure Track_Paths_for_Subdivision
              ( file : in file_type; lq : in Laur_Sys;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                sub : in Mixed_Subdivision; qsols : out Solution_List ) is

    tmp : Mixed_Subdivision := sub;
    mic : Mixed_Cell;
    ctm : Standard_Complex_VecVecs.VecVec(coeffv'range);
    cnt : natural := 0;
    last : Solution_List := qsols;

  begin
    for i in coeffv'range loop
      ctm(i) := new Standard_Complex_Vectors.Vector'
                      (coeffv(i).all'range => Create(0));
    end loop;
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      cnt := cnt + 1;
      Track_Paths_for_Cell
        (file,lq,ls,hq,ctm,coeffv,expvec,jacmat,mulfac,mix,cnt,mic,
         qsols,last);
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Complex_VecVecs.Clear(ctm);
  end Track_Paths_for_Subdivision;

  procedure Track_Paths_for_Subdivision
              ( infile,outfile : in file_type; m : in natural;
                lq : in Laur_Sys;
                vs : in Standard_Floating_VecVecs.Array_of_VecVecs;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                qsols : out Solution_List ) is

    mlb : Mixed_Labels;
    mic : Mixed_Cell;
    ctm : Standard_Complex_VecVecs.VecVec(coeffv'range);
    cnt : natural := 0;
    last : Solution_List := qsols;
    fail : boolean;

  begin
    for i in coeffv'range loop
      ctm(i) := new Standard_Complex_Vectors.Vector'
                      (coeffv(i).all'range => Create(0));
    end loop;
    while cnt < m loop
      cnt := cnt + 1;
      Read_Next(infile,lq'last,mix'last,cnt,mlb,fail);
      if fail then
        put("Failed to read cell "); put(cnt,1); put_line("...");
      else
        mic := Create_Coordinates(vs,mlb);
        Track_Paths_for_Cell
          (outfile,lq,ls,hq,ctm,coeffv,expvec,jacmat,mulfac,mix,cnt,mic,
           qsols,last);
        Deep_Clear(mic);
      end if;
      Clear(mlb);
    end loop;
    Standard_Complex_VecVecs.Clear(ctm);
  exception
     when others
       => put("Exception raised when processing cell "); put(cnt,1);
          new_line; return;
  end Track_Paths_for_Subdivision;

  procedure Polyhedral_Continuation
              ( file : in file_type; n : in natural; p : in Poly_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                sub : in Mixed_Subdivision;
                q : out Poly_Sys; qsols : out Solution_List ) is

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

  begin
    q := Random_Coefficient_Systems.Create(n,mix,ls);
    put_line(file,q);
    new_line(file);
    put_line(file,"TITLE : a random coefficient system");
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
    Track_Paths_for_Subdivision
      (file,lq,ls,hq,coeffv,expvec,jacmat,mulfac,mix,sub,qsols);
  end Polyhedral_Continuation;

  procedure Polyhedral_Continuation
              ( infile,outfile : in file_type; 
                n : in natural; p : in Poly_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                q : out Poly_Sys; qsols : out Solution_List ) is

    vs : Standard_Floating_VecVecs.Array_of_VecVecs(mix'range);
    ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    m : natural;
    lq : Laur_Sys(1..n);
    hq : Eval_Coeff_Laur_Sys(1..n);
    expvec : Exponent_Vectors_Array(1..n);
    coeffv : Standard_Complex_VecVecs.VecVec(1..n);
    jacmat : Eval_Coeff_Jaco_Mat(q'range,q'first..q'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));
    fail : boolean;

  begin
    Read_Lifted_Supports(infile,n,mix'last,vs,fail);
    if fail then
      put_line("Failed to read the lifted supports.");
    else
      ls := Arrays_of_Floating_Vector_Lists.Shallow_Create(vs);
      q := Random_Coefficient_Systems.Create(n,mix,ls);
      put_line(outfile,q);
      new_line(outfile);
      put_line(outfile,"TITLE : a random coefficient system");
      new_line(outfile);
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
      Track_Paths_for_Subdivision
        (infile,outfile,m,lq,vs,ls,hq,coeffv,expvec,jacmat,mulfac,mix,qsols);
    end if;
  end Polyhedral_Continuation;

end Polyhedral_Coefficient_Tracking;
