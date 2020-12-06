with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_JacoMats;
with DoblDobl_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_JacoMats;
with QuadDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with Standard_Poly_Laur_Convertors;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Poly_Laur_Convertors;
with QuadDobl_Poly_Laur_Convertors;
with Standard_Complex_Poly_Randomizers; 
with Standard_Complex_Laur_Randomizers;  
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Continuation_Parameters;
with Exponent_Vectors;                   use Exponent_Vectors;
with Random_Coefficient_Systems;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Integer_Polyhedral_Continuation;    use Integer_Polyhedral_Continuation;
with Floating_Polyhedral_Continuation;   use Floating_Polyhedral_Continuation;
with DoblDobl_Polyhedral_Continuation;   use DoblDobl_Polyhedral_Continuation;
with QuadDobl_Polyhedral_Continuation;   use QuadDobl_Polyhedral_Continuation;
with Stable_Polyhedral_Continuation;     use Stable_Polyhedral_Continuation;
with Black_Mixed_Volume_Computations;    use Black_Mixed_Volume_Computations;
with Multitasking_Polyhedral_Trackers;   use Multitasking_Polyhedral_Trackers;

package body Black_Polyhedral_Continuations is

  procedure Black_Box_Polyhedral_Continuation
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mix : in Standard_Integer_Vectors.Vector;
                 lifsup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Laur_Functions;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Standard_Laur_Poly_Convertors;
    use Standard_Poly_Laur_Convertors;
    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;

    n : constant integer32 := p'length;
    lq,llq : Laur_Sys(p'range);
    h : Eval_Coeff_Laur_Sys(q'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));

  begin
    if verbose > 0 then
      put("-> in black_polyhedral_continuations.");
      put_line("Black_Box_Polyhedral_Continuation 1,");
      put_line("for a polynomial system in double precision ...");
    end if;
    q := Standard_Complex_Poly_Randomizers.Complex_Randomize1(p);
    lq := Polynomial_to_Laurent_System(q);
    llq := Perform_Lifting(n,mix,lifsup,lq);
    Clear(lq); Clear(q);
    lq := Eval(llq,Create(1.0),n+1);
    q := Laurent_to_Polynomial_System(lq);
   -- Mixed_Solve(llq,mix,mixsub,qsols);   too expensive !!!!
    h := Create(lq);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Create(lq);
    Create(lq,j,m);
    Continuation_Parameters.Tune(0);
    Continuation_Parameters.start_end_game := 0.0;
    Mixed_Solve(llq,lifsup,h,c,e,j,m,mix,mixsub,qsols);
    Set_Continuation_Parameter(qsols,Create(0.0));
    Continuation_Parameters.start_end_game := 0.1;
    Clear(lq); Clear(llq);
    Clear(h); Clear(j); Clear(m);
    Standard_Complex_VecVecs.Clear(c);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 mix : in Standard_Integer_Vectors.Vector;
                 lifsup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Laur_Functions;
    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;

    n : constant integer32 := p'length;
    lq,llq : Laur_Sys(p'range);
    h : Eval_Coeff_Laur_Sys(q'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));

  begin
    if verbose > 0 then
      put("-> in black_polyhedral_continuations.");
      put_line("Black_Box_Polyhedral_Continuation 2,");
      put_line("for a Laurent system in double precision ...");
    end if;
    lq := Standard_Complex_Laur_Randomizers.Complex_Randomize1(p);
    llq := Perform_Lifting(n,mix,lifsup,lq);
    Clear(lq); Clear(q);
    q := Eval(llq,Create(1.0),n+1);
   -- Mixed_Solve(llq,mix,mixsub,qsols);   too expensive !!!!
    h := Create(q);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector := Coeff(q(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Create(q);
    Create(q,j,m);
    Continuation_Parameters.Tune(0);
    Continuation_Parameters.start_end_game := 0.0;
    Mixed_Solve(llq,lifsup,h,c,e,j,m,mix,mixsub,qsols);
    Set_Continuation_Parameter(qsols,Create(0.0));
    Continuation_Parameters.start_end_game := 0.1;
    Clear(llq); Clear(h); Clear(j); Clear(m);
    Standard_Complex_VecVecs.Clear(c);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Laur_Functions;
    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;

    hq : Eval_Coeff_Laur_Sys(p'range);
    expvec : Exponent_Vectors_Array(p'range);
    coeffv : Standard_Complex_VecVecs.VecVec(p'range);
    jacmat : Eval_Coeff_Jaco_Mat(p'range,p'first..p'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));

  begin
    if verbose > 0 then
      put("-> in black_polyhedral_continuations.");
      put_line("Black_Box_Polyhedral_Continuation 3");
      put_line("for a Laurent system in double precision ...");
    end if;
    q := Random_Coefficient_Systems.Create(natural32(p'last),mix.all,lifsup);
    hq := Create(q);
    expvec := Create(q);
    for i in q'range loop
      declare
        c : constant Standard_Complex_Vectors.Vector := Coeff(q(i));
      begin
        coeffv(i) := new Standard_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    Create(q,jacmat,mulfac);
    Continuation_Parameters.Tune(0);
    Continuation_Parameters.start_end_game := 0.0;
    if nt = 0 then
      Mixed_Solve(q,lifsup,hq,coeffv,expvec,jacmat,mulfac,
                  mix.all,mcc,qsols);
    else
      Silent_Multitasking_Path_Tracker
        (q,nt,p'last,mix'last,mix.all,lifsup,mcc,
         hq,coeffv,expvec,jacmat,mulfac,qsols);
    end if;
    Set_Continuation_Parameter(qsols,Create(0.0));
    Continuation_Parameters.start_end_game := 0.1;
    Clear(hq); Clear(jacmat); Clear(mulfac);
    Standard_Complex_VecVecs.Clear(coeffv);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out DoblDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Laur_Functions;
    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;

    hq : Eval_Coeff_Laur_Sys(p'range);
    expvec : Exponent_Vectors_Array(p'range);
    coeffv : DoblDobl_Complex_VecVecs.VecVec(p'range);
    jacmat : Eval_Coeff_Jaco_Mat(p'range,p'first..p'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));
    zero : constant double_double := create(0.0);

  begin
    if verbose > 0 then
      put("-> in black_polyhedral_continuations.");
      put_line("Black_Box_Polyhedral_Continuation 4,");
      put_line("for a Laurent system in double double precision ...");
    end if;
    q := Random_Coefficient_Systems.Create(natural32(p'last),mix.all,lifsup);
    hq := Create(q);
    expvec := Create(q);
    for i in q'range loop
      declare
        c : constant DoblDobl_Complex_Vectors.Vector := Coeff(q(i));
      begin
        coeffv(i) := new DoblDobl_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    Create(q,jacmat,mulfac);
    Continuation_Parameters.Tune(0); --,32);
    Continuation_Parameters.start_end_game := 0.0;
    if nt = 0 then
      Mixed_Solve(q,lifsup,hq,coeffv,expvec,jacmat,mulfac,
                  mix.all,mcc,qsols);
    else
      Silent_Multitasking_Path_Tracker
        (q,nt,p'last,mix'last,mix.all,lifsup,mcc,
         hq,coeffv,expvec,jacmat,mulfac,qsols);
    end if;
    Set_Continuation_Parameter(qsols,Create(zero));
    Continuation_Parameters.start_end_game := 0.1;
    Clear(hq); Clear(jacmat); Clear(mulfac);
    DoblDobl_Complex_VecVecs.Clear(coeffv);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out QuadDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Laur_Functions;
    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;

    hq : Eval_Coeff_Laur_Sys(p'range);
    expvec : Exponent_Vectors_Array(p'range);
    coeffv : QuadDobl_Complex_VecVecs.VecVec(p'range);
    jacmat : Eval_Coeff_Jaco_Mat(p'range,p'first..p'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));
    zero : constant quad_double := create(0.0);

  begin
    if verbose > 0 then
      put("-> in black_polyhedral_continuations.");
      put_line("Black_Box_Polyhedral_Continuation 5");
      put_line("for a Laurent system in quad double precision ...");
    end if;
    q := Random_Coefficient_Systems.Create(natural32(p'last),mix.all,lifsup);
    hq := Create(q);
    expvec := Create(q);
    for i in q'range loop
      declare
        c : constant QuadDobl_Complex_Vectors.Vector := Coeff(q(i));
      begin
        coeffv(i) := new QuadDobl_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    Create(q,jacmat,mulfac);
    Continuation_Parameters.Tune(0); --,64);
    Continuation_Parameters.start_end_game := 0.0;
    if nt = 0 then
      Mixed_Solve(q,lifsup,hq,coeffv,expvec,jacmat,mulfac,
                  mix.all,mcc,qsols);
    else
      Silent_Multitasking_Path_Tracker
        (q,nt,p'last,mix'last,mix.all,lifsup,mcc,
         hq,coeffv,expvec,jacmat,mulfac,qsols);
    end if;
    Set_Continuation_Parameter(qsols,Create(zero));
    Continuation_Parameters.start_end_game := 0.1;
    Clear(hq); Clear(jacmat); Clear(mulfac);
    QuadDobl_Complex_VecVecs.Clear(coeffv);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 stlb : in double_float;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 orgmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0
                   : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Laur_Functions;
    use Standard_Complex_Laur_Systems;
    use Standard_Poly_Laur_Convertors;
    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;

    lq : Laur_Sys(p'range);
    hq : Eval_Coeff_Laur_Sys(p'range);
    expvec : Exponent_Vectors_Array(p'range);
    coeffv : Standard_Complex_VecVecs.VecVec(p'range);
    jacmat : Eval_Coeff_Jaco_Mat(p'range,p'first..p'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));

  begin
    if verbose > 0 then
      put("-> in black_polyhedral_continuations.");
      put_line("Black_Box_Polyhedral_Continuation 6");
      put_line("for a polynomial system in double precision ...");
    end if;
    q := Random_Coefficient_Systems.Create(natural32(p'last),mix.all,lifsup);
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
    Continuation_Parameters.Tune(0);
    Continuation_Parameters.start_end_game := 0.0;
    if nt = 0 then
      Mixed_Solve(lq,lifsup,hq,coeffv,expvec,jacmat,mulfac,
                  mix.all,orgmcc,qsols);
    else
      Silent_Multitasking_Path_Tracker
        (lq,nt,p'last,mix'last,mix.all,lifsup,orgmcc,
         hq,coeffv,expvec,jacmat,mulfac,qsols);
    end if;
    Set_Continuation_Parameter(qsols,Create(0.0));
    if not Floating_Mixed_Subdivisions.Is_Null(stbmcc) then
      Silent_Polyhedral_Continuation
        (lq,stlb,mix,lifsup,stbmcc,qsols0,verbose-1);
     -- put_line("looking at the stable mixed cells ...");
     -- Reporting_Polyhedral_Continuation
     --   (standard_output,lq,stlb,mix,lifsup,stbmcc,qsols0);
      Set_Continuation_Parameter(qsols0,Create(0.0));
     -- put("Length_Of(qsols0) = "); put(Length_Of(qsols0),1); new_line;
   -- else
   --   put_line("no stable mixed cells");
    end if;
    Set_Continuation_Parameter(qsols,Create(0.0));
    Continuation_Parameters.start_end_game := 0.1;
    Clear(lq); Clear(hq); Clear(jacmat); Clear(mulfac);
    Standard_Complex_VecVecs.Clear(coeffv);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 stlb : in double_float;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 orgmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0
                   : in out DoblDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Laur_Functions;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Poly_Laur_Convertors;
    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;

    lq : Laur_Sys(p'range);
    hq : Eval_Coeff_Laur_Sys(p'range);
    expvec : Exponent_Vectors_Array(p'range);
    coeffv : DoblDobl_Complex_VecVecs.VecVec(p'range);
    jacmat : Eval_Coeff_Jaco_Mat(p'range,p'first..p'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));
    zero : constant double_double := create(0.0);

  begin
    if verbose > 0 then
      put("-> in black_polyhedral_continuations.");
      put_line("Black_Box_Polyhedral_Continuation 7");
      put_line("for a polynomial system in double double precision ...");
    end if;
    q := Random_Coefficient_Systems.Create(natural32(p'last),mix.all,lifsup);
    lq := Polynomial_to_Laurent_System(q);
    hq := Create(lq);
    expvec := Create(q);
    for i in q'range loop
      declare
        c : constant DoblDobl_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        coeffv(i) := new DoblDobl_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    Create(lq,jacmat,mulfac);
    Continuation_Parameters.Tune(0); --,32);
    Continuation_Parameters.start_end_game := 0.0;
    if nt = 0 then
      Mixed_Solve(lq,lifsup,hq,coeffv,expvec,jacmat,mulfac,
                  mix.all,orgmcc,qsols);
    else
      Silent_Multitasking_Path_Tracker
        (lq,nt,p'last,mix'last,mix.all,lifsup,orgmcc,
         hq,coeffv,expvec,jacmat,mulfac,qsols);
    end if;
    Set_Continuation_Parameter(qsols,Create(zero));
    if not Floating_Mixed_Subdivisions.Is_Null(stbmcc) then
      Silent_Polyhedral_Continuation
        (lq,stlb,mix,lifsup,stbmcc,qsols0,verbose-1);
     -- put_line("looking at the stable mixed cells ...");
     -- Reporting_Polyhedral_Continuation
     --   (standard_output,lq,stlb,mix,lifsup,stbmcc,qsols0);
      Set_Continuation_Parameter(qsols0,Create(zero));
     -- put("Length_Of(qsols0) = "); put(Length_Of(qsols0),1); new_line;
   -- else
   --   put_line("no stable mixed cells");
    end if;
    Set_Continuation_Parameter(qsols,Create(zero));
    Continuation_Parameters.start_end_game := 0.1;
    Clear(lq); Clear(hq); Clear(jacmat); Clear(mulfac);
    DoblDobl_Complex_VecVecs.Clear(coeffv);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 stlb : in double_float;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 orgmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0
                   : in out QuadDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Laur_Functions;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Poly_Laur_Convertors;
    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;

    lq : Laur_Sys(p'range);
    hq : Eval_Coeff_Laur_Sys(p'range);
    expvec : Exponent_Vectors_Array(p'range);
    coeffv : QuadDobl_Complex_VecVecs.VecVec(p'range);
    jacmat : Eval_Coeff_Jaco_Mat(p'range,p'first..p'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));
    zero : constant quad_double := create(0.0);

  begin
    if verbose > 0 then
      put("-> in black_polyhedral_continuations.");
      put_line("Black_Box_Polyhedral_Continuation 8");
      put_line("for a polynomial system in quad double precision ...");
    end if;
    q := Random_Coefficient_Systems.Create(natural32(p'last),mix.all,lifsup);
    lq := Polynomial_to_Laurent_System(q);
    hq := Create(lq);
    expvec := Create(q);
    for i in q'range loop
      declare
        c : constant QuadDobl_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        coeffv(i) := new QuadDobl_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    Create(lq,jacmat,mulfac);
    Continuation_Parameters.Tune(0); --,64);
    Continuation_Parameters.start_end_game := 0.0;
    if nt = 0 then
      Mixed_Solve(lq,lifsup,hq,coeffv,expvec,jacmat,mulfac,
                  mix.all,orgmcc,qsols);
    else
      Silent_Multitasking_Path_Tracker
        (lq,nt,p'last,mix'last,mix.all,lifsup,orgmcc,
         hq,coeffv,expvec,jacmat,mulfac,qsols);
    end if;
    Set_Continuation_Parameter(qsols,Create(zero));
    if not Floating_Mixed_Subdivisions.Is_Null(stbmcc) then
      Silent_Polyhedral_Continuation
        (lq,stlb,mix,lifsup,stbmcc,qsols0,verbose-1);
     -- put_line("looking at the stable mixed cells ...");
     -- Reporting_Polyhedral_Continuation
     --   (standard_output,lq,stlb,mix,lifsup,stbmcc,qsols0);
      Set_Continuation_Parameter(qsols0,Create(zero));
     -- put("Length_Of(qsols0) = "); put(Length_Of(qsols0),1); new_line;
   -- else
   --   put_line("no stable mixed cells");
    end if;
    Set_Continuation_Parameter(qsols,Create(zero));
    Continuation_Parameters.start_end_game := 0.1;
    Clear(lq); Clear(hq); Clear(jacmat); Clear(mulfac);
    QuadDobl_Complex_VecVecs.Clear(coeffv);
  end Black_Box_Polyhedral_Continuation;

  procedure Read_System
              ( file : in out file_type; filename : in string;
                lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

  -- DESCRIPTION :
  --   Attempts to read a polynomial system from file.

  begin
    if filename /= "" then
      Open_Input_File(file,filename);
      get(file,lp);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   lp := null; return;
  end Read_System;

  procedure Call_MixedVol
      ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
        mivo,stmv : out natural32;
        stlb : out double_float;
        lifsup : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
        mix,perm,iprm : out Standard_Integer_Vectors.Link_to_Vector;
        orgmcc,stbmcc : out Floating_Mixed_Subdivisions.Mixed_Subdivision ) is

    mixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv,smv,tmv,orgcnt,stbcnt : natural32;

  begin
    Black_Box_Mixed_Volume_Computation
      (p,mix,perm,iprm,stlb,lifsup,mixsub,orgmcc,stbmcc,
       mv,smv,tmv,orgcnt,stbcnt);
    mivo := mv; stmv := smv;
  end Call_MixedVol;

  procedure Polyhedral_Solver
              ( file : in file_type; nt : in natural32;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    q : Poly_Sys(p'range);
    qsols,qsols0 : Solution_List;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
    stlb : double_float;
    orgmcc,stbmcc : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv,stmv : natural32;

  begin
    tstart(timer);
    Call_MixedVol(p,mv,stmv,stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,mv,1); new_line(file);
    put(file,"stable mixed volume : ");
    put(file,stmv,1); new_line(file);
    new_line(file);
    print_times(file,timer,"Mixed-Volume Computation");
    if mv > 0 then
      tstart(timer);
      Black_Box_Polyhedral_Continuation
        (integer32(nt),p,mix,stlb,lifsup.all,orgmcc,stbmcc,q,qsols,qsols0);
      tstop(timer);
      new_line(file);
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      new_line(file);
      put_line(file,q);
      new_line(file);
      put_line(file,"SOLUTIONS with nonzero coordinates:");
      new_line(file);
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      new_line(file);
      if not Is_Null(qsols0) then
        put_line(file,"SOLUTIONS with zero coordinates :");
        new_line(file);
        put(file,Length_Of(qsols0),natural32(Head_Of(qsols0).n),qsols0);
      end if;
      print_times(file,timer,"Polyhedral Continuation");
    end if;
  end Polyhedral_Solver;

  procedure Main ( nt : in natural32; infilename,outfilename : in string ) is

    use Standard_Complex_Poly_Systems;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;

  begin
    Read_System(infile,infilename,lp);
    if lp = null then
      new_line;
      get(lp);
    end if;
    Close(infile);
    Create_Output_File(outfile,outfilename);
    put(outfile,natural32(lp'last),lp.all);
    Polyhedral_Solver(outfile,nt,lp.all);
    ended_moment := Ada.Calendar.Clock;
    new_line(outfile);
    put(outfile,"PHC ran from ");
    Write_Time_Stamp(outfile,start_moment);
    put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
    put_line(outfile,".");
    Write_Elapsed_Time(outfile,start_moment,ended_moment);
    Close(outfile);
  end Main;

end Black_Polyhedral_Continuations;
