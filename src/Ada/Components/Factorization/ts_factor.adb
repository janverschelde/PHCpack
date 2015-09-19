with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Natural_VecVecs;
with Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Random_Polynomials;        use Standard_Random_Polynomials;
with Sensitivity_of_Factorization;
with Monodromy_Partitions;               use Monodromy_Partitions;
with Interpolate_Multivariate_Factor;    use Interpolate_Multivariate_Factor;
with Multivariate_Factorization;         use Multivariate_Factorization;
with Drivers_to_Factor_Polynomials;      use Drivers_to_Factor_Polynomials;

procedure ts_factor is

-- DESCRIPTION :
--   Test on factorization of complex multivariate polynomials.

-- GENERATING THE INPUT :

  function Generate_Polynomial
              ( n : natural32; d,nt,mu : Standard_Natural_Vectors.Vector )
              return Poly is

  -- DESCRIPTION :
  --   Returns a polynomial with random complex coefficients,
  --   in n variables, as the expanded product of m factors,
  --   each factor of degree d, occurring mu times.

  -- ON ENTRY :
  --   n        number of variables in the polynomial on return;
  --   m        number of factors in the polynomial on return;
  --   d        vector of range 1..m with degrees of the factors;
  --   nt       nt(k) is the number of terms in the k-th factor;
  --   mu       mu(k) is the multiplicity of the k-th factor.

  -- REQUIRED : d'range = nt'range = mu'range.

    res : Poly;
   -- f : Poly := Random(n,d(d'first),nt(nt'first));
    f : Poly := Random_Sparse_Poly(n,d(d'first),nt(nt'first),0);

  begin
    while Degree(f) < integer32(d(d'first)) loop
      Clear(f);
     -- f:= Random(n,d(d'first),nt(nt'first));
      f:= Random_Sparse_Poly(n,d(d'first),nt(nt'first),0);
    end loop;
    Copy(f,res);
    for j in 1..mu(mu'first)-1 loop
      Mul(res,f);
    end loop;
    for i in d'first+1..d'last loop
      declare
       -- f : Poly := Random(n,d(i),nt(i));
        f : Poly := Random_Sparse_Poly(n,d(i),nt(i),0);
      begin
        while Degree(f) < integer32(d(i)) loop
          Clear(f);
         -- f := Random(n,d(i),nt(i));
          f := Random_Sparse_Poly(n,d(i),nt(i),0);
        end loop;
        for j in 1..mu(i) loop
          Mul(res,f);
        end loop;
        Clear(f);
      end;
    end loop;
    return res;
  end Generate_Polynomial;

  procedure Random_Reducible_Template
              ( n,m : out natural32; ldeg,lnbt,
                lmul : out Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Generates a template for a random multivariate reducible polynomial.
  --   The template consists in n, the number of variables, m, the number
  --   of factors, the degrees, number of monomials, and multiplicity of
  --   each factor.

  begin
    n := 0; m := 0;
    new_line;
    put_line("Generating a random multivariate complex polynomial...");
    put("  Give the number of variables : "); get(n);
    put("  Give the number of factors : "); get(m);
    declare
      deg,nbt,mul : Standard_Natural_Vectors.Vector(1..integer32(m));
    begin
      for i in 1..integer32(m) loop
        put("Generating factor "); put(i,1); put_line(" :");
        deg(i) := 0; nbt(i) := 0; mul(i) := 0;
        put("  Give its degree : ");          get(deg(i));
        put("  Give its number of terms : "); get(nbt(i));
        put("  Give its multiplicity : ");    get(mul(i));
      end loop;
      ldeg := new Standard_Natural_Vectors.Vector'(deg);
      lnbt := new Standard_Natural_Vectors.Vector'(nbt);
      lmul := new Standard_Natural_Vectors.Vector'(mul);
    end;
  end Random_Reducible_Template;

  procedure Random_Reducible_Polynomial ( n : out natural32; p : out Poly ) is

  -- DESCRIPTION :
  --   Interactive generation of a random polynomial in n variables.

    m : natural32;
    ldeg,lnbt,lmul : Standard_Natural_Vectors.Link_to_Vector;

  begin
    Random_Reducible_Template(n,m,ldeg,lnbt,lmul);
    p := Generate_Polynomial(n,ldeg.all,lnbt.all,lmul.all);
  end Random_Reducible_Polynomial;

-- REPORTING THE OUTPUT :

  procedure Write_Factors ( file : in file_type; factors : in Poly_Sys;
                            mu : in Standard_Natural_Vectors.Vector ) is
  begin
    for i in factors'range loop
      new_line(file);
      put(file,"factor "); put(file,i,1);
      put_line(file," :"); put_line(file,factors(i));
      put(file,"with multiplicity = "); put(file,mu(i),1);
      new_line(file);
    end loop;
  end Write_Factors;

-- THE MAIN DRIVERS :

  procedure Driver_to_Factor 
              ( file : in file_type; output,monodromy : in boolean;
                n : in natural32; p : in Poly ) is

    timer,total_timer : Timing_Widget;
    d : constant natural32 := natural32(Degree(p));
    b,v : Standard_Complex_Vectors.Vector(1..integer32(n));
    fail : boolean;
    deco : Standard_Natural_VecVecs.Link_to_VecVec;
    wp : Standard_Complex_Vectors.Link_to_Vector;
    mf,mw : Standard_Natural_Vectors.Link_to_Vector;
    rad,dst : Standard_Floating_Vectors.Vector(1..integer32(d));
    rdp,factors : Link_to_Poly_Sys;
    maxdif,maxres : double_float;
    mongrp_time,lintrc_time,itrpol_time,mulval_time,total : duration;

  begin
    tstart(total_timer);
    tstart(timer);
    if monodromy
     then Factor(file,output,p,n,d,deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
     else Trace_Factor(file,p,n,d,deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
    end if;
    tstop(timer);
    if fail then
      put_line(file,"Failed to compute witness points.");
    else
      mongrp_time := Elapsed_User_Time(timer);
      put(file,"r = "); put(file,rad,3); new_line(file);
      put(file,"R = "); put(file,dst,3); new_line(file);
      new_line(file);
      print_times(file,timer,"finding the factorization");
      new_line(file);
      flush(file);
      tstart(timer);
      Certify(file,p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,maxdif);
      tstop(timer);
      lintrc_time := Elapsed_User_Time(timer);
      put(file,"Maximal difference certificate : ");
      put(file,maxdif,3); put_line(file,".");
      new_line(file);
      print_times(file,timer,"linear traces certification");
      new_line(file);
      tstart(timer);
      Interpolate(file,p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,factors);
      tstop(timer);
      itrpol_time := Elapsed_User_Time(timer);
      Write_Factors(file,factors.all,mf.all);
      new_line(file);
      print_times(file,timer,"interpolation of factors");
      new_line(file);
      tstart(timer);
      Multiply_Factors(file,p,factors.all,mf.all,maxres);
      tstop(timer);
      mulval_time := Elapsed_User_Time(timer);
      put(file,"Validation by multiplied factors :");
      put(file,maxres,3); put_line(file,".");
      new_line(file);
      print_times(file,timer,"multiplication of factors");
      new_line(file);
      tstop(total_timer);
      total := Elapsed_User_Time(total_timer);
      Write_Timing_Summary
        (file,mongrp_time,lintrc_time,itrpol_time,mulval_time,total);
    end if;
  end Driver_to_Factor;

  procedure Driver_to_Random_Factor
              ( file : in file_type; output,monodromy : in boolean;
                n : in natural32; p : in Poly; succnt : in out natural32;
                maxdif,maxres : out double_float;
                grp,trc,ipl,val,tot : out duration ) is

  -- DESCRIPTION :
  --   Factors a given polynomial and records timings.

  -- ON ENTRY :
  --   file     file for intermediate output and diagnostics;
  --   output   flag to request output during continuation;
  --   monodromy decides whether monodromy is to be used;
  --   n        number of variables of the polynomial;
  --   p        a random polynomial;
  --   succnt   current number of successful cases.

  -- ON RETURN :
  --   succnt   updated counter if succes;
  --   maxdif   maximal difference in linear trace certificate;
  --   maxres   maximal residual in validation by multiplication;
  --   grp      elapsed user time for monodromy grouping;
  --   trc      elapsed user time for linear trace certificate;
  --   ipl      elapsed user time for interpolation;
  --   val      elapsed user time for validation by multiplication;
  --   tot      total elapsed user time.

    tol : constant double_float := 1.0E-8;
    timer,total_timer : Timing_Widget;
    d : constant natural32 := natural32(Degree(p));
    b,v : Standard_Complex_Vectors.Vector(1..integer32(n));
    fail : boolean;
    deco : Standard_Natural_VecVecs.Link_to_VecVec;
    wp : Standard_Complex_Vectors.Link_to_Vector;
    mf,mw : Standard_Natural_Vectors.Link_to_Vector;
    rad,dst : Standard_Floating_Vectors.Vector(1..integer32(d));
    rdp,factors : Link_to_Poly_Sys;
    mongrp_time,lintrc_time,itrpol_time,mulval_time,total : duration;

  begin
    tstart(total_timer);
    tstart(timer);
    if monodromy
     then Factor(file,output,p,n,d,deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
     else Trace_Factor(file,p,n,d,deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
    end if;
    tstop(timer);
    if fail then
      put_line(file,"Failed to compute witness points.");
      put_line(file,"The polynomial causing the failure : ");
      put_line(file,p); 
      raise PROGRAM_ERROR;
    else
      mongrp_time := Elapsed_User_Time(timer);
      new_line(file);
      put(file,"r = "); put(file,rad,3); new_line(file);
      put(file,"R = "); put(file,dst,3); new_line(file);
      new_line(file);
      print_times(file,timer,"finding the factorization");
      new_line(file);
      flush(file);
      tstart(timer);
      Certify(file,p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,maxdif);
      tstop(timer);
      lintrc_time := Elapsed_User_Time(timer);
      put(file,"Maximal difference certificate : ");
      put(file,maxdif,3); put_line(file,".");
      new_line(file);
      print_times(file,timer,"linear traces certification");
      new_line(file);
      tstart(timer);
      Interpolate(file,p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,factors);
      tstop(timer);
      itrpol_time := Elapsed_User_Time(timer);
      Write_Factors(file,factors.all,mf.all);
      new_line(file);
      print_times(file,timer,"interpolation of factors");
      new_line(file);
      tstart(timer);
      Multiply_Factors(file,p,factors.all,mf.all,maxres);
      tstop(timer);
      mulval_time := Elapsed_User_Time(timer);
      put(file,"Validation by multiplied factors :");
      put(file,maxres,3); put_line(file,".");
      new_line(file);
      print_times(file,timer,"multiplication of factors");
      new_line(file);
      tstop(total_timer);
      total := Elapsed_User_Time(total_timer);
      Write_Timing_Summary
        (file,mongrp_time,lintrc_time,itrpol_time,mulval_time,total);
      put(mf.all); put(maxdif,3); put(maxres,3);
      put(" "); print_hms(Standard_Output,total);
      if maxdif < tol and maxres < tol then
        put_line("  Success");
        succnt := succnt + 1;
      else 
        put_line("  Failure");
      end if;
    end if;
    grp := mongrp_time;
    trc := lintrc_time;
    ipl := itrpol_time;
    val := mulval_time;
    tot := total;
  end Driver_to_Random_Factor;

  procedure Driver_to_Factor ( n : in natural32; p : in Poly ) is

    d : constant natural32 := natural32(Degree(p));
    b,v : Standard_Complex_Vectors.Vector(1..integer32(n));
    fail : boolean;
    deco : Standard_Natural_VecVecs.Link_to_VecVec;
    wp : Standard_Complex_Vectors.Link_to_Vector;
    mf,mw : Standard_Natural_Vectors.Link_to_Vector;
    rad,dst : Standard_Floating_Vectors.Vector(1..integer32(d));
    rdp,factors : Link_to_Poly_Sys;
    maxdif,maxres : double_float;
    timer,total_timer : Timing_Widget;
    mongrp_time,lintrc_time,itrpol_time,mulval_time,total : duration;

  begin
    tstart(total_timer);
    tstart(timer);
    Factor(p,n,d,deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
    tstop(timer);
    if fail then
      put_line("Failed to compute witness points.");
    else
      mongrp_time := Elapsed_User_Time(timer);
      put("r = "); put(rad,3); new_line;
      put("R = "); put(dst,3); new_line;
      put_line("The factorization, with multiplicities : ");
      Write_Factors(Standard_Output,deco.all,mw.all);
      put_line("The witness points : "); put_line(wp.all);
      put("with multiplicities : "); put(mw.all); new_line;
      tstart(timer);
      Certify(p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,maxdif);
      tstop(timer);
      lintrc_time := Elapsed_User_Time(timer);
      put("Maximal difference certificate : ");
      put(maxdif,3); put_line(".");
      tstart(timer);
      Interpolate(p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,factors);
      tstop(timer);
      itrpol_time := Elapsed_User_Time(timer);
      Write_Factors(Standard_Output,factors.all,mf.all);
      tstart(timer);
      Multiply_Factors(p,factors.all,mf.all,maxres);
      tstop(timer);
      mulval_time := Elapsed_User_Time(timer);
      put("Validation by multiplied factors :");
      put(maxres,3); put_line(".");
      tstop(total_timer);
      total := Elapsed_User_Time(total_timer);
      Write_Timing_Summary
        (Standard_Output,mongrp_time,lintrc_time,
         itrpol_time,mulval_time,total);
    end if;
  end Driver_to_Factor;

  procedure Driver_to_Random_Factor
               ( n : in natural32; p : in Poly; succnt : in out natural32;
                 maxdif,maxres : out double_float;
                 grp,trc,ipl,val,tot : out duration ) is

  -- DESCRIPTION :
  --   Writes only one summary line about the factorization of
  --   a random polynomial p.  Counts the number of successes.

  -- ON ENTRY :
  --   n        number of variables of the polynomial;
  --   p        a random polynomial;
  --   succnt   current number of successful cases.

  -- ON RETURN :
  --   succnt   updated counter if succes;
  --   maxdif   maximal difference of linear trace certificates;
  --   maxres   maximal residual of validation by multiplication;
  --   grp      elapsed user time for monodromy grouping;
  --   trc      elapsed user time for linear trace certificate;
  --   ipl      elapsed user time for interpolation;
  --   val      elapsed user time for validation by multiplication;
  --   tot      total elapsed user time.

    tol : constant double_float := 1.0E-8;
    d : constant natural32 := natural32(Degree(p));
    b,v : Standard_Complex_Vectors.Vector(1..integer32(n));
    fail : boolean;
    deco : Standard_Natural_VecVecs.Link_to_VecVec;
    wp : Standard_Complex_Vectors.Link_to_Vector;
    mf,mw : Standard_Natural_Vectors.Link_to_Vector;
    rad,dst : Standard_Floating_Vectors.Vector(1..integer32(d));
    rdp,factors : Link_to_Poly_Sys;
    timer,total_timer : Timing_Widget;
    total : duration;

  begin
    tstart(total_timer);
    tstart(timer);
    Factor(p,n,d,deco,mf,b,v,wp,mw,rdp,rad,dst,fail);
    tstop(timer);
    grp := Elapsed_User_Time(timer);
    if fail then
      put_line("Failed to compute witness points.");
      put_line("The polynomial causing the failure : ");
      put_line(p);
      put_line("The computed witness points : ");
      put_line(wp);
      raise PROGRAM_ERROR;
    else
      put("r = "); put(rad,3); new_line;
      put("R = "); put(dst,3); new_line;
      put(mf.all);
      tstart(timer);
      Certify(p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,maxdif);
      tstop(timer);
      trc := Elapsed_User_Time(timer);
      put(maxdif,3);
      tstart(timer);
      Interpolate(p,b,v,wp.all,mw.all,deco.all,mf.all,rdp,factors);
      tstop(timer);
      ipl := Elapsed_User_Time(timer);
      tstart(timer);
      Multiply_Factors(p,factors.all,mf.all,maxres);
      tstop(timer);
      val := Elapsed_User_Time(timer);
      put(maxres,3);
      tstop(total_timer);
      total := Elapsed_User_Time(total_timer);
      tot := total;
      put(" "); print_hms(Standard_Output,total);
      if maxdif < tol and maxres < tol then
        put_line("  Success");
        succnt := succnt + 1;
      else
        put_line("  Failure");
      end if;
    end if;
  end Driver_to_Random_Factor;

  procedure Factor ( n : in natural32; p : in Poly ) is

    ans : character;
    outfile : file_type;
    output,monodromy : boolean;

  begin
    new_line;
    put("Do you wish intermediate output on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
      new_line;
      put_line("Which factorization method you like ?");
      put_line("  1. using monodromy to group witness points;");
      put_line("  2. combinatorially exploiting linear traces.");
      put("Type 1 or 2 to choose : ");
      Ask_Alternative(ans,"12");
      monodromy := (ans = '1');
      if ans = '1' then
        new_line;
        put("Do you wish intermediate output"
          & " during continuation ? (y/n) ");
        Ask_Yes_or_No(ans);
        output := (ans = 'y');
      end if;
      new_line;
      put_line("See the output file for results...");
      new_line;
      Driver_to_Factor(outfile,output,monodromy,n,p);
    else
      Driver_to_Factor(n,p);
    end if;
  end Factor;

  procedure Write_Summary ( file : in file_type; nb,cnt : in natural32;
                            dif,res : double_float;
                            grp,trc,ipl,val,tot : duration ) is

    dnb : constant duration := duration(nb);
    fnb : constant double_float := double_float(nb);

  begin
    put(file,"Tested "); put(file,nb,1);
    put(file," cases of which ");
    put(file,cnt,1); put_line(file," succeeded.");
    put_line(file,"------------------------------------------------");
    put(file,"Average linear trace certificate     :");
    put(file,dif/fnb,2); new_line(file);
    put(file,"Average validation by multiplication :");
    put(file,res/fnb,2); new_line(file);
    put_line(file,"------------------------------------------------");
    put_line(file,"Average user times : ");
    put(file,"  monodromy grouping          : ");
    print_hms(file,grp/dnb); new_line(file);
    put(file,"  linear traces certification : ");
    print_hms(file,trc/dnb); new_line(file);
    put(file,"  interpolation at factors    : ");
    print_hms(file,ipl/dnb); new_line(file);
    put(file,"  multiplication validation   : ");
    print_hms(file,val/dnb); new_line(file);
    put_line(file,"------------------------------------------------");
    put(file,"  total for all stages        : ");
    print_hms(file,tot/dnb); new_line(file);               
  end Write_Summary;

  procedure Factor_Random_Polynomials
              ( file : in file_type; output,monodromy : in boolean;
                n,m,nb : in natural32;
                deg,nbt,mul : in Standard_Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given the template (n,m,deg,nbt,mul) of a random polynomial,
  --   nb experiments will be performed, writing output to file.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   output   flag to request intermediate output during continuation;
  --   monodromy decides whether to use monodromy or not;
  --   n        number of variables in the template;
  --   m        number of monomials in the template;
  --   nb       number of experiments;
  --   deg      degrees of the factors in the polynomial;
  --   nbt      number of terms in every factor of the polynomial;
  --   mul      multiplicities of the factors.

    cnt : natural32 := 0;
    p : Poly;
    p_dif,p_res : double_float;
    t_dif,t_res : double_float := 0.0;
    p_grp,p_trc,p_ipl,p_val,p_tot : duration;
    t_grp,t_trc,t_ipl,t_val,t_tot : duration := 0.0;

  begin
    for i in 1..nb loop
      p := Generate_Polynomial(n,deg,nbt,mul);
      Normalize(p);
      Driver_to_Random_Factor
        (file,output,monodromy,
         n,p,cnt,p_dif,p_res,p_grp,p_trc,p_ipl,p_val,p_tot);
      t_dif := t_dif + p_dif;
      t_res := t_res + p_res;
      t_grp := t_grp + p_grp;
      t_trc := t_trc + p_trc;
      t_ipl := t_ipl + p_ipl;
      t_val := t_val + p_val;
      t_tot := t_tot + p_tot;
      Clear(p);
    end loop;
    Write_Summary
      (Standard_Output,nb,cnt,t_dif,t_res,t_grp,t_trc,t_ipl,t_val,t_tot);
  end Factor_Random_Polynomials;

  procedure Factor_Random_Polynomials
              ( n,nb : in natural32;
                deg,nbt,mul : in Standard_Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given the template (n,m,deg,nbt,mul) of a random polynomial,
  --   nb experiments will be performed, writing output to screen.

  -- ON ENTRY :
  --   n        number of variables in the template;
  --   m        number of monomials in the template;
  --   nb       number of experiments;
  --   deg      degrees of the factors in the polynomial;
  --   nbt      number of terms in every factor of the polynomial;
  --   mul      multiplicities of the factors.

    cnt : natural32 := 0;
    p : Poly;
    p_dif,p_res : double_float;
    t_dif,t_res : double_float := 0.0;
    p_grp,p_trc,p_ipl,p_val,p_tot : duration;
    t_grp,t_trc,t_ipl,t_val,t_tot : duration := 0.0;

  begin
    for i in 1..nb loop
      p := Generate_Polynomial(n,deg,nbt,mul);
      Normalize(p);
      Driver_to_Random_Factor
        (n,p,cnt,p_dif,p_res,p_grp,p_trc,p_ipl,p_val,p_tot);
      t_dif := t_dif + p_dif;
      t_res := t_res + p_res;
      t_grp := t_grp + p_grp;
      t_trc := t_trc + p_trc;
      t_ipl := t_ipl + p_ipl;
      t_val := t_val + p_val;
      t_tot := t_tot + p_tot;
      Clear(p);
    end loop;
    Write_Summary
      (Standard_Output,nb,cnt,t_dif,t_res,t_grp,t_trc,t_ipl,t_val,t_tot);
  end Factor_Random_Polynomials;

  procedure Factor_Random_Polynomials is

    n,m,nb : natural32 := 0;
    ldeg,lnbt,lmul : Standard_Natural_Vectors.Link_to_Vector;
    p : Poly;
    ans : character;
    outfile : file_type;
    output,monodromy : boolean;

  begin
    Random_Reducible_Template(n,m,ldeg,lnbt,lmul);
    new_line;
    put("How many experiments with this template ? "); get(nb);
    new_line;
    put("Do you wish intermediate output on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
      new_line;
      put_line("Which factorization method you like ?");
      put_line("  1. using monodromy to group witness points;");
      put_line("  2. combinatorially exploiting linear traces.");
      put("Type 1 or 2 to choose : ");
      Ask_Alternative(ans,"12");
      monodromy := (ans = '1');
      if ans = '1' then
        new_line;
        put("Do you wish intermediate output"
          & " during continuation ? (y/n) ");
        Ask_Yes_or_No(ans);
        output := (ans = 'y');
      end if;
      new_line;
      put_line("See the output file for results...");
      new_line;
      Factor_Random_Polynomials
        (outfile,output,monodromy,n,m,nb,ldeg.all,lnbt.all,lmul.all);
    else
      if nb = 1 then
        p := Generate_Polynomial(n,ldeg.all,lnbt.all,lmul.all);
        Normalize(p);
        Driver_to_Factor(n,p);
      else
        Factor_Random_Polynomials(n,nb,ldeg.all,lnbt.all,lmul.all);
      end if;
    end if;
  end Factor_Random_Polynomials;

  procedure Main is

    n : natural32;
    p : Poly;
    ans : character;

  begin
    new_line;
    put_line("Numerical Factorization of Complex Multivariate Polynomials");
    new_line;
    put_line("Choose one of the following options : ");
    put_line("  1. Conduct a random sensitivity experiment.");
    put_line("  2. Factor randomly generated complex polynomials.");
    put_line("  3. Give your own polynomial to factor.");
    put("Type 1, 2, or 3 to choose : "); 
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Sensitivity_of_Factorization; return;
      when '2' => Factor_Random_Polynomials;
      when '3' => Read_Polynomial(n,p);
                  Normalize(p);
                  Factor(n,p);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_factor;
