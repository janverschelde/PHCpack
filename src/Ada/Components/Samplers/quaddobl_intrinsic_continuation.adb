with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers_Polar;    use QuadDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Singular_Values;  use QuadDobl_Complex_Singular_Values;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with QuadDobl_Intrinsic_Newton;         use QuadDobl_Intrinsic_Newton;
with QuadDobl_Intrinsic_Trackers;       use QuadDobl_Intrinsic_Trackers;

package body QuadDobl_Intrinsic_Continuation is

  tol : constant double_float := 1.0E-8;

-- AUXILIARIES FOR VALIDATION :

  procedure Write_Banner ( file : in file_type ) is
  begin
    for i in 1..78 loop
      put(file,"=");
    end loop;
    new_line(file);
  end Write_Banner;

  function Is_Clustered
             ( s : Solu_Info_Array; i,j : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if solutions i and j are clustered w.r.t. tol.

  begin
    for k in s(i).sol.v'range loop
      if AbsVal(s(i).sol.v(k) - s(j).sol.v(k)) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Is_Clustered;

  function Is_Real ( x : Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if all imaginary parts are smaller than tol.

  begin
    for i in x'range loop
      if AbsVal(IMAG_PART(x(i))) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Is_Real;

  procedure Root_Accounting
              ( file : in file_type; s : in out Solu_Info_Array;
                k : in integer32; fail : in boolean;
                nbregu,nbreal,nbcomp,nbsing,nbclus,nbfail
                  : in out natural32 ) is

  -- DESCRIPTION :
  --   The main purpose of this routine is to detect clustering,
  --   but it also scans for other characteristics of roots.

    isclus : boolean := false;

  begin
    if fail then
      put_line(file," no solution ==");
      nbfail := nbfail + 1;
    else
      for i in s'first..k-1 loop
        if Is_Clustered(s,i,k) then
          put(file," clustered with ");
          put(file,i,1); put_line(" ==");
          s(k).sol.m := i;
          isclus := true;
        end if;
        exit when isclus;
      end loop;
      if isclus then
        nbclus := nbclus + 1;
        if s(k).sol.rco > tol
         then nbregu := nbregu + 1;
         else nbsing := nbsing + 1;
        end if;
        if Is_Real(s(k).sol.v)
         then nbreal := nbreal + 1;
         else nbcomp := nbcomp + 1;
        end if;
      elsif s(k).sol.rco > tol then
        nbregu := nbregu + 1;
        if Is_Real(s(k).sol.v)
         then put_line(file," real regular ==");
              nbreal := nbreal + 1;
         else put_line(file," complex regular ==");
              nbcomp := nbcomp + 1;
        end if;
      else
        nbsing := nbsing + 1;
        if Is_Real(s(k).sol.v)
         then put_line(file," real singular ==");
              nbreal := nbreal + 1;
         else put_line(file," complex singular ==");
              nbcomp := nbcomp + 1;
        end if;
        s(k).sol.m := 0;
      end if;
    end if;
  end Root_Accounting;

  procedure Root_Accounting
              ( s : in out Solu_Info_Array; k : in integer32;
                fail : in boolean;
                nbregu,nbreal,nbcomp,nbsing,nbclus,nbfail
                  : in out natural32 ) is

  -- DESCRIPTION :
  --   The main purpose of this routine is to detect clustering,
  --   but it also scans for other characteristics of roots.

    isclus : boolean := false;

  begin
    if fail then
      nbfail := nbfail + 1;
    else
      for i in s'first..k-1 loop
        if Is_Clustered(s,i,k) then
          s(k).sol.m := i;
          isclus := true;
        end if;
        exit when isclus;
      end loop;
      if isclus then
        nbclus := nbclus + 1;
        if s(k).sol.rco > tol
         then nbregu := nbregu + 1;
         else nbsing := nbsing + 1;
        end if;
        if Is_Real(s(k).sol.v)
         then nbreal := nbreal + 1;
         else nbcomp := nbcomp + 1;
        end if;
      elsif s(k).sol.rco > tol then
        nbregu := nbregu + 1;
        if Is_Real(s(k).sol.v)
         then nbreal := nbreal + 1;
         else nbcomp := nbcomp + 1;
        end if;
      else
        nbsing := nbsing + 1;
        if Is_Real(s(k).sol.v)
         then nbreal := nbreal + 1;
         else nbcomp := nbcomp + 1;
        end if;
        s(k).sol.m := 0;
      end if;
    end if;
  end Root_Accounting;

  procedure Write_Report
              ( file : in file_type; nb,nbregu,nbreal,nbcomp,
                nbsing,nbclus,nbfail : in natural32 ) is
  begin
    put(file,"A list of "); put(file,nb,1);
    put_line(file," solutions has been refined :");
    put(file,"Number of regular solutions   : ");
    put(file,nbregu,1); put_line(file,".");
    put(file,"Number of singular solutions  : ");
    put(file,nbsing,1); put_line(file,".");
    put(file,"Number of real solutions      : ");
    put(file,nbreal,1); put_line(file,".");
    put(file,"Number of complex solutions   : ");
    put(file,nbcomp,1); put_line(file,".");
    put(file,"Number of clustered solutions : ");
    put(file,nbclus,1); put_line(file,".");
    put(file,"Number of failures            : ");
    put(file,nbfail,1); put_line(file,".");
  end Write_Report;

  procedure Report ( file : in file_type;
                     i : in natural32; s : in Solu_Info ) is

  -- DESCRIPTION :
  --   Reports on the continuation of the i-th solution path.

  begin
    put(file,"== "); put(file,i,1); put(file," = ");
    put(file,"#step : "); put(file,s.nstep,1); put(file," = ");
    put(file,"#fail : "); put(file,s.nfail,1); put(file," = ");
    put(file,"#iter : "); put(file,s.niter,1); put(file," = ");
    if REAL_PART(s.sol.t) < 1.0
     then put_line(file," failure ==");
     else put_line(file," success ==");
    end if;
    s.sol.err := create(s.cora);
    s.sol.rco := create(s.rcond);
    s.sol.res := create(s.resa);
    put(file,s.sol.all); new_line(file);  -- 01/13 uncommented
    put_diagnostics(file,s.sol.all); new_line(file);
  end Report;

-- CONTINUATION ROUTINES :

  procedure Silent_Affine_LU_Continue
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

    procedure Track is new Silent_Affine_LU_Track(Path);

  begin
    for i in s'range loop
      Track(f,jf,s(i),p,c);
    end loop;
  end Silent_Affine_LU_Continue;

  procedure Silent_Projective_LU_Continue
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

    procedure Track is new Silent_Projective_LU_Track(Path);
    k : constant natural32 := 0;

  begin
    for i in s'range loop
      Track(f,jf,s(i),k,p,c);
    end loop;
  end Silent_Projective_LU_Continue;

  procedure G_Silent_LU_Continue
               ( ne,nv : in natural32; s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

    procedure Track is new G_Silent_LU_Track(f,jf,Path);

  begin
    for i in s'range loop
      Track(ne,nv,s(i),p,c);
    end loop;
  end G_Silent_LU_Continue;

  procedure Reporting_Affine_LU_Continue
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

    procedure Track is new Reporting_Affine_LU_Track(Path);
   -- procedure Track is new Silent_Affine_LU_Track(Path);

  begin
    for i in s'range loop
      Track(file,f,jf,s(i),p,c);
     -- Track(f,jf,s(i),p,c);
      Report(file,natural32(i),s(i));
    end loop;
  end Reporting_Affine_LU_Continue;

  procedure Reporting_Projective_LU_Continue
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

    procedure Track is new Reporting_Projective_LU_Track(Path);
   -- procedure Track is new Silent_Projective_LU_Track(Path);
    k : constant natural32 := 0;

  begin
    for i in s'range loop
      Track(file,f,jf,s(i),k,p,c);
     -- Track(f,jf,s(i),k,p,c);
      Report(file,natural32(i),s(i));
    end loop;
  end Reporting_Projective_LU_Continue;

  procedure G_Reporting_LU_Continue
               ( file : in file_type; 
                 ne,nv : in natural32; s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

    procedure Track is new G_Reporting_LU_Track(f,jf,Path);
   -- procedure Track is new G_Silent_LU_Track(f,jf,Path);

  begin
    for i in s'range loop
      Track(file,ne,nv,s(i),p,c);
     -- Track(ne,nv,s(i),p,c);
      Report(file,natural32(i),s(i));
    end loop;
  end G_Reporting_LU_Continue;

  procedure Silent_QR_Continue
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

    procedure Track is new Silent_QR_Track(Path);

  begin
    for i in s'range loop
      Track(f,jf,s(i),p,c);
    end loop;
  end Silent_QR_Continue;

  procedure G_Silent_QR_Continue
               ( ne,nv : in natural32; s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

    procedure Track is new G_Silent_QR_Track(f,jf,Path);

  begin
    for i in s'range loop
      Track(ne,nv,s(i),p,c);
    end loop;
  end G_Silent_QR_Continue;

  procedure Reporting_QR_Continue
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

   -- procedure Track is new Reporting_QR_Track(Path);
    procedure Track is new Silent_QR_Track(Path);

  begin
    for i in s'range loop
     -- Track(file,f,jf,s(i),p,c);
      Track(f,jf,s(i),p,c);
      Report(file,natural32(i),s(i));
    end loop;
  end Reporting_QR_Continue;

  procedure G_Reporting_QR_Continue
               ( file : in file_type; 
                 ne,nv : in natural32; s : in out Solu_Info_Array;
                 p : in Pred_Pars; c : in Corr_Pars ) is

   -- procedure Track is new G_Reporting_QR_Track(f,jf,Path);
    procedure Track is new G_Silent_QR_Track(f,jf,Path);

  begin
    for i in s'range loop
     -- Track(file,ne,nv,s(i),p,c);
      Track(ne,nv,s(i),p,c);
      Report(file,natural32(i),s(i));
    end loop;
  end G_Reporting_QR_Continue;

-- VALIDATION ROUTINES :

  procedure LU_Validate
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 ) is

    nbit : natural32;
    fail : boolean;
    nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail : natural32 := 0;
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    for i in s'range loop
      Affine_LU_Newton
        (f,jf,plane,s(i).sol.v,
         create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
         scora,scorr,sresa,sresr,nbit,c.maxit,srcond,fail);
      s(i).cora := to_double(scora);
      s(i).corr := to_double(scorr);
      s(i).resa := to_double(sresa);
      s(i).resr := to_double(sresr);
      s(i).rcond := to_double(srcond);
      Root_Accounting(s,i,fail,nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    end loop;
    rg := nbregu; sn := nbsing;
    rl := nbreal; cm := nbcomp;
    cl := nbclus; fl := nbfail;
  end LU_Validate;

  procedure Silent_LU_Validate
               ( n : in natural32;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 ) is

    nbit : natural32;
    fail : boolean;
    nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail : natural32 := 0;
    scora,scorr,sresa,sresr,srcond : quad_double;

    procedure Newton is new Silent_Affine_LU_RCO_Newton(f,jf);

  begin
    for i in s'range loop
      Newton
        (n,plane,s(i).sol.v,
         create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
         scora,scorr,sresa,sresr,nbit,c.maxit,srcond,fail);
      s(i).cora := to_double(scora);
      s(i).corr := to_double(scorr);
      s(i).resa := to_double(sresa);
      s(i).resr := to_double(sresr);
      s(i).rcond := to_double(srcond);
      Root_Accounting(s,i,fail,nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    end loop;
    rg := nbregu; sn := nbsing;
    rl := nbreal; cm := nbcomp;
    cl := nbclus; fl := nbfail;
  end Silent_LU_Validate;

  procedure LU_Validate
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 ) is

    nbit : natural32;
    fail : boolean;
    nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail : natural32 := 0;
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,s'last,1); put(file," ");
    put(file,plane'last,1); new_line(file);
    Write_Banner(file);
    for i in s'range loop
      put(file,"solution "); put(file,i,1); put(file," : ");
      put(file,"  start residual :"); put(file,s(i).resa,3);
      Affine_LU_Newton
        (f,jf,plane,s(i).sol.v,
         create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
         scora,scorr,sresa,sresr,nbit,c.maxit,srcond,fail);
      s(i).cora := to_double(scora);
      s(i).corr := to_double(scorr);
      s(i).resa := to_double(sresa);
      s(i).resr := to_double(sresr);
      s(i).rcond := to_double(srcond);
      put(file,"  #iterations : "); put(file,s(i).niter+nbit,1);
      if fail
       then put_line(file,"  failure");
       else put_line(file,"  success");
      end if;
      s(i).sol.err := scora;
      s(i).sol.rco := srcond;
      s(i).sol.res := sresa;
      put(file,s(i).sol.all);
      Root_Accounting(file,s,i,fail,nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    end loop;
    Write_Banner(file);
    Write_Report(file,natural32(s'last),
                 nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    Write_Banner(file);
    rg := nbregu; sn := nbsing;
    rl := nbreal; cm := nbcomp;
    cl := nbclus; fl := nbfail;
  end LU_Validate;

  procedure Reporting_LU_Validate
               ( file : in file_type; n : in natural32;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 ) is

    nbit : natural32;
    fail : boolean;
    nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail : natural32 := 0;
    scora,scorr,sresa,sresr,srcond : quad_double;

    procedure Newton is new Silent_Affine_LU_RCO_Newton(f,jf);

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,s'last,1); put(file," ");
    put(file,plane'last,1); new_line(file);
    Write_Banner(file);
    for i in s'range loop
      put(file,"solution "); put(file,i,1); put(file," : ");
      put(file,"  start residual :"); put(file,s(i).resa,3);
      Newton
        (n,plane,s(i).sol.v,
         create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
         scora,scorr,sresa,sresr,nbit,c.maxit,srcond,fail);
      s(i).cora := to_double(scora);
      s(i).corr := to_double(scorr);
      s(i).resa := to_double(sresa);
      s(i).resr := to_double(sresr);
      s(i).rcond := to_double(srcond);
      put(file,"  #iterations : "); put(file,s(i).niter+nbit,1);
      if fail
       then put_line(file,"  failure");
       else put_line(file,"  success");
      end if;
      s(i).sol.err := scora;
      s(i).sol.rco := srcond;
      s(i).sol.res := sresa;
      put(file,s(i).sol.all);
      Root_Accounting(file,s,i,fail,nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    end loop;
    Write_Banner(file);
    Write_Report(file,natural32(s'last),
                 nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    Write_Banner(file);
    rg := nbregu; sn := nbsing;
    rl := nbreal; cm := nbcomp;
    cl := nbclus; fl := nbfail;
  end Reporting_LU_Validate;

  procedure SV_Validate
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 ) is

    mm : constant integer32 := Min0(f'last+1,plane'last);
    sv : Vector(1..mm);
    nbit : natural32;
    fail : boolean;
    nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail : natural32 := 0;
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    for i in s'range loop
      Affine_SV_Newton
        (f,jf,plane,s(i).sol.v,
         create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
         scora,scorr,sresa,sresr,nbit,c.maxit,sv,fail);
      s(i).cora := to_double(scora);
      s(i).corr := to_double(scorr);
      s(i).resa := to_double(sresa);
      s(i).resr := to_double(sresr);
      srcond := Radius(sv(s(i).sol.v'last)/sv(sv'first));
      s(i).rcond := to_double(srcond);
      Root_Accounting(s,i,fail,nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    end loop;
    rg := nbregu; sn := nbsing;
    rl := nbreal; cm := nbcomp;
    cl := nbclus; fl := nbfail;
  end SV_Validate;

  procedure Silent_SV_Validate
               ( n : in natural32;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 ) is

    mm : constant integer32 := Min0(integer32(n)+1,plane'last);
    sv : Vector(1..mm);
    nbit : natural32;
    fail : boolean;
    nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail : natural32 := 0;
    scora,scorr,sresa,sresr,srcond : quad_double;

    procedure Newton is new Silent_Affine_SV_Newton(f,jf);

  begin
    for i in s'range loop
      Newton(n,plane,s(i).sol.v,
             create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
             scora,scorr,sresa,sresr,nbit,c.maxit,sv,fail);
      s(i).cora := to_double(scora);
      s(i).corr := to_double(scorr);
      s(i).resa := to_double(sresa);
      s(i).resr := to_double(sresr);
      srcond := Radius(sv(s(i).sol.v'last)/sv(sv'first));
      s(i).rcond := to_double(srcond);
      Root_Accounting(s,i,fail,nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    end loop;
    rg := nbregu; sn := nbsing;
    rl := nbreal; cm := nbcomp;
    cl := nbclus; fl := nbfail;
  end Silent_SV_Validate;

  procedure SV_Validate
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 ) is

    mm : constant integer32 := Min0(f'last+1,plane'last);
    sv : Vector(1..mm);
    nbit : natural32;
    fail : boolean;
    nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail : natural32 := 0;
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,s'last,1); put(file," ");
    put(file,plane'last,1); new_line(file);
    Write_Banner(file);
    for i in s'range loop
      put(file,"solution "); put(file,i,1); put(file," : ");
      put(file,"  start residual :"); put(file,s(i).resa,3);
      Affine_SV_Newton
        (f,jf,plane,s(i).sol.v,
         create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
         scora,scorr,sresa,sresr,nbit,c.maxit,sv,fail);
      s(i).cora := to_double(scora);
      s(i).corr := to_double(scorr);
      s(i).resa := to_double(sresa);
      s(i).resr := to_double(sresr);
      srcond := Radius(sv(s(i).sol.v'last)/sv(sv'first));
      s(i).rcond := to_double(srcond);
      put(file,"  #iterations : "); put(file,s(i).niter+nbit,1);
      if fail
       then put_line(file,"  failure");
       else put_line(file,"  success");
      end if;
      s(i).sol.err := scora;
      s(i).sol.rco := srcond;
      s(i).sol.res := sresa;
      put(file,s(i).sol.all);
      Root_Accounting(file,s,i,fail,nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    end loop;
    Write_Banner(file);
    Write_Report(file,natural32(s'last),
                 nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    Write_Banner(file);
    rg := nbregu; sn := nbsing;
    rl := nbreal; cm := nbcomp;
    cl := nbclus; fl := nbfail;
  end SV_Validate;

  procedure Reporting_SV_Validate
               ( file : in file_type; n : in natural32;
                 plane : in Matrix; s : in out Solu_Info_Array;
                 c : in Corr_Pars; rg,sn,rl,cm,cl,fl : out natural32 ) is

    mm : constant integer32 := Min0(integer32(n)+1,plane'last);
    sv : Vector(1..mm);
    nbit : natural32;
    fail : boolean;
    nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail : natural32 := 0;
    scora,scorr,sresa,sresr,srcond : quad_double;

    procedure Newton is new Silent_Affine_SV_Newton(f,jf);

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,s'last,1); put(file," ");
    put(file,plane'last,1); new_line(file);
    Write_Banner(file);
    for i in s'range loop
      put(file,"solution "); put(file,i,1); put(file," : ");
      put(file,"  start residual :"); put(file,s(i).resa,3);
      Newton(n,plane,s(i).sol.v,
             create(c.epsax),create(c.epsrx),create(c.epsaf),create(c.epsrf),
             scora,scorr,sresa,sresr,nbit,c.maxit,sv,fail);
      s(i).cora := to_double(scora);
      s(i).corr := to_double(scorr);
      s(i).resa := to_double(sresa);
      s(i).resr := to_double(sresr);
      srcond := Radius(sv(s(i).sol.v'last)/sv(sv'first));
      s(i).rcond := to_double(srcond);
      put(file,"  #iterations : "); put(file,s(i).niter+nbit,1);
      if fail
       then put_line(file,"  failure");
       else put_line(file,"  success");
      end if;
      s(i).sol.err := scora;
      s(i).sol.rco := srcond;
      s(i).sol.res := sresa;
      put(file,s(i).sol.all);
      Root_Accounting(file,s,i,fail,nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    end loop;
    Write_Banner(file);
    Write_Report(file,natural32(s'last),
                 nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail);
    Write_Banner(file);
    rg := nbregu; sn := nbsing;
    rl := nbreal; cm := nbcomp;
    cl := nbclus; fl := nbfail;
  end Reporting_SV_Validate;

-- using LOCAL COORDINATES :

  procedure Silent_Local_LU_Continue
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 start,target : in Matrix; reoriented : in boolean;
                 s : in out Solu_Info_Array;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    fail : boolean;

  begin
    for i in s'range loop
      Silent_Recentered_LU_Track 
        (f,jf,start,target,reoriented,s(i),pp,cp,fail);
    end loop;
  end Silent_Local_LU_Continue;

  procedure Reporting_Local_LU_Continue
               ( file : in file_type;
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 start,target : in Matrix; reoriented : in boolean;
                 s : in out Solu_Info_Array;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    fail : boolean;
    fail_cnt : natural32 := 0;
    n : constant integer32 := f'last;
   -- mm : constant natural := Min0(n+1,target'last);
   -- sv : Vector(1..mm);
   -- nbit : natural;
    nbregu,nbsing,nbreal,nbcomp,nbclus,nbfail : natural32 := 0;

  begin
    for i in s'range loop
      put(file,"tracking solution path ");
      put(file,i,1); put_line(file," ...");
      Reporting_Recentered_LU_Track
        (file,f,jf,start,target,reoriented,s(i),pp,cp,fail);
       -- Affine_SV_Newton
       --   (file,f,jf,p,sol.sol.v,cp.epsax,cp.epsrx,cp.epsaf,cp.epsrf,
       --    sol.cora,sol.corr,sol.resa,sol.resr,nbit,cp.maxit,sv,fail);
       -- sol.rcond := Radius(sv(sol.sol.v'last)/sv(sv'first));
      put(file,i,3);
      put(file," #step "); put(file,s(i).nstep,1);
      put(file," #fail "); put(file,s(i).nfail,1);
      put(file," #iter "); put(file,s(i).niter,1);
      put(file," err "); put(file,s(i).cora,2);
      put(file," rco "); put(file,s(i).rcond,2);
      put(file," res "); put(file,s(i).resa,2);
      if fail
       then put_line(file," FAIL"); fail_cnt := fail_cnt + 1;
       else put_line(file," ok");
      end if;
    end loop;
  end Reporting_Local_LU_Continue;

end QuadDobl_Intrinsic_Continuation;
