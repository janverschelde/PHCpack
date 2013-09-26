with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Process_io;                         use Process_io;

package body QuadDobl_Correctors is

-- THE LOOSE CORRECTORS ARE STUBS !

  procedure Affine_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Affine_Single_Loose_Normal_Silent_Corrector;

  procedure Affine_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Affine_Single_Loose_Normal_Reporting_Corrector;

  procedure Affine_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Affine_Single_Loose_Conditioned_Silent_Corrector;

  procedure Affine_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Affine_Single_Loose_Conditioned_Reporting_Corrector;

  procedure Affine_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := Fun(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop        -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := Jef(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when (info /= 0);                  -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(j,y'last,ipvt,y); ncora := hihi_part(Norm(y));
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := Fun(s.sol.v,s.sol.t);  nresa := hihi_part(Norm(y));
      normv := hihi_part(Norm(s.sol.v));
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
      end if;
      nit := nit + 1;
      if nit > 1 then
        stop := (((ncora > s.cora) and then (nresa > s.resa))
          or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                     -- STOP WHEN DIVERGENCE
        exit when stop;
      end if;
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                         -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  --exception
   -- when numeric_error => return;
  --  when others => --put_line("exception in corrector");
  --                 Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Severe_Normal_Silent_Corrector;

  procedure Affine_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := Fun(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    info : integer32;
    normv,ncora,nresa,ncorr,nresr : double_float;
    normj,rcond : quad_double;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := Jef(s.sol.v,s.sol.t);
      wrkj := orgj;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(wrkj,y'last,ipvt,y); ncora := hihi_part(Norm(y));
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := Fun(s.sol.v,s.sol.t);  nresa := hihi_part(Norm(y));
      normv := hihi_part(Norm(s.sol.v));
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
      end if;
      nit := nit + 1;
      if nit > 1 then
        stop := (((ncora > s.cora) and then (nresa > s.resa))
          or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                      -- STOP WHEN DIVERGENCE
        exit when stop;
      end if;
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
    if info = 0 then
      normj := Norm1(orgj);
      estco(wrkj,y'last,ipvt,normj,rcond);       -- ESTIMATE CONDITION NUMBER
      s.rcond := hihi_part(rcond);
    else
      s.rcond := 0.0;
    end if;
  -- exception
   -- when numeric_error => return;
   -- when others => --put_line("exception in corrector");
   --                Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Severe_Conditioned_Silent_Corrector;

  procedure Affine_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := Fun(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,ncorr,nresa,nresr : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := Jef(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;                     -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(j,y'last,ipvt,y); ncora := hihi_part(Norm(y));
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := Fun(s.sol.v,s.sol.t);  nresa := hihi_part(Norm(y));
      normv := hihi_part(Norm(s.sol.v));
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
      end if;
      cWrite(file,ncora,ncorr,nresa,nresr);               -- WRITE STATISTICS
      nit := nit + 1;
      if nit > 1 then
        stop := (((ncora > s.cora) and then (nresa > s.resa))
          or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                      -- STOP WHEN DIVERGENCE
        exit when stop;
      end if;
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
 -- exception
   -- when numeric_error => return;
 --   when others => --put_line("exception in corrector");
 --                  Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Severe_Normal_Reporting_Corrector;

  procedure Affine_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := Fun(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    info : integer32;
    normv,ncora,nresa,ncorr,nresr : double_float;
    normj,rcond : quad_double;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := Jef(s.sol.v,s.sol.t);
      wrkj := orgj;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(wrkj,y'last,ipvt,y); ncora := hihi_part(Norm(y));
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := Fun(s.sol.v,s.sol.t);  nresa := hihi_part(Norm(y));
      normv := hihi_part(Norm(s.sol.v));
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
      end if;
      cWrite(file,ncora,ncorr,nresa,nresr);               -- WRITE STATISTICS
      nit := nit + 1;
      if nit > 1 then
        stop := (((ncora > s.cora) and then (nresa > s.resa))
          or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                      -- STOP WHEN DIVERGENCE
        exit when stop;
      end if;
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
    if info = 0 then
      normj := Norm1(orgj);
      estco(wrkj,y'last,ipvt,normj,rcond);       -- ESTIMATE CONDITION NUMBER
      s.rcond := hihi_part(rcond);
    else
      s.rcond := 0.0;
    end if;
 -- exception
   -- when numeric_error => return;
 --   when others => --put_line("exception in corrector");
 --                  Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Severe_Conditioned_Reporting_Corrector;

-- THE ROUTINES BELOW ARE STUBS :

  procedure Affine_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Affine_Multiple_Loose_Normal_Silent_Corrector;

  procedure Affine_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Affine_Multiple_Loose_Normal_Reporting_Corrector;

  procedure Affine_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Affine_Multiple_Loose_Conditioned_Silent_Corrector;

  procedure Affine_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Affine_Multiple_Loose_Conditioned_Reporting_Corrector;

  procedure Affine_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Affine_Multiple_Severe_Normal_Silent_Corrector;

  procedure Affine_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Affine_Multiple_Severe_Normal_Reporting_Corrector;

  procedure Affine_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Affine_Multiple_Severe_Conditioned_Silent_Corrector;

  procedure Affine_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Affine_Multiple_Severe_Conditioned_Reporting_Corrector;

-- MORE STUBS :

  procedure Projective_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Projective_Single_Loose_Normal_Silent_Corrector;

  procedure Projective_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Projective_Single_Loose_Normal_Reporting_Corrector;

  procedure Projective_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Projective_Single_Loose_Conditioned_Silent_Corrector;

  procedure Projective_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Projective_Single_Loose_Conditioned_Reporting_Corrector;

  procedure Projective_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Projective_Single_Severe_Normal_Silent_Corrector;

  procedure Projective_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Projective_Single_Severe_Normal_Reporting_Corrector;

  procedure Projective_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Projective_Single_Severe_Conditioned_Silent_Corrector;

  procedure Projective_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Projective_Single_Severe_Conditioned_Reporting_Corrector;

  procedure Projective_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Projective_Multiple_Loose_Normal_Silent_Corrector;

  procedure Projective_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Projective_Multiple_Loose_Normal_Reporting_Corrector;

  procedure Projective_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Projective_Multiple_Loose_Conditioned_Silent_Corrector;

  procedure Projective_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Projective_Multiple_Loose_Conditioned_Reporting_Corrector;

  procedure Projective_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Projective_Multiple_Severe_Normal_Silent_Corrector;

  procedure Projective_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Projective_Multiple_Severe_Normal_Reporting_Corrector;

  procedure Projective_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Projective_Multiple_Severe_Conditioned_Silent_Corrector;

  procedure Projective_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is
  begin
    fail := false;
  end Projective_Multiple_Severe_Conditioned_Reporting_Corrector;

end QuadDobl_Correctors;
