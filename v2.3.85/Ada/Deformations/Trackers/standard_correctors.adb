with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Handle_Underflow_Gracefully;        use Handle_Underflow_Gracefully;
with Process_io;                         use Process_io;

package body Standard_Correctors is

-- AUXILIARIES FOR IMPLEMENTING MULTIPLE CORRECTORS :

  procedure Equals ( s : in Solu_Info_Array; x : in Vector; i : in integer32;
                     d : in double_float; j : in out integer32 ) is

    eq : boolean;

  begin
    while j < i loop
      eq := true;
      for k in x'range loop
        eq := Equal(s(j).sol.v(k),x(k),d);
        exit when not eq;
      end loop;
      exit when eq;
      j := j + 1;
    end loop;
  end Equals;

  generic
    with procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars );
  procedure Multiple_Silent_Corrector 
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  -- DESCRIPTION :
  --   Allows to build any multiple silent corrector,
  --   depending on the corrector supplied as generic parameter.

  generic
    with procedure Corrector ( file : in file_type;
                               s : in out Solu_Info; c : in Corr_Pars );
  procedure Multiple_Reporting_Corrector 
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  -- DESCRIPTION :
  --   Allows to build any multiple reporting corrector,
  --   depending on the corrector supplied as generic parameter.

  procedure Multiple_Silent_Corrector 
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    i,j : integer32;

  begin
    fail := false;
    i := pivot;
    loop
      Corrector(s(i),c);
      fail := (((s(i).resa > c.epsaf) and then (s(i).cora > c.epsax))
        or else ((s(i).resr > c.epsrf) and then (s(i).corr > c.epsrx)));
      if not fail then
        j := 1;
        Equals(s,s(i).sol.v,i,dist_sols,j);
        if j /= i
         then fail := true;
        end if;
      end if;
      exit when fail;
      i := i + 1;
      if i > s'last
       then i := s'first;
      end if;
      exit when (i = pivot);
    end loop;
    if fail
     then pivot := i;
    end if;
  end Multiple_Silent_Corrector;

  procedure Multiple_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    i,j : integer32;

  begin
    fail := false;
    i := pivot;
    loop
      Write_path(file,natural32(i));
      Corrector(file,s(i),c);
      sWrite(file,s(i).sol.all);
      fail := (((s(i).resa > c.epsaf) and then (s(i).cora > c.epsax))
        or else ((s(i).resr > c.epsrf) and then (s(i).corr > c.epsrx)));
      if not fail then
        j := 1;
        Equals(s,s(i).sol.v,i,dist_sols,j);
        if j /= i
         then fail := true;
        end if;
      end if;
      exit when fail;
      i := i + 1;
      if i > s'last
       then i := s'first;
      end if;
      exit when (i = pivot);
    end loop;
    if fail
     then pivot := i;
    end if;
  end Multiple_Reporting_Corrector;

-- TARGET PROCEDURES :

  procedure Affine_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : double_float;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;                     -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                         -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return; 
    -- suppress all underflow errors ...
  end Affine_Single_Loose_Normal_Silent_Corrector;

  procedure Affine_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : double_float;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;                     -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      cWrite(file,s.cora,s.corr,s.resa,s.resr);           -- WRITE STATISTICS 
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Loose_Normal_Reporting_Corrector;

  procedure Affine_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nit : natural32 := 0;
    normv,normj : double_float;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := dH(s.sol.v,s.sol.t);
      wrkj := orgj;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(wrkj,y'last,ipvt,y); s.cora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
    if info = 0 then
      normj := Norm1(orgj);
      estco(wrkj,y'last,ipvt,normj,s.rcond); -- ESTIMATE CONDITION NUMBER
    else 
      s.rcond := 0.0;
    end if;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Loose_Conditioned_Silent_Corrector;

  procedure Affine_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    info : integer32;
    normv,normj : double_float;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := dH(s.sol.v,s.sol.t);
      wrkj := orgj;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(wrkj,y'last,ipvt,y); s.cora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      cWrite(file,s.cora,s.corr,s.resa,s.resr);           -- WRITE STATISTICS
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
    if info = 0 then
      normj := Norm1(orgj);
      estco(wrkj,y'last,ipvt,normj,s.rcond); -- ESTIMATE CONDITION NUMBER
    else
      s.rcond := 0.0;
    end if;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Loose_Conditioned_Reporting_Corrector;

  procedure Affine_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop        -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when (info /= 0);                  -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  nresa := Norm(y);
      normv := Norm(s.sol.v);
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
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Severe_Normal_Silent_Corrector;

  procedure Affine_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;                     -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  nresa := Norm(y);
      normv := Norm(s.sol.v);
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
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Severe_Normal_Reporting_Corrector;

  procedure Affine_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    info : integer32;
    normv,ncora,nresa,ncorr,nresr,normj : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := dH(s.sol.v,s.sol.t);
      wrkj := orgj;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(wrkj,y'last,ipvt,y); ncora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  nresa := Norm(y);
      normv := Norm(s.sol.v);
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
      estco(wrkj,y'last,ipvt,normj,s.rcond); -- ESTIMATE CONDITION NUMBER
    else
      s.rcond := 0.0;
    end if;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Severe_Conditioned_Silent_Corrector;

  procedure Affine_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    info : integer32;
    normv,ncora,nresa,ncorr,nresr,normj : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := dH(s.sol.v,s.sol.t);
      wrkj := orgj;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      lusolve(wrkj,y'last,ipvt,y); ncora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  nresa := Norm(y);
      normv := Norm(s.sol.v);
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
      estco(wrkj,y'last,ipvt,normj,s.rcond); -- ESTIMATE CONDITION NUMBER
    else
      s.rcond := 0.0;
    end if;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Severe_Conditioned_Reporting_Corrector;

  procedure Affine_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Normal_Silent_Corrector;

  procedure Affine_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Normal_Reporting_Corrector;

  procedure Affine_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Conditioned_Silent_Corrector;

  procedure Affine_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Conditioned_Reporting_Corrector;

  procedure Affine_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Normal_Silent_Corrector;

  procedure Affine_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Normal_Reporting_Corrector;

  procedure Affine_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Conditioned_Silent_Corrector;

  procedure Affine_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Conditioned_Reporting_Corrector;

  procedure Projective_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : double_float;

  begin
    while nit < c.maxit loop        -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);             -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;                    -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      y(y'last) := Create(0.0);                   -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := Create(0.0);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0 then
        s.corr := s.cora / normv; s.resr := s.resa / normv;
        for jj in s.sol.v'range loop               -- SCALE THE SOLUTION
          s.sol.v(jj) := s.sol.v(jj)/normv;
        end loop;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
         and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                         -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Projective_Single_Loose_Normal_Silent_Corrector;

  procedure Projective_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : double_float;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);              -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;                     -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      y(y'last) := Create(0.0);                    -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := Create(0.0);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0 then
        s.corr := s.cora / normv; s.resr := s.resa / normv;
        for jj in s.sol.v'range loop                -- SCALE THE SOLUTION
          s.sol.v(jj) := s.sol.v(jj)/normv;
        end loop;
      end if;
      cWrite(file,s.cora,s.corr,s.resa,s.resr);           -- WRITE STATISTICS
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Projective_Single_Loose_Normal_Reporting_Corrector;

  procedure Projective_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    info : integer32;
    normv,normj : double_float;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        orgj(orgj'last(1),jj) := s.sol.v(jj);        -- CORRECT PERPENDICULAR
      end loop;
      wrkj := orgj;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      y(y'last) := Create(0.0);                    -- IGNORE SCALING EQUATION
      lusolve(wrkj,y'last,ipvt,y); s.cora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := Create(0.0);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0 then
        s.corr := s.cora / normv; s.resr := s.resa / normv;
        for jj in s.sol.v'range loop                -- SCALE THE SOLUTION
          s.sol.v(jj) := s.sol.v(jj)/normv;
        end loop;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
    if info = 0 then
      normj := Norm1(orgj);
      estco(wrkj,y'last,ipvt,normj,s.rcond); -- ESTIMATE CONDITION NUMBER
    else
      s.rcond := 0.0;
    end if;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Projective_Single_Loose_Conditioned_Silent_Corrector;

  procedure Projective_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nit : natural32 := 0;
    normv,normj : double_float;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        orgj(orgj'last(1),jj) := s.sol.v(jj);        -- CORRECT PERPENDICULAR
      end loop;
      wrkj := orgj;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      y(y'last) := Create(0.0);                    -- IGNORE SCALING EQUATION
      lusolve(wrkj,y'last,ipvt,y); s.cora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := Create(0.0);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0 then
        s.corr := s.cora / normv; s.resr := s.resa / normv;
        for jj in s.sol.v'range loop                -- SCALE THE SOLUTION
          s.sol.v(jj) := s.sol.v(jj)/normv;
        end loop;
      end if;
      cWrite(file,s.cora,s.corr,s.resa,s.resr);           -- WRITE STATISTICS
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
    if info = 0 then
      normj := Norm1(orgj);
      estco(wrkj,y'last,ipvt,normj,s.rcond); -- ESTIMATE CONDITION NUMBER
    else
      s.rcond := 0.0;
    end if;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Projective_Single_Loose_Conditioned_Reporting_Corrector;

  procedure Projective_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);              -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;                     -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      y(y'last) := Create(0.0);                    -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := Create(0.0);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0 then
        ncorr := ncora / normv; nresr := nresa / normv;
        for jj in s.sol.v'range loop                -- SCALE THE SOLUTION
          s.sol.v(jj) := s.sol.v(jj)/normv;
        end loop;
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
    end loop;                   -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Projective_Single_Severe_Normal_Silent_Corrector;

  procedure Projective_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);              -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;                     -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      y(y'last) := Create(0.0);                    -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := Create(0.0);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0 then
        ncorr := ncora / normv; nresr := nresa / normv;
        for jj in s.sol.v'range loop                -- SCALE THE SOLUTION
          s.sol.v(jj) := s.sol.v(jj)/normv;
        end loop;
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
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Projective_Single_Severe_Normal_Reporting_Corrector;

  procedure Projective_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr,normj : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        orgj(orgj'last(1),jj) := s.sol.v(jj);        -- CORRECT PERPENDICULAR
      end loop;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      y(y'last) := Create(0.0);                    -- IGNORE SCALING EQUATION
      lusolve(wrkj,y'last,ipvt,y); ncora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := Create(0.0);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0 then
        ncorr := ncora / normv; nresr := nresa / normv;
        for jj in s.sol.v'range loop                -- SCALE THE SOLUTION
          s.sol.v(jj) := s.sol.v(jj)/normv;
        end loop;
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
      estco(wrkj,y'last,ipvt,normj,s.rcond); -- ESTIMATE CONDITION NUMBER
    else
      s.rcond := 0.0;
    end if;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Projective_Single_Severe_Conditioned_Silent_Corrector;

  procedure Projective_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    wrkj,orgj : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    info : integer32;
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr,normj : double_float;
    stop : boolean := false;

  begin
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      orgj := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        orgj(orgj'last(1),jj) := s.sol.v(jj);        -- CORRECT PERPENDICULAR
      end loop;
      wrkj := orgj;
      lufac(wrkj,y'last,ipvt,info);
      exit when (info /= 0);                   -- STOP WHEN SINGULAR JACOBIAN
      Min(y);
      y(y'last) := Create(0.0);                    -- IGNORE SCALING EQUATION
      lusolve(wrkj,y'last,ipvt,y); ncora := Norm(y);
      Add(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := Create(0.0);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0 then
        ncorr := ncora / normv; nresr := nresa / normv;
        for jj in s.sol.v'range loop                -- SCALE THE SOLUTION
          s.sol.v(jj) := s.sol.v(jj)/normv;
        end loop;
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
      estco(wrkj,y'last,ipvt,normj,s.rcond); -- ESTIMATE CONDITION NUMBER
    else
      s.rcond := 0.0;
    end if;
  exception
   -- when numeric_error => return;
    when others => --put_line("exception in corrector");
                   Underflow_to_Zero(s.sol.v); return;
  end Projective_Single_Severe_Conditioned_Reporting_Corrector;

  procedure Projective_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Normal_Silent_Corrector;

  procedure Projective_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Normal_Reporting_Corrector;

  procedure Projective_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Conditioned_Silent_Corrector;

  procedure Projective_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Conditioned_Reporting_Corrector;

  procedure Projective_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Normal_Silent_Corrector;

  procedure Projective_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Normal_Reporting_Corrector;

  procedure Projective_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Conditioned_Silent_Corrector;

  procedure Projective_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Conditioned_Reporting_Corrector;

end Standard_Correctors;
