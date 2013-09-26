with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;
with Process_io;                         use Process_io;

package body Multprec_Correctors is

-- AUXILIARIES FOR IMPLEMENTING MULTIPLE CORRECTORS :

  procedure Equals ( s : in Solu_Info_Array; x : in Vector; i : in integer32;
                     d : in Floating_Number; j : in out integer32 ) is

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
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  -- DESCRIPTION :
  --   Allows to build any multiple silent corrector,
  --   depending on the corrector supplied as generic parameter.

  generic
    with procedure Corrector ( file : in file_type;
                               s : in out Solu_Info; c : in Corr_Pars );
  procedure Multiple_Reporting_Corrector 
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean );

  -- DESCRIPTION :
  --   Allows to build any multiple reporting corrector,
  --   depending on the corrector supplied as generic parameter.

  procedure Multiple_Silent_Corrector 
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    i,j : integer32;
    ffail : boolean := false;

  begin
    i := pivot;
    loop
      Corrector(s(i),c);
      ffail := (((s(i).resa > c.epsaf) and then (s(i).cora > c.epsax))
        or else ((s(i).resr > c.epsrf) and then (s(i).corr > c.epsrx)));
      if not ffail then
        j := 1;
        Equals(s,s(i).sol.v,i,dist_sols,j);
        if j /= i
         then ffail := true;
        end if;
      end if;
      exit when ffail;
      i := i + 1;
      if i > s'last
       then i := s'first;
      end if;
      exit when (i = pivot);
    end loop;
    if ffail
     then pivot := i;
    end if;
    fail := ffail;
  end Multiple_Silent_Corrector;

  procedure Multiple_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    i,j : integer32;
    ffail : boolean := false;

  begin
    i := pivot;
    loop
      Write_path(file,natural32(i));
      Corrector(file,s(i),c);
      sWrite(file,s(i).sol.all);
      ffail := (((s(i).resa > c.epsaf) and then (s(i).cora > c.epsax))
        or else ((s(i).resr > c.epsrf) and then (s(i).corr > c.epsrx)));
      if not ffail then
        j := 1;
        Equals(s,s(i).sol.v,i,dist_sols,j);
        if j /= i
         then ffail := true;
        end if;
      end if;
      exit when ffail;
      i := i + 1;
      if i > s'last
       then i := s'first;
      end if;
      exit when (i = pivot);
    end loop;
    if ffail
     then pivot := i;
    end if;
    fail := ffail;
  end Multiple_Reporting_Corrector;

-- TARGET PROCEDURES :

  procedure Affine_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : Floating_Number;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));       -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop    -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      if info /= 0
       then Clear(j); exit;               -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(s.cora); s.cora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(s.resa); s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(s.corr); s.corr := s.cora/normv;
            Clear(s.resr); s.resr := s.resa/normv;
      end if;
      Clear(normv);
      nit := nit + 1;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                       -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Affine_Single_Loose_Normal_Silent_Corrector;

  procedure Affine_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : Floating_Number;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));         -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop      -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      if info /= 0
       then Clear(j); exit;                 -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(s.cora); s.cora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(s.resa); s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(s.corr); s.corr := s.cora/normv;
            Clear(s.resr); s.resr := s.resa/normv;
      end if;
      Clear(normv);
      cWrite(file,s.cora,s.corr,s.resa,s.resr);       -- WRITE STATISTICS 
      nit := nit + 1;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                      -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Affine_Single_Loose_Normal_Reporting_Corrector;

  procedure Affine_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : Floating_Number;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));           -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufco(j,y'last,ipvt,s.rcond);
      if Equal(s.rcond,0.0)
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(s.cora); s.cora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(s.resa); s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(s.corr); s.corr := s.cora/normv;
            Clear(s.resr); s.resr := s.resa/normv;
      end if;
      Clear(normv);
      nit := nit + 1;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                       -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Affine_Single_Loose_Conditioned_Silent_Corrector;

  procedure Affine_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : Floating_Number;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufco(j,y'last,ipvt,s.rcond);
      if Equal(s.rcond,0.0)
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(s.cora); s.cora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(s.resa); s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(s.corr); s.corr := s.cora/normv;
            Clear(s.resr); s.resr := s.resa/normv;
      end if;
      Clear(normv);
      cWrite(file,s.cora,s.corr,s.resa,s.resr);        -- WRITE STATISTICS
      cWrite(file,s.rcond,natural32(s.sol.m));
      nit := nit + 1;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                       -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Affine_Single_Loose_Conditioned_Reporting_Corrector;

  procedure Affine_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : Floating_Number;
    stop : boolean;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      if info /= 0
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(ncora); ncora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(nresa); nresa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(ncorr); ncorr := ncora/normv;
            Clear(nresr); nresr := nresa/normv;
      end if;
      Clear(normv);
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
               or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                 -- STOP WHEN DIVERGENCE
      Copy(ncora,s.cora); Clear(ncora);
      Copy(nresa,s.resa); Clear(nresa);
      Copy(ncorr,s.corr); Clear(ncorr);
      Copy(nresr,s.resr); Clear(nresr);
      exit when stop;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                     -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Affine_Single_Severe_Normal_Silent_Corrector;

  procedure Affine_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : Floating_Number;
    stop : boolean;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      if info /= 0
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(ncora); ncora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(nresa); nresa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(ncorr); ncorr := ncora/normv;
            Clear(nresr); nresr := nresa/normv;
      end if;
      Clear(normv);
      cWrite(file,ncora,ncorr,nresa,nresr);           -- WRITE STATISTICS
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
              or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                  -- STOP WHEN DIVERGENCE
      Copy(ncora,s.cora); Clear(ncora);
      Copy(nresa,s.resa); Clear(nresa);
      Copy(ncorr,s.corr); Clear(ncorr); 
      Copy(nresr,s.resr); Clear(nresr);
      exit when stop;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                      -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Affine_Single_Severe_Normal_Reporting_Corrector;

  procedure Affine_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : Floating_Number;
    stop : boolean;

  begin
    Clear(s.resa); s.resa := Norm(y); 
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop        -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufco(j,y'last,ipvt,s.rcond);
      if Equal(s.rcond,0.0)
       then Clear(j); exit;                   -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(ncora); ncora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(nresa); nresa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(ncorr); ncorr := ncora/normv;
            Clear(nresr); nresr := nresa/normv;
      end if;
      Clear(normv);
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
            or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                    -- STOP WHEN DIVERGENCE
      Copy(ncora,s.cora); Clear(ncora);
      Copy(nresa,s.resa); Clear(nresa);
      Copy(ncorr,s.corr); Clear(ncorr);
      Copy(nresr,s.resr); Clear(nresr);
      exit when stop;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                        -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Affine_Single_Severe_Conditioned_Silent_Corrector;

  procedure Affine_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : Floating_Number;
    stop : boolean;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));           -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop        -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufco(j,y'last,ipvt,s.rcond);
      if Equal(s.rcond,0.0)
       then Clear(j); exit;                   -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(ncora); ncora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(nresa); nresa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(ncorr); ncorr := ncora/normv;
            Clear(nresr); nresr := nresa/normv;
      end if;
      Clear(normv);
      cWrite(file,ncora,ncorr,nresa,nresr);             -- WRITE STATISTICS
      cWrite(file,s.rcond,natural32(s.sol.m));
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
            or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                  -- STOP WHEN DIVERGENCE
      Copy(ncora,s.cora); Clear(ncora);
      Copy(nresa,s.resa); Clear(nresa);
      Copy(ncorr,s.corr); Clear(ncorr);
      Copy(nresr,s.resr); Clear(nresr);
      exit when stop;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                      -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Affine_Single_Severe_Conditioned_Reporting_Corrector;

  procedure Affine_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Normal_Silent_Corrector;

  procedure Affine_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Normal_Reporting_Corrector;

  procedure Affine_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Conditioned_Silent_Corrector;

  procedure Affine_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Conditioned_Reporting_Corrector;

  procedure Affine_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Normal_Silent_Corrector;

  procedure Affine_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Normal_Reporting_Corrector;

  procedure Affine_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Conditioned_Silent_Corrector;

  procedure Affine_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
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
    normv : Floating_Number;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        Copy(s.sol.v(jj),j(j'last(1),jj));         -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      if info /= 0
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      Clear(y(y'last)); y(y'last) := Create(integer(0));  -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(s.cora); s.cora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(y(y'last)); y(y'last) := Create(integer(0));
      Clear(s.resa); s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(s.corr); s.corr := s.cora/normv;
            Clear(s.resr); s.resr := s.resa/normv;
            for jj in s.sol.v'range loop            -- SCALE THE SOLUTION
              Div(s.sol.v(jj),normv);
            end loop;
      end if;
      Clear(normv);
      nit := nit + 1;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
         and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                       -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Projective_Single_Loose_Normal_Silent_Corrector;

  procedure Projective_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : Floating_Number;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        Copy(s.sol.v(jj),j(j'last(1),jj));         -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      if info /= 0
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      Clear(y(y'last)); y(y'last) := Create(integer(0));   -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(s.cora); s.cora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(y(y'last)); y(y'last) := Create(integer(0));
      Clear(s.resa); s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(s.corr); s.corr := s.cora/normv;
            Clear(s.resr); s.resr := s.resa/normv;
            for jj in s.sol.v'range loop             -- SCALE THE SOLUTION
              Div(s.sol.v(jj),normv);
            end loop;
      end if;
      Clear(normv);
      cWrite(file,s.cora,s.corr,s.resa,s.resr);        -- WRITE STATISTICS
      nit := nit + 1;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                       -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Projective_Single_Loose_Normal_Reporting_Corrector;

  procedure Projective_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : Floating_Number;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        Copy(s.sol.v(jj),j(j'last(1),jj));         -- CORRECT PERPENDICULAR
      end loop;
      lufco(j,y'last,ipvt,s.rcond);
      if Equal(s.rcond,0.0)
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      Clear(y(y'last)); y(y'last) := Create(integer(0));  -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(s.cora); s.cora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(y(y'last)); y(y'last) := Create(integer(0));
      Clear(s.resa); s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(s.corr); s.corr := s.cora/normv;
            Clear(s.resr); s.resr := s.resa/normv;
            for jj in s.sol.v'range loop              -- SCALE THE SOLUTION
              Div(s.sol.v(jj),normv);
            end loop;
      end if;
      Clear(normv);
      nit := nit + 1;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                        -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Projective_Single_Loose_Conditioned_Silent_Corrector;

  procedure Projective_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv : Floating_Number;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));           -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop        -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);             -- CORRECT PERPENDICULAR
      end loop;
      lufco(j,y'last,ipvt,s.rcond);
      if Equal(s.rcond,0.0)
       then Clear(j); exit;                   -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      Clear(y(y'last)); y(y'last) := Create(integer(0));   -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(s.cora); s.cora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(y(y'last)); y(y'last) := Create(integer(0));
      Clear(s.resa); s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(s.corr); s.corr := s.cora/normv;
            Clear(s.resr); s.resr := s.resa/normv;
            for jj in s.sol.v'range loop               -- SCALE THE SOLUTION
              Div(s.sol.v(jj),normv);
            end loop;
      end if;
      Clear(normv);
      cWrite(file,s.cora,s.corr,s.resa,s.resr);          -- WRITE STATISTICS
      cWrite(file,s.rcond,natural32(s.sol.m));
      nit := nit + 1;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                         -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Projective_Single_Loose_Conditioned_Reporting_Corrector;

  procedure Projective_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : Floating_Number;
    stop : boolean;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv; 
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);            -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      if info /= 0
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      Clear(y(y'last)); y(y'last) := Create(integer(0));  -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y); 
      Clear(ncora); ncora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(y(y'last)); y(y'last) := Create(integer(0));
      Clear(nresa); nresa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(ncorr); ncorr := ncora/normv;
            Clear(nresr); nresr := nresa/normv;
            for jj in s.sol.v'range loop              -- SCALE THE SOLUTION
              Div(s.sol.v(jj),normv);
            end loop;
      end if;
      Clear(normv);
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
             or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                    -- STOP WHEN DIVERGENCE
      Copy(ncora,s.cora); Clear(ncora);
      Copy(nresa,s.resa); Clear(nresa);
      Copy(ncorr,s.corr); Clear(ncorr);
      Copy(nresr,s.resr); Clear(nresr);
      exit when stop;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                        -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Projective_Single_Severe_Normal_Silent_Corrector;

  procedure Projective_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer32;
    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : Floating_Number;
    stop : boolean;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION 
    normv := Norm(s.sol.v);
    if Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        Copy(s.sol.v(jj),j(j'last(1),jj));         -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      if info /= 0
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      Clear(y(y'last)); y(y'last) := Create(integer(0));  -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(ncora); ncora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(y(y'last)); y(y'last) := Create(integer(0));
      Clear(nresa); nresa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(ncorr); ncorr := ncora/normv;
            Clear(nresr); nresr := nresa/normv;
            for jj in s.sol.v'range loop             -- SCALE THE SOLUTION
              Div(s.sol.v(jj),normv);
            end loop;
      end if;
      Clear(normv);
      cWrite(file,ncora,ncorr,nresa,nresr);            -- WRITE STATISTICS
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
            or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                   -- STOP WHEN DIVERGENCE
      Copy(ncora,s.cora); Clear(ncora);
      Copy(nresa,s.resa); Clear(nresa);
      Copy(ncorr,s.corr); Clear(ncorr);
      Copy(nresr,s.resr); Clear(nresr);
      exit when stop;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
         and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                      -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Projective_Single_Severe_Normal_Reporting_Corrector;

  procedure Projective_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : Floating_Number;
    stop : boolean;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));          -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop       -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        Copy(s.sol.v(jj),j(j'last(1),jj));         -- CORRECT PERPENDICULAR
      end loop;
      lufco(j,y'last,ipvt,s.rcond);
      if Equal(s.rcond,0.0)
       then Clear(j); exit;                  -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      Clear(y(y'last)); y(y'last) := Create(integer(0));  -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(ncora); ncora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      Clear(y(y'last)); y(y'last) := Create(integer(0));
      Clear(nresa); nresa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(ncorr); ncorr := ncora/normv;
            Clear(nresr); nresr := nresa/normv;
            for jj in s.sol.v'range loop              -- SCALE THE SOLUTION
              Div(s.sol.v(jj),normv);
            end loop;
      end if;
      Clear(normv);
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
             or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                    -- STOP WHEN DIVERGENCE
      Copy(ncora,s.cora); Clear(ncora);
      Copy(nresa,s.resa); Clear(nresa);
      Copy(ncorr,s.corr); Clear(ncorr);
      Copy(nresr,s.resr); Clear(nresr);
      exit when stop;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                     -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Projective_Single_Severe_Conditioned_Silent_Corrector;

  procedure Projective_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : Vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Standard_Integer_Vectors.Vector(y'range);
    nit : natural32 := 0;
    normv,ncora,nresa,ncorr,nresr : Floating_Number;
    stop : boolean;

  begin
    Clear(s.resa); s.resa := Norm(y);
    Clear(s.cora); s.cora := Create(integer(1));            -- INITIALIZATION
    normv := Norm(s.sol.v);
    if not Equal(normv,0.0)
     then Clear(s.corr); s.corr := s.cora/normv;
          Clear(s.resr); s.resr := s.resa/normv;
    end if;
    Clear(normv);
    while nit < c.maxit loop         -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        Copy(s.sol.v(jj),j(j'last(1),jj));           -- CORRECT PERPENDICULAR
      end loop;
      lufco(j,y'last,ipvt,s.rcond);
      if Equal(s.rcond,0.0)
       then Clear(j); exit;                    -- STOP WHEN SINGULAR JACOBIAN
      end if;
      Min(y);
      Clear(y(y'last)); y(y'last) := Create(integer(0));    -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); Clear(j);
      Add(s.sol.v,y);
      Clear(ncora); ncora := Norm(y); Clear(y);
      Add(s.length_path,s.cora);
      y := H(s.sol.v,s.sol.t);
      y(y'last) := Create(integer(0));
      Clear(nresa); nresa := Norm(y);
      normv := Norm(s.sol.v);
      if not Equal(normv,0.0)
       then Clear(ncorr); ncorr := ncora/normv;
            Clear(nresr); nresr := nresa/normv;
            for jj in s.sol.v'range loop                -- SCALE THE SOLUTION
              Div(s.sol.v(jj),normv);
            end loop;
      end if;
      Clear(normv);
      cWrite(file,ncora,ncorr,nresa,nresr);               -- WRITE STATISTICS
      cWrite(file,s.rcond,natural32(s.sol.m));
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
              or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                      -- STOP WHEN DIVERGENCE
      Copy(ncora,s.cora); Clear(ncora);
      Copy(nresa,s.resa); Clear(nresa);
      Copy(ncorr,s.corr); Clear(ncorr);
      Copy(nresr,s.resr); Clear(nresr);
      exit when stop;
      exit when (((s.cora < c.epsax) or else (s.resa < c.epsaf))
        and then ((s.corr < c.epsrx) or else (s.resr < c.epsrf)));
    end loop;                          -- STOP WHEN DESIRED PRECISION REACHED
    Clear(y);
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when constraint_error => return;
  end Projective_Single_Severe_Conditioned_Reporting_Corrector;

  procedure Projective_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Normal_Silent_Corrector;

  procedure Projective_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Normal_Reporting_Corrector;

  procedure Projective_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Conditioned_Silent_Corrector;

  procedure Projective_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Conditioned_Reporting_Corrector;

  procedure Projective_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Normal_Silent_Corrector;

  procedure Projective_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Normal_Reporting_Corrector;

  procedure Projective_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Conditioned_Silent_Corrector;

  procedure Projective_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out integer32; dist_sols : in Floating_Number;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Conditioned_Reporting_Corrector;

end Multprec_Correctors;
