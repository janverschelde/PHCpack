with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_QR_Least_Squares;  use DoblDobl_Complex_QR_Least_Squares;
with DoblDobl_Complex_Singular_Values;   use DoblDobl_Complex_Singular_Values;
--with Handle_Underflow_Gracefully;        use Handle_Underflow_Gracefully;
with Process_io;

package body DoblDobl_Orthogonal_Correctors is

  procedure Silent_QRLS_Corrector
              ( n : in integer32; s : in out Solu_Info; c : in Corr_Pars ) is

    nit : natural32 := 0;
    y : Vector(1..n) := H(s.sol.v,s.sol.t);
    j : Matrix(1..n,s.sol.v'range);
    m : constant integer32 := s.sol.v'last;
    dx : Vector(1..m);
    rsd,dum,dum2,dum3 : Vector(1..n);
    qraux : Vector(1..m) := (1..m => Create(integer(0)));
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    info : integer32;
    normv : double_float;

  begin
    while nit < c.maxit loop
      j := dH(s.sol.v,s.sol.t);
      QRD(j,qraux,jpvt,false);
      Min(y);
      QRLS(j,n,m,qraux,y,dum,dum2,dx,rsd,dum3,110,info);
      Add(s.sol.v,dx);
      s.cora := hi_part(Norm(dx)); s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      s.resa := hi_part(Norm(y));
      normv := hi_part(Norm(s.sol.v));
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  --exception -- suppress all underflow errors ...
  --  when others => Underflow_to_Zero(s.sol.v); return; 
  end Silent_QRLS_Corrector;

  procedure Reporting_QRLS_Corrector
              ( file : in file_type;
                n : in integer32; s : in out Solu_Info; c : in Corr_Pars ) is

    nit : natural32 := 0;
    y : Vector(1..n) := H(s.sol.v,s.sol.t);
    j : Matrix(1..n,s.sol.v'range);
    m : constant integer32 := s.sol.v'last;
    dx : Vector(1..m);
    rsd,dum,dum2,dum3 : Vector(1..n);
    qraux : Vector(1..m) := (1..m => Create(integer(0)));
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    info : integer32;
    normv : double_float;

  begin
    while nit < c.maxit loop
      j := dH(s.sol.v,s.sol.t);
      QRD(j,qraux,jpvt,false);
      Min(y);
      QRLS(j,n,m,qraux,y,dum,dum2,dx,rsd,dum3,110,info);
      Add(s.sol.v,dx);
      s.cora := hi_part(Norm(dx)); s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      s.resa := hi_part(Norm(y));
      normv := hi_part(Norm(s.sol.v));
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      Process_io.cWrite(file,s.cora,s.corr,s.resa,s.resr);
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  --exception -- suppress all underflow errors ...
  --  when others => Underflow_to_Zero(s.sol.v); return; 
  end Reporting_QRLS_Corrector;

  procedure Silent_SVD_Corrector
              ( n : in integer32; s : in out Solu_Info; c : in Corr_Pars ) is

    nit : natural32 := 0;
    y : Vector(1..n) := H(s.sol.v,s.sol.t);
    jac : Matrix(1..n,s.sol.v'range);
    p : constant integer32 := s.sol.v'last;
    dx : Vector(1..p);
    mm : constant integer32 := DoblDobl_Complex_Singular_Values.Min0(n+1,p);
    sv : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;
    normv : double_float;

  begin
    while nit < c.maxit loop
      jac := dH(s.sol.v,s.sol.t);
      SVD(jac,n,p,sv,e,u,v,job,info);
      s.rcond := hi_part(Inverse_Condition_Number(sv));
      Min(y);
      dx := Solve(u,v,sv,y);
      Add(s.sol.v,dx);
      s.cora := hi_part(Norm(dx)); s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      s.resa := hi_part(Norm(y));
      normv := hi_part(Norm(s.sol.v));
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  --exception -- suppress all underflow errors ...
  --  when others => Underflow_to_Zero(s.sol.v); return; 
  end Silent_SVD_Corrector;

  procedure Reporting_SVD_Corrector
              ( file : in file_type;
                n : in integer32; s : in out Solu_Info; c : in Corr_Pars ) is

    nit : natural32 := 0;
    y : Vector(1..n) := H(s.sol.v,s.sol.t);
    jac : Matrix(1..n,s.sol.v'range);
    p : constant integer32 := s.sol.v'last;
    dx : Vector(1..p);
    mm : constant integer32 := DoblDobl_Complex_Singular_Values.Min0(n+1,p);
    sv : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;
    normv : double_float;

  begin
    while nit < c.maxit loop
      jac := dH(s.sol.v,s.sol.t);
      SVD(jac,n,p,sv,e,u,v,job,info);
      s.rcond := hi_part(Inverse_Condition_Number(sv));
      Min(y);
      dx := Solve(u,v,sv,y);
      Add(s.sol.v,dx);
      s.cora := hi_part(Norm(dx)); s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      s.resa := hi_part(Norm(y));
      normv := hi_part(Norm(s.sol.v));
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      Process_io.cWrite(file,s.cora,s.corr,s.resa,s.resr);
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  --exception -- suppress all underflow errors ...
  --  when others => Underflow_to_Zero(s.sol.v); return; 
  end Reporting_SVD_Corrector;

end DoblDobl_Orthogonal_Correctors;
