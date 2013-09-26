with Standard_Natural_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Handle_Underflow_Gracefully;        use Handle_Underflow_Gracefully;
with Process_io;                         use Process_io;

package body Standard_Monomial_Correctors is

  procedure Affine_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    vx : constant Vector := V(s.sol.v,s.sol.t);
    y : Vector(s.sol.v'range) := F(vx);
    A : Matrix(y'range,y'range);
    ipvt : Standard_Natural_Vectors.Vector(y'range);
    info : integer;
    nit : natural := 0;
    ncora,nresa,rcorr,nresr : double_float;

  begin
    while nit < c.maxit loop            -- stop when max #iterations reached
      A := J(vx);
      lufac(A,A'last(1),ipvt,info);
      exit when (info /= 0);            -- stop when singular Jacobian
      Min(y);
      lusolve(A,A'last(1),ipvt,y); ncora := Norm(y);
      Add(s.sol.v,y); s.length_path := s.length_path + s.cora;
    end loop;
  exception 
    when others => Underflow_to_Zero(s.sol.v); return;
  end Affine_Single_Severe_Normal_Silent_Corrector;

  procedure Affine_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Affine_Single_Severe_Conditioned_Silent_Corrector;

  procedure Affine_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Affine_Single_Severe_Normal_Reporting_Corrector;

  procedure Affine_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is
  begin
    null;
  end Affine_Single_Severe_Conditioned_Reporting_Corrector;

end Standard_Monomial_Correctors;
