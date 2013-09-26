with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;

package body DoblDobl_Continuation_Data_io is

  procedure Write_Statistics
               ( file : in file_type;
                 i,nstep,nfail,niter,nsyst : in natural32 ) is
  begin
    put(file,"== "); put(file,i,1); put(file," = ");
    put(file," #step : "); put(file,nstep,3);
    put(file," #fail : "); put(file,nfail,2);
    put(file," #iter : "); put(file,niter,3);
    if nsyst /= niter
     then put(file," #syst : "); put(file,nsyst,3);
    end if;
    put(file," = ");
  end Write_Statistics;

  procedure Write_Diagnostics
               ( file : in file_type; s : in Solution;
                 tol_zero,tol_sing : in double_float;
                 nbfail,nbregu,nbsing : in out natural32;
                 kind : out natural32 ) is
  begin
    if (abs(REAL_PART(s.t) - 1.0) > tol_zero)
      or else ((s.err > tol_zero) and (s.res > tol_zero)) then
      put_line(file,"failure");
      nbfail := nbfail + 1;
      kind := 0;
    elsif s.rco < tol_sing then
      put_line(file,"singular solution");
      nbsing := nbsing + 1;
      kind := 2;
    else
      put_line(file,"regular solution");
      nbregu := nbregu + 1;
      kind := 1;
    end if;
  end Write_Diagnostics;

  procedure Write_Diagnostics
               ( file : in file_type; s : in Solu_Info;
                 tol_zero,tol_sing : in double_float;
                 nbfail,nbregu,nbsing : in out natural32;
                 kind : out natural32 ) is
  begin
    if (abs(REAL_PART(s.sol.t) - 1.0) > tol_zero) or else
      ((s.cora > tol_zero) and (s.resa > tol_zero)) then
      put_line(file,"failure");
      nbfail := nbfail + 1;
      kind := 0;
    elsif s.rcond < tol_sing then
      put_line(file,"singular solution");
      nbsing := nbsing + 1;
      kind := 2;
    else
      put_line(file,"regular solution");
      nbregu := nbregu + 1;
      kind := 1;
    end if;
  end Write_Diagnostics;

  procedure Write_Solution ( file : in file_type; s : in Solution;
                             length_path : in double_double ) is
  begin
    put(file,"t : "); put(file,s.t); new_line(file);
    put(file,"m : "); put(file,s.m,1);
    put(file,"                  Length of path : ");
    put(file,length_path,3);
    new_line(file);
    put_line(file,"the solution for t : ");
    put_vector(file,s);
    put(file,"==");
    put(file," err : "); put(file,s.err,3);  put(file," =");
    put(file," rco : "); put(file,s.rco,3); put(file," =");
    put(file," res : "); put(file,s.res,3);  put_line(file," ==");
  end Write_Solution;

  procedure Write_Solution ( file : in file_type; s : in Solu_Info ) is
  begin
    put(file,"t : "); put(file,s.sol.t); new_line(file);
    put(file,"m : "); put(file,s.sol.m,1);
    put(file,"                  Length of path : ");
    put(file,s.length_path);
    new_line(file);
    put_line(file,"the solution for t : ");
    put_vector(file,s.sol.all);
    put(file,"==");
    put(file," err : "); put(file,s.cora,3);  put(file," =");
    put(file," rco : "); put(file,s.rcond,3); put(file," =");
    put(file," res : "); put(file,s.resa,3);  put_line(file," ==");
  end Write_Solution;

  procedure Write_Next_Solution
               ( file : in file_type;
                 cnt : in out natural32; s : in Solu_Info ) is
  begin
    cnt := cnt + 1;
    Write_Statistics(file,cnt,s.nstep,s.nfail,s.niter,s.nsyst);
    new_line(file);
    Write_Solution(file,s);
  end Write_Next_Solution;

  procedure Write_Next_Solution
               ( file : in file_type;
                 cnt : in out natural32; s : in Solu_Info;
                 tol_zero,tol_sing : in double_float;
                 nbfail,nbregu,nbsing : in out natural32;
                 kind : out natural32 ) is
  begin
    cnt := cnt + 1;
    Write_Statistics(file,cnt,s.nstep,s.nfail,s.niter,s.nsyst);
    Write_Diagnostics(file,s,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    Write_Solution(file,s);
  end Write_Next_Solution;

end DoblDobl_Continuation_Data_io;
