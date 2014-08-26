with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;

package body QuadDobl_Continuation_Data_io is

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
                             length_path : in quad_double ) is
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

  procedure Path_Tracking_Statistics
               ( file : in file_type; s : in Solu_Info_Array ) is

    min_nstep,max_nstep,min_nfail,max_nfail,min_niter,max_niter : natural32;
    avg_nstep,avg_nfail,avg_niter : double_float;
    len : constant double_float := double_float(s'last);

  begin
    min_nstep := s(s'first).nstep; max_nstep := s(s'first).nstep;
    min_nfail := s(s'first).nfail; max_nfail := s(s'first).nfail;
    min_niter := s(s'first).niter; max_niter := s(s'first).niter;
    avg_nstep := double_float(min_nstep)/len;
    avg_nfail := double_float(min_nfail)/len;
    avg_niter := double_float(min_niter)/len;
    for i in s'first+1..s'last loop
      if s(i).nstep < min_nstep then
        min_nstep := s(i).nstep;
      elsif s(i).nstep > max_nstep then
        max_nstep := s(i).nstep;
      end if;
      if s(i).nfail < min_nfail then
        min_nstep := s(i).nfail;
      elsif s(i).nfail > max_nfail then
        max_nfail := s(i).nfail;
      end if;
      if s(i).niter < min_niter then
        min_niter := s(i).niter;
      elsif s(i).niter > max_niter then
        max_niter := s(i).niter;
      end if;
      avg_nstep := avg_nstep + double_float(s(i).nstep)/len;
      avg_nfail := avg_nfail + double_float(s(i).nfail)/len;
      avg_niter := avg_niter + double_float(s(i).niter)/len;
    end loop;
    put_line(file,"path tracking statistics :");
    put(file,"           min       max       avg"); new_line(file);
    put(file,"#step : "); put(file,min_nstep,6); put(file,max_nstep,10);
    put(file,avg_nstep,8,3,0); new_line(file);
    put(file,"#fail : "); put(file,min_nfail,6); put(file,max_nfail,10);
    put(file,avg_nfail,8,3,0); new_line(file);
    put(file,"#iter : "); put(file,min_niter,6); put(file,max_niter,10);
    put(file,avg_niter,8,3,0); new_line(file);
  end Path_Tracking_Statistics;

end QuadDobl_Continuation_Data_io;
