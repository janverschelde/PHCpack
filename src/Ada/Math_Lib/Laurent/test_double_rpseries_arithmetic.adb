with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Test_Real_Powered_Series;

package body Test_Double_rpSeries_Arithmetic is

  function Random_Series
             ( size : integer32; vrblvl : integer32 := 0 )
             return Double_Real_Power_Series.Series is

    cff : Standard_Complex_Vectors.Vector(0..size);
    pwt : Standard_Floating_Vectors.Vector(1..size);

  begin
    if vrblvl > 0 then
      put("-> generating a series of size "); put(size,1); put_line(" ...");
    end if;
    Test_Real_Powered_Series.Random_Series(size,cff,pwt);
    if vrblvl > 0
     then Test_Real_Powered_Series.Write(cff,pwt);
    end if;
    return Double_Real_Power_Series.make(cff,pwt);
  end Random_Series;

  function Random_Series
             ( size : integer32; vrblvl : integer32 := 0 )
             return Double_Real_Power_Series.Link_to_Series is

    res_rep : constant Double_Real_Power_Series.Series(size)
            := Random_Series(size,vrblvl);
    res : constant Double_Real_Power_Series.Link_to_Series
        := new Double_Real_Power_Series.Series'(res_rep);

  begin
    return res;
  end Random_Series;

  procedure Write ( s : in Double_Real_Power_Series.Series ) is
  begin
    Test_Real_Powered_Series.Write(s.cff,s.pwt);
  end Write;

  procedure Write ( s : in Double_Real_Power_Series.Link_to_Series ) is
  begin
    Test_Real_Powered_Series.Write(s.cff,s.pwt);
  end Write;

  procedure SQRT_Newton_Step
              ( x,A : in Double_Real_Power_Series.Link_to_Series;
                fx,dx : in out Double_Real_Power_Series.Link_to_Series ) is

    two : constant complex_number := create(2.0);

  begin
    Double_Real_Power_Series.Copy(x,fx);
    Double_Real_Power_Series.Mul(fx,x);    -- fx = x^2
    Double_Real_Power_Series.Sub(fx,A);    -- fx = x^2 - A
    Double_Real_Power_Series.Copy(fx,dx);
    Double_Real_Power_Series.Div(dx,two);  -- dx = (x^2 - A)/2
    Double_Real_Power_Series.Div(dx,x);    -- dx = (x^2 - A)/(2*x)
  end SQRT_Newton_Step;

  procedure Test_SQRT ( size : in integer32 ) is

    rps : constant Double_Real_Power_Series.Link_to_Series
        := Random_Series(size);
    wrk,fx,dx : Double_Real_Power_Series.Link_to_Series;
    maxit : constant integer32 := 10;
    nbrit : integer32 := maxit;
    pwtwrk,pwtfx : double_float;
    done : boolean := false;

  begin
    put_line("a random power series :");
    Write(rps);
    Double_Real_Power_Series.Copy(rps,wrk);
    wrk.cff(0) := Standard_Complex_Numbers_Polar.Root(wrk.cff(0),2,1);
    for i in 1..maxit loop
      SQRT_Newton_Step(wrk,rps,fx,dx);
      Double_Real_Power_Series.Sub(wrk,dx);
      put("At step "); put(i,1); put_line(" dx :"); Write(dx);
      put("At step "); put(i,1); put_line(" fx :"); Write(fx);
      put("At step "); put(i,1); put_line(" x :"); Write(wrk);
      pwtwrk := wrk.pwt(wrk.pwt'last);
      pwtfx := fx.pwt(1);
      put(" last exponent x : "); put(pwtwrk); new_line;
      put("first exponent fx :"); put(pwtfx);
      if pwtwrk < pwtfx
       then put_line("  converged, done!"); nbrit := i; done := true;
       else put_line("  more steps needed ...");
      end if;
      exit when done;
    end loop;
    if done 
     then put("Converged in "); put(nbrit,1);
     else put("Did not converged after "); put(maxit,1);
    end if;
    put_line(" steps.");
  end Test_SQRT;

  procedure Main is

    size : integer32 := 0;

  begin
    new_line;
    put("Give the size of the series : "); get(size);
    Test_SQRT(size);
  end Main;

end Test_Double_rpSeries_Arithmetic;
