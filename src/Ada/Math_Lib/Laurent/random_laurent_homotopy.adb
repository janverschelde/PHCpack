with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Double_Real_Powered_Series;

package body Random_Laurent_Homotopy is

  function Random_Monomial_Support
             ( dim,low,upp : integer32 )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..dim);

  begin
    for i in 1..dim loop
      res(i) := Standard_Random_Numbers.Random(low,upp);
    end loop;
    return res;
  end Random_Monomial_Support;

  function Random_Polynomial_Support
             ( nbr,dim,low,upp : integer32 )
             return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(1..nbr);

  begin
    for i in 1..nbr loop
      declare
        mon : constant Standard_Integer_Vectors.Vector(1..dim)
            := Random_Monomial_Support(dim,low,upp);
      begin
        res(i) := new Standard_Integer_Vectors.Vector'(mon);
      end;
    end loop;
    return res;
  end Random_Polynomial_Support;

  procedure Random_Laurent_System
              ( nbp,dim,low,upp : in integer32;
                nbm : in Standard_Integer_Vectors.Vector;
                deg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : out Standard_Complex_VecVecs.VecVec;
                tpw : out Standard_Floating_VecVecs.VecVec ) is
  begin
    for i in 1..nbp loop
      declare
        dpi : constant Standard_Integer_VecVecs.VecVec(1..nbm(i))
            := Random_Polynomial_Support(nbm(i),dim,low,upp);
        cfi : constant Standard_Complex_Vectors.Vector(1..nbm(i))
            := Standard_Random_Vectors.Random_Vector(1,nbm(i));
        cti : constant Standard_Floating_Vectors.Vector(1..nbm(i))
            := Double_Real_Powered_Series.Random_Leading_Powers(nbm(i));
      begin
        deg(i) := new Standard_Integer_VecVecs.VecVec'(dpi);
        cff(i) := new Standard_Complex_Vectors.Vector'(cfi);
        tpw(i) := new Standard_Floating_Vectors.Vector'(cti);
      end;
    end loop;
  end Random_Laurent_System;

  procedure Random_Homotopy
              ( pdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pcf : in Standard_Complex_VecVecs.VecVec;
                ptp : in Standard_Floating_VecVecs.VecVec;
                scf : in Standard_Complex_VecVecs.VecVec;
                spw : in Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                htp : out Standard_Floating_VecVecs.VecVec ) is

    dim : constant integer32 := pdg'last;
    nbt,nbm,size,idx : integer32;

  begin
    for i in pdg'range loop
      nbt := spw(i)'last;  -- number of powers in i-th series
      nbm := pdg(i)'last;  -- number of monomials in i-th polynomial
      size := (nbt+2)*nbm; -- size of i-th homotopy polynomial
      declare
        hdi : Standard_Integer_VecVecs.VecVec(1..size);
        hci : Standard_Complex_Vectors.Vector(1..size);
        hti : Standard_Floating_Vectors.Vector(1..size);
      begin
        idx := 0;
        for j in pdg(i)'range loop -- consider j-th monomial
          declare
            deg : Standard_Integer_Vectors.Vector(1..dim) := pdg(i)(j).all;
          begin
            deg(i) := deg(i) + 1; -- multiply by x(i)
            idx := idx + 1;
            hdi(idx) := new Standard_Integer_Vectors.Vector'(deg);
            hci(idx) := pcf(i)(j);
            hti(idx) := ptp(i)(j);
            deg(i) := deg(i) - 1; -- multiply by series
            idx := idx + 1;
            hdi(idx) := new Standard_Integer_Vectors.Vector'(deg);
            hci(idx) := -scf(i)(0)*pcf(i)(j); -- constant of series
            hti(idx) := ptp(i)(j);
            for k in 1..nbt loop -- run over all terms of series
              idx := idx + 1;
              hdi(idx) := new Standard_Integer_Vectors.Vector'(deg);
              hci(idx) := -scf(i)(k)*pcf(i)(j); -- k-th term of series
              hti(idx) := spw(i)(k) + ptp(i)(j);
            end loop;
          end;
        end loop;
        hdg(i) := new Standard_Integer_VecVecs.VecVec'(hdi);
        hcf(i) := new Standard_Complex_Vectors.Vector'(hci);
        htp(i) := new Standard_Floating_Vectors.Vector'(hti);
      end;
    end loop;
  end Random_Homotopy;

  function Evaluate_Monomial
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Complex_Number; tpw : double_float;
               zpt : Standard_Complex_Vectors.Vector; tpt : double_float;
               vrblvl : integer32 := 0 ) return Complex_Number is

  -- DESCRIPTION :
  --   Evaluates the monomial with exponents in deg, coefficient in cff,
  --   and leading power of t in tpw at point zpt and parameter tpt.

    pwr : constant double_float := tpt**tpw;
    res : Complex_Number := pwr*cff;

  begin
    if vrblvl > 0 then
      put_line("-> in Test_Leading_Evaluations.evaluate_monomial ...");
      put("deg : "); put(deg);
      put(", cff : "); put(cff); new_line;
      put("zpt : "); put(zpt); new_line;
      put("tpw : "); put(tpw); new_line;
      put("tpt : "); put(tpt); new_line;
    end if;
    for i in zpt'range loop
      if deg(i) > 0 then
        for j in 1..deg(i) loop
          res := res*zpt(i);
        end loop;
      elsif deg(i) < 0 then
        for j in 1..(-deg(i)) loop
          res := res/zpt(i);
        end loop;
      end if;
    end loop;
    if vrblvl > 0 then
      put("res : "); put(res); new_line;
    end if;
    return res;
  end Evaluate_Monomial;

  function Evaluate_Homotopy
             ( deg : Standard_Integer_VecVecs.Array_of_VecVecs;
               cff : Standard_Complex_VecVecs.VecVec;
               tpw : Standard_Floating_VecVecs.VecVec;
               zpt : Standard_Complex_Vectors.Vector; tpt : double_float )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(cff'range);

  begin
    for i in res'range loop
      res(i) := create(0.0);
      for j in cff(i)'range loop
        res(i) := res(i)
          + Evaluate_Monomial(deg(i)(j).all,cff(i)(j),tpw(i)(j),zpt,tpt);
      end loop;
    end loop;
    return res;
  end Evaluate_Homotopy;

  procedure Test_Random_Homotopy
              ( hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : in Standard_Complex_VecVecs.VecVec;
                htp : in Standard_Floating_VecVecs.VecVec;
                scf : in Standard_Complex_VecVecs.VecVec;
                spw : in Standard_Floating_VecVecs.VecVec ) is

    tpt : constant double_float := abs(Standard_Random_Numbers.Random);
    zpt : constant Standard_Complex_Vectors.Vector(scf'range)
        := Double_Real_Powered_Series.Evaluate_Series(scf,spw,tpt);
    hpt : constant Standard_Complex_Vectors.Vector(hcf'range)
        := Evaluate_Homotopy(hdg,hcf,htp,zpt,tpt);

  begin
    put("random t : "); put(tpt); new_line;
    put_line("evaluated series :"); put_line(zpt);
    put_line("evaluated homotopy :"); put_line(hpt);
  end Test_Random_Homotopy;

end Random_Laurent_Homotopy;
