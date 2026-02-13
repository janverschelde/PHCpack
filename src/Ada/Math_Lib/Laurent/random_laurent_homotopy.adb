with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Double_Real_Powered_Series;
with Test_Real_Powered_Series;

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
                tpw : out Standard_Floating_VecVecs.VecVec;
                intpow : in boolean := false ) is
  begin
    for i in 1..nbp loop
      declare
        dpi : constant Standard_Integer_VecVecs.VecVec(1..nbm(i))
            := Random_Polynomial_Support(nbm(i),dim,low,upp);
        cfi : constant Standard_Complex_Vectors.Vector(1..nbm(i))
            := Standard_Random_Vectors.Random_Vector(1,nbm(i));
        cti : Standard_Floating_Vectors.Vector(1..nbm(i))
            := Double_Real_Powered_Series.Random_Leading_Powers(nbm(i));
        pwt : integer32;
      begin
        deg(i) := new Standard_Integer_VecVecs.VecVec'(dpi);
        cff(i) := new Standard_Complex_Vectors.Vector'(cfi);
        if intpow then
          for k in cti'range loop
            pwt := integer32(cti(k));
            if pwt = 0
             then pwt := 1;
            end if;
            cti(k) := double_float(pwt); 
          end loop;
        end if;
        tpw(i) := new Standard_Floating_Vectors.Vector'(cti);
      end;
    end loop;
  end Random_Laurent_System;

  procedure Random_Laurent_System
              ( nbp,dim,low,upp,size : in integer32;
                nbm : in Standard_Integer_Vectors.Vector;
                deg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : out Standard_Complex_VecVecs.Array_of_VecVecs;
                tpw : out Standard_Floating_VecVecs.Array_of_VecVecs;
                intpow : in boolean := false ) is
  begin
    for i in 1..nbp loop
      declare
        dpi : constant Standard_Integer_VecVecs.VecVec(1..nbm(i))
            := Random_Polynomial_Support(nbm(i),dim,low,upp);
        cfi : Standard_Complex_VecVecs.VecVec(1..nbm(i));
        pwi : Standard_Floating_VecVecs.VecVec(1..nbm(i));
      begin
        deg(i) := new Standard_Integer_VecVecs.VecVec'(dpi);
        for j in 1..nbm(i) loop
          declare
            jcf : Standard_Complex_Vectors.Vector(0..size);
            jpw : Standard_Floating_Vectors.Vector(1..size);
          begin
            Test_Real_Powered_Series.Random_Series(size,jcf,jpw);
            if intpow then
              for k in jpw'range loop
                jpw(k) := Standard_Floating_Numbers.create(k);
              end loop;
            end if;
            cfi(j) := new Standard_Complex_Vectors.Vector'(jcf);
            pwi(j) := new Standard_Floating_Vectors.Vector'(jpw);
          end;
        end loop;
        cff(i) := new Standard_Complex_VecVecs.VecVec'(cfi);
        tpw(i) := new Standard_Floating_VecVecs.VecVec'(pwi);
      end;
    end loop;
  end Random_Laurent_System;

  procedure Product_Homotopy_Polynomial
              ( pdg : in Standard_Integer_VecVecs.VecVec;
                pcf : in Standard_Complex_Vectors.Vector;
                ptp : in Standard_Floating_Vectors.Vector;
                scf : in Standard_Complex_VecVecs.VecVec;
                spw : in Standard_Floating_VecVecs.VecVec;
                idxfac : in integer32;
                hdg : out Standard_Integer_VecVecs.VecVec;
                hcf : out Standard_Complex_Vectors.Vector;
                htp : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := scf'last;
    nbt : constant integer32
        := spw(idxfac)'last;  -- number of powers in the series
    deg : Standard_Integer_Vectors.Vector(1..dim);
    idx : integer32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in Random_Laurent_Homotopy.");
      put_line("random_homotopy_polynomial ...");
    end if;
    for j in pdg'range loop -- consider j-th monomial
      if vrblvl > 0
       then put("processing monomial "); put(j,1); put_line(" ...");
      end if;
      deg := pdg(j).all;
      deg(idxfac) := deg(idxfac) + 1; -- multiply by x(i)
      idx := idx + 1;
      hdg(idx) := new Standard_Integer_Vectors.Vector'(deg);
      hcf(idx) := pcf(j);
      htp(idx) := ptp(j);
      deg(idxfac) := deg(idxfac) - 1; -- multiply by series
      idx := idx + 1;
      hdg(idx) := new Standard_Integer_Vectors.Vector'(deg);
      hcf(idx) := -scf(idxfac)(0)*pcf(j); -- constant of series
      htp(idx) := ptp(j);
      for k in 1..nbt loop -- run over all terms of series
        idx := idx + 1;
        hdg(idx) := new Standard_Integer_Vectors.Vector'(deg);
        hcf(idx) := -scf(idxfac)(k)*pcf(j); -- k-th term of series
        htp(idx) := spw(idxfac)(k) + ptp(j);
      end loop;
      if vrblvl > 0 
       then put("idx : "); put(idx,1); new_line;
      end if;
    end loop;
  end Product_Homotopy_Polynomial;

  procedure Product_Homotopy
              ( pdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pcf : in Standard_Complex_VecVecs.VecVec;
                ptp : in Standard_Floating_VecVecs.VecVec;
                scf : in Standard_Complex_VecVecs.VecVec;
                spw : in Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                htp : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    nbt,nbm,size : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in Random_Laurent_Homotopy.random_homotopy ...");
    end if;
    for i in pdg'range loop
      nbt := spw(i)'last;  -- number of powers in i-th series
      nbm := pdg(i)'last;  -- number of monomials in i-th polynomial
      size := (nbt+2)*nbm; -- size of i-th homotopy polynomial
      if vrblvl > 0 then
        put("nbt : "); put(nbt,1);
        put(", nbm : "); put(nbm,1); new_line;
        put("polynomial "); put(i,1);
        put(" has size "); put(size,1); new_line;
      end if;
      declare
        hdi : Standard_Integer_VecVecs.VecVec(1..size);
        hci : Standard_Complex_Vectors.Vector(1..size);
        hti : Standard_Floating_Vectors.Vector(1..size);
      begin
        Product_Homotopy_Polynomial
          (pdg(i).all,pcf(i).all,ptp(i).all,scf,spw,i,hdi,hci,hti,vrblvl-1);
        hdg(i) := new Standard_Integer_VecVecs.VecVec'(hdi);
        hcf(i) := new Standard_Complex_Vectors.Vector'(hci);
        htp(i) := new Standard_Floating_Vectors.Vector'(hti);
      end;
    end loop;
  end Product_Homotopy;

  procedure Canonical_Binomial
              ( idx,dim : in integer32;
                bdg : out Standard_Integer_VecVecs.VecVec;
                bcf : out Standard_Complex_Vectors.Vector;
                btp : out Standard_Floating_Vectors.Vector ) is

    mon : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);

  begin
    mon(idx) := 1;
    bdg(1) := new Standard_Integer_Vectors.Vector'(mon);
    mon(idx) := 0;
    bdg(2) := new Standard_Integer_Vectors.Vector'(mon);
    bcf(1) := create(1.0);
    bcf(2) := create(-1.0);
    btp(1) := 0.0;
    btp(2) := 0.0;
  end Canonical_Binomial;

  procedure Binomial_Homotopy
              ( pdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pcf : in Standard_Complex_VecVecs.VecVec;
                ptp : in Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                htp : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    bdg : Standard_Integer_VecVecs.VecVec(1..2);
    bcf : Standard_Complex_Vectors.Vector(1..2);
    btp : Standard_Floating_Vectors.Vector(1..2);
    dim : constant integer32 := pdg'last;
    size : integer32;

  begin
    for i in pdg'range loop
      Canonical_Binomial(i,dim,bdg,bcf,btp);
      size := pdg(i)'last+2;
      if vrblvl > 0 then
        put("nbm : "); put(pdg(i)'last,1);
        put(", polynomial "); put(i,1);
        put(" has size "); put(size,1); new_line;
      end if;
      declare
        hdi : Standard_Integer_VecVecs.VecVec(1..size);
        hci : Standard_Complex_Vectors.Vector(1..size);
        hti : Standard_Floating_Vectors.Vector(1..size);
      begin
        hdi(1) := bdg(1); hdi(2) := bdg(2);
        hci(1) := bcf(1); hci(2) := bcf(2);
        hti(1) := btp(1); hti(2) := btp(2);
        for j in pdg(i)'range loop
          hdi(2+j) := pdg(i)(j);
          hci(2+j) := pcf(i)(j);
          hti(2+j) := ptp(i)(j);
        end loop;
        hdg(i) := new Standard_Integer_VecVecs.VecVec'(hdi);
        hcf(i) := new Standard_Complex_Vectors.Vector'(hci);
        htp(i) := new Standard_Floating_Vectors.Vector'(hti);
      end;
    end loop;
  end Binomial_Homotopy;

  procedure Scale_Homotopy_Powers
              ( hct : in out Standard_Floating_Vectors.Vector ) is

    minpwr : double_float := hct(hct'first);

  begin
    for j in hct'first+1..hct'last loop
      if hct(j) < minpwr
       then minpwr := hct(j);
      end if;
    end loop;
    for j in hct'range loop
       hct(j) := hct(j) - minpwr;
    end loop;
  end Scale_Homotopy_Powers;

  procedure Scale_Homotopy_Powers
              ( hct : in Standard_Floating_VecVecs.VecVec ) is
  begin
    for i in hct'range loop
      Scale_Homotopy_Powers(hct(i).all);
    end loop;
  end Scale_Homotopy_Powers;

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

  procedure Test_Product_Homotopy
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
  end Test_Product_Homotopy;

end Random_Laurent_Homotopy;
