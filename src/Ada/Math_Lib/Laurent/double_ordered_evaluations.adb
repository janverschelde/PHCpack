with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Floating_Vectors_IO;      use Standard_Floating_Vectors_IO;
with Double_Leading_Evaluations;

package body Double_Ordered_Evaluations is

  procedure Normalize
              ( cf : in out Standard_Complex_Vectors.Vector;
                dg : in Standard_Floating_Vectors.Vector ) is

    tol : constant double_float := 1.0E-12;
    dif : double_float;

  begin
    for i in cf'first..cf'last-1 loop
      for j in i+1..cf'last loop
        dif := abs(dg(i) - dg(j));
        exit when (dif > tol);
        cf(i) := cf(i) + cf(j);
        cf(j) := create(0.0);
      end loop;
    end loop;
  end Normalize;

  procedure First_Order_Evaluation
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    lc0 : Standard_Complex_Vectors.Vector(cff'range);
    lc1 : Standard_Complex_Vectors.Vector(cff'range);
    lpw : Standard_Floating_Vectors.Vector(pwr'range);
    idx : integer32 := 0;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_order_evaluation ...");
    end if;
    for i in lc0'range loop
      lc0(i) := cff(i)(cff(i)'first);
      lc1(i) := cff(i)(cff(i)'first+1);
      lpw(i) := pwr(i)(pwr(i)'first);
    end loop;
    if vrblvl > 0
     then put_line("lpw :"); put_line(lpw);
    end if;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,lc0,0,vrblvl-1);
      for j in cff'range loop -- all first order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + lpw(j);
        ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*lc1(j);
      end loop;
    end loop;
    Sort(ydg,ycf); 
    Normalize(ycf,ydg);
  end First_Order_Evaluation;

  procedure Second_Order_Evaluation
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    lc0 : Standard_Complex_Vectors.Vector(cff'range);
    lc1 : Standard_Complex_Vectors.Vector(cff'range);
    lpw : Standard_Floating_Vectors.Vector(pwr'range);
    idx : integer32 := 0;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("second_order_evaluation ...");
    end if;
    for i in lc0'range loop
      lc0(i) := cff(i)(cff(i)'first);
      lc1(i) := cff(i)(cff(i)'first+1);
      lpw(i) := pwr(i)(pwr(i)'first);
    end loop;
    if vrblvl > 0 then
      put_line("lpw : ");
      put_line(lpw);
    end if;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,lc0,0,vrblvl-1);
      for j in cff'range loop -- all first order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + lpw(j);
        ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*lc1(j);
      end loop;
      for j in cff'range loop -- all pure second order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*lpw(j);
        ycf(idx) := pcf(i)*Second_Derivative(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*lc1(j)*lc1(j)/create(2.0); -- Taylor series
      end loop;
      for j in cff'range loop -- all mixed second order terms
        for k in j+1..cff'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + lpw(j) + lpw(k);
          ycf(idx)
            := pcf(i)*Second_Mixed_Derivative(pdg(i).all,lc0,j,k,vrblvl-1);
          ycf(idx) := ycf(idx)*lc1(j)*lc1(k); -- /create(2.0);
        end loop;
      end loop;
    end loop;
    if vrblvl > 0 then
      put_line("Before sorting and normalizing ...");
      for i in ycf'range loop
        put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
      end loop;
      put("idx : "); put(idx,1); new_line;
      put("ycf'last : "); put(ycf'last,1); new_line;
    end if;
    Sort(ydg,ycf); 
    Normalize(ycf,ydg);
  end Second_Order_Evaluation;

end Double_Ordered_Evaluations;
