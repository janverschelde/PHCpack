with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
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
    for i in lc0'range loop
      lc0(i) := cff(i)(cff(i)'first);
      lc1(i) := cff(i)(cff(i)'first+1);
      lpw(i) := pwr(i)(pwr(i)'first);
    end loop;
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

end Double_Ordered_Evaluations;
