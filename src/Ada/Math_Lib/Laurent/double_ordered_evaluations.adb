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

  function Positive_Minimum_Index
             ( c : Standard_Complex_Vectors.Vector;
               v : Standard_Floating_Vectors.Vector ) return integer32 is

    res,idx : integer32;
    psm : double_float;
    tol : constant double_float := 1.0E-12;

  begin
    for i in v'range loop -- find first positive number
      if AbsVal(c(i)) > tol then
        if v(i) > tol then
          idx := i;
          psm := v(i);
          exit;
        end if;
      end if;
    end loop;
    res := idx;
    for i in idx+1..v'last loop
      if AbsVal(c(i)) > tol then
        if v(i) > tol and then v(i) < psm
         then psm := v(i); res := i;
        end if;
      end if;
    end loop;
    return res;
  end Positive_Minimum_Index;

  procedure First_Derivative_First_Order
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
  end First_Derivative_First_Order;

  procedure Second_Derivative_First_Order
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
  end Second_Derivative_First_Order;

  procedure Third_Derivative_First_Order
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
      put_line("third_order_evaluation ...");
    end if;
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
          ycf(idx) := pcf(i)*Second_Mixed_Derivative
                               (pdg(i).all,lc0,j,k,vrblvl-1);
          ycf(idx) := ycf(idx)*lc1(j)*lc1(k); -- /create(2.0);
        end loop;
      end loop;
      for j in cff'range loop -- all pure third order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 3.0*lpw(j);
        ycf(idx) := pcf(i)*Third_Derivative(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*lc1(j)*lc1(j)/create(6.0); -- Taylor series
      end loop;
      for j in cff'range loop -- all semi mixed third order terms
        for k in j+1..cff'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + lpw(j) + lpw(k);
          ycf(idx) := pcf(i)*Third_Semi_Mixed_Derivative
                               (pdg(i).all,lc0,j,k,vrblvl-1);
          ycf(idx) := ycf(idx)*lc1(j)*lc1(k)/create(2.0);
          idx := idx + 1; -- flip role of j and k
          ydg(idx) := pct(i) + lpw(k) + lpw(j);
          ycf(idx) := pcf(i)*Third_Semi_Mixed_Derivative
                               (pdg(i).all,lc0,k,j,vrblvl-1);
          ycf(idx) := ycf(idx)*lc1(k)*lc1(j)/create(2.0);
        end loop;
      end loop;
      for j in cff'range loop -- all fully mixed third order terms
        for k in j+1..cff'last loop
          for L in k+1..cff'last loop
            idx := idx + 1;
            ydg(idx) := pct(i) + lpw(j) + lpw(k) + lpw(L);
            ycf(idx) := pcf(i)*Third_Fully_Mixed_Derivative
                                 (pdg(i).all,lc0,j,k,L,vrblvl-1);
            ycf(idx) := ycf(idx)*lc1(j)*lc1(k)*lc1(L);
          end loop;
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
  end Third_Derivative_First_Order;

-- ON A POLYNOMIAL HOMOTOPY :

  procedure First_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    nbr,idx : integer32;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_order_evaluation 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      declare
        ycf : Standard_Complex_Vectors.Vector(1..(dim+1)*nbr);
        ydg : Standard_Floating_Vectors.Vector(1..(dim+1)*nbr);
      begin
        First_Derivative_First_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cff,pwr,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the first order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
             put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Positive_Minimum_Index(ycf,ydg);
        psm(i) := ydg(idx);
        csm(i) := ycf(idx);
      end;
    end loop;
  end First_Derivative_First_Order;

  procedure Second_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    nbr,size,idx : integer32;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("second_order_evaluation 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := (1 + dim + dim*(dim+1)/2)*nbr;
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Second_Derivative_First_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cff,pwr,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the second order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
            put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Positive_Minimum_Index(ycf,ydg);
        psm(i) := ydg(idx);
        csm(i) := ycf(idx);
      end;
    end loop;
  end Second_Derivative_First_Order;

  procedure Third_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                psm : out Standard_Floating_Vectors.Vector;
                csm : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    nbr,size,idx : integer32;
    dd1 : constant integer32 := dim*(dim-1);

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("third_order_evaluation 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := (1 + dim + dim*(dim+1)/2 + dim + dd1 + dd1*(dim-2)/6)*nbr;
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Third_Derivative_First_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cff,pwr,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the third order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
            put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Positive_Minimum_Index(ycf,ydg);
        psm(i) := ydg(idx);
        csm(i) := ycf(idx);
      end;
    end loop;
  end Third_Derivative_First_Order;

end Double_Ordered_Evaluations;
