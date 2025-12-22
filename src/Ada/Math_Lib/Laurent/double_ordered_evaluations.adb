with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Floating_Vectors_IO;      use Standard_Floating_Vectors_IO;
with Double_Real_Powered_Series;
with Double_Leading_Evaluations;

package body Double_Ordered_Evaluations is

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
      put_line("first_derivative_first_order 1 ...");
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
    Double_Real_Powered_Series.Sort(ydg,ycf); 
    Double_Real_Powered_Series.Normalize(ycf,ydg);
  end First_Derivative_First_Order;

  procedure First_Derivative_Second_Order
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
    lc2 : Standard_Complex_Vectors.Vector(cff'range);
    pw1 : Standard_Floating_Vectors.Vector(pwr'range);
    pw2 : Standard_Floating_Vectors.Vector(pwr'range);
    idx : integer32 := 0;
    shared : Complex_Number;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_derivative_second_order 1 ...");
    end if;
    for i in lc0'range loop
      lc0(i) := cff(i)(cff(i)'first);
      lc1(i) := cff(i)(cff(i)'first+1);
      lc2(i) := cff(i)(cff(i)'first+2);
      pw1(i) := pwr(i)(pwr(i)'first);
      pw2(i) := pwr(i)(pwr(i)'first+1);
    end loop;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,lc0,0,vrblvl-1);
      for j in cff'range loop -- all first order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j);
        shared := pcf(i)*Leading_Coefficient(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc1(j);
        idx := idx + 1;              -- second order terms
        ydg(idx) := pct(i) + pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc2(j);
      end loop;
    end loop;
    if vrblvl > 0 then
      put("idx : "); put(idx,1); new_line;
      put("ycf'last : "); put(ycf'last,1); new_line;
    end if;
    Double_Real_Powered_Series.Sort(ydg,ycf); 
    Double_Real_Powered_Series.Normalize(ycf,ydg);
  end First_Derivative_Second_Order;

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
      put_line("second_derivative_first_order 1 ...");
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
      for j in cff'range loop -- all pure second derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*lpw(j);
        ycf(idx) := pcf(i)*Second_Derivative(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*lc1(j)*lc1(j)/create(2.0); -- Taylor series
      end loop;
      for j in cff'range loop -- all mixed second derivative terms
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
    Double_Real_Powered_Series.Sort(ydg,ycf); 
    Double_Real_Powered_Series.Normalize(ycf,ydg);
  end Second_Derivative_First_Order;

  procedure Second_Derivative_Second_Order
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
    lc2 : Standard_Complex_Vectors.Vector(cff'range);
    pw1 : Standard_Floating_Vectors.Vector(pwr'range);
    pw2 : Standard_Floating_Vectors.Vector(pwr'range);
    idx : integer32 := 0;
    shared : Complex_Number;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("second_derivative_second_order 1 ...");
    end if;
    for i in lc0'range loop
      lc0(i) := cff(i)(cff(i)'first);
      lc1(i) := cff(i)(cff(i)'first+1);
      lc2(i) := cff(i)(cff(i)'first+2);
      pw1(i) := pwr(i)(pwr(i)'first);
      pw2(i) := pwr(i)(pwr(i)'first+1);
    end loop;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,lc0,0,vrblvl-1);
      for j in cff'range loop -- all first order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j);
        shared := pcf(i)*Leading_Coefficient(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc1(j);
        idx := idx + 1;
        ydg(idx) := pct(i) + pw2(j); -- second order term
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc2(j);
      end loop;
      for j in cff'range loop -- all pure second derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*pw1(j);
        shared := pcf(i)*Second_Derivative(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc1(j)*lc1(j)/create(2.0); -- Taylor series
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc2(j)*lc2(j)/create(2.0); -- Taylor series
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j) + pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc1(j)*lc2(j); -- double compensates /2
      end loop;
      for j in cff'range loop -- all mixed second derivative terms
        for k in j+1..cff'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + pw1(j) + pw1(k);
          shared
            := pcf(i)*Second_Mixed_Derivative(pdg(i).all,lc0,j,k,vrblvl-1);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc1(j)*lc1(k); -- /create(2.0);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw2(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc2(j)*lc2(k); -- /create(2.0);
        end loop;
      end loop;
    end loop;
    if vrblvl > 0 then
      put("idx : "); put(idx,1); new_line;
      put("ycf'last : "); put(ycf'last,1); new_line;
    end if;
    Double_Real_Powered_Series.Sort(ydg,ycf); 
    Double_Real_Powered_Series.Normalize(ycf,ydg);
  end Second_Derivative_Second_Order;

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
      put_line("third_derivative_first_order 1 ...");
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
    Double_Real_Powered_Series.Sort(ydg,ycf); 
    Double_Real_Powered_Series.Normalize(ycf,ydg);
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
      put_line("first_derivative_first_order 2 ...");
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
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        psm(i) := ydg(idx);
        csm(i) := ycf(idx);
      end;
    end loop;
  end First_Derivative_First_Order;

  procedure First_Derivative_Second_Order
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
      put_line("first_derivative_second_order 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      declare
        ycf : Standard_Complex_Vectors.Vector(1..(2*dim+1)*nbr);
        ydg : Standard_Floating_Vectors.Vector(1..(2*dim+1)*nbr);
      begin
        First_Derivative_Second_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cff,pwr,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the first order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
             put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        psm(i) := ydg(idx);
        csm(i) := ycf(idx);
      end;
    end loop;
  end First_Derivative_Second_Order;

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
      put_line("second_derivative_first_order 2 ...");
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
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        psm(i) := ydg(idx);
        csm(i) := ycf(idx);
      end;
    end loop;
  end Second_Derivative_First_Order;

  procedure Second_Derivative_Second_Order
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
      put_line("second_derivative_second_order 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := (1 + 3*dim + dim*(dim+1))*nbr;
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Second_Derivative_Second_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cff,pwr,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the second order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
            put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        psm(i) := ydg(idx);
        csm(i) := ycf(idx);
      end;
    end loop;
  end Second_Derivative_Second_Order;

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
      put_line("third_derivative_first_order 2 ...");
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
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        psm(i) := ydg(idx);
        csm(i) := ycf(idx);
      end;
    end loop;
  end Third_Derivative_First_Order;

end Double_Ordered_Evaluations;
