with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Double_Real_Powered_Series;
with Double_Leading_Evaluations;

package body Double_Ordered_Evaluations is

  function Size_Evaluation
             ( dim,ndf,ord,nbr : integer32 ) return integer32 is

    res : integer32 := -1;

  begin
    if ord = 1 then
      if ndf = 1 then -- first derivative first order
        res := (dim+1)*nbr;
      elsif ndf = 2 then -- second derivative first order
        res := (1 + dim + dim*(dim+1)/2)*nbr;
      elsif ndf = 3 then -- third derivative first order
        res := 1 + dim + dim*(dim+1)/2 
                 + dim + dim*(dim-1) + dim*(dim-1)*(dim-2)/6;
        res := res*nbr;
      end if;
    end if;
    return res;
  end Size_Evaluation;

  procedure First_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_derivative_first_order 0 ...");
    end if;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,0,vrblvl-1);
      for j in cf0'range loop -- all first order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j);
        ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*cf1(j);
      end loop;
    end loop;
    Double_Real_Powered_Series.Sort(ydg,ycf); 
    Double_Real_Powered_Series.Normalize(ycf,ydg);
  end First_Derivative_First_Order;

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
    First_Derivative_First_Order(pcf,pct,pdg,lc0,lc1,lpw,ycf,ydg,vrblvl);
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
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("second_derivative_first_order 0 ...");
    end if;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,0,vrblvl-1);
      for j in cf0'range loop -- all first order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j);
        ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*cf1(j);
      end loop;
      for j in cf0'range loop -- all pure second derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*pw1(j);
        ycf(idx) := pcf(i)*Second_Derivative(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*cf1(j)*cf1(j)/create(2.0); -- Taylor series
      end loop;
      for j in cf0'range loop -- all mixed second derivative terms
        for k in j+1..cf0'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + pw1(j) + pw1(k);
          ycf(idx)
            := pcf(i)*Second_Mixed_Derivative(pdg(i).all,cf0,j,k,vrblvl-1);
          ycf(idx) := ycf(idx)*cf1(j)*cf1(k);
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
    Second_Derivative_First_Order(pcf,pct,pdg,lc0,lc1,lpw,ycf,ydg,vrblvl);
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
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("third_derivative_first_order 0 ...");
    end if;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,0,vrblvl-1);
      for j in cf0'range loop -- all first order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j);
        ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*cf1(j);
      end loop;
      for j in cf0'range loop -- all pure second derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*pw1(j);
        ycf(idx) := pcf(i)*Second_Derivative(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*cf1(j)*cf1(j)/create(2.0); -- Taylor series
      end loop;
      for j in cf0'range loop -- all mixed second derivative terms
        for k in j+1..cf0'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + pw1(j) + pw1(k);
          ycf(idx) := pcf(i)*Second_Mixed_Derivative
                               (pdg(i).all,cf0,j,k,vrblvl-1);
          ycf(idx) := ycf(idx)*cf1(j)*cf1(k);
        end loop;
      end loop;
      for j in cf0'range loop -- all pure third derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 3.0*pw1(j);
        ycf(idx) := pcf(i)*Third_Derivative(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := ycf(idx)*cf1(j)*cf1(j)/create(6.0); -- Taylor series
      end loop;
      for j in cf0'range loop -- all semi mixed third derivative terms
        for k in j+1..cf0'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + 2.0*pw1(j) + pw1(k);
          ycf(idx) := pcf(i)*Third_Semi_Mixed_Derivative
                               (pdg(i).all,cf0,j,k,vrblvl-1);
          ycf(idx) := ycf(idx)*cf1(j)*cf1(k)/create(2.0);
          idx := idx + 1; -- flip role of j and k
          ydg(idx) := pct(i) + 2.0*pw1(k) + pw1(j);
          ycf(idx) := pcf(i)*Third_Semi_Mixed_Derivative
                               (pdg(i).all,cf0,k,j,vrblvl-1);
          ycf(idx) := ycf(idx)*cf1(k)*cf1(j)/create(2.0);
        end loop;
      end loop;
      for j in cf0'range loop -- all fully mixed third derivative terms
        for k in j+1..cf0'last loop
          for L in k+1..cf0'last loop
            idx := idx + 1;
            ydg(idx) := pct(i) + pw1(j) + pw1(k) + pw1(L);
            ycf(idx) := pcf(i)*Third_Fully_Mixed_Derivative
                                 (pdg(i).all,cf0,j,k,L,vrblvl-1);
            ycf(idx) := ycf(idx)*cf1(j)*cf1(k)*cf1(L);
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
    Third_Derivative_First_Order(pcf,pct,pdg,lc0,lc1,lpw,ycf,ydg,vrblvl);
  end Third_Derivative_First_Order;

  procedure Third_Derivative_Second_Order
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
      put_line("third_derivative_second_order 1 ...");
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
      for j in cff'range loop -- all first derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j);
        shared := pcf(i)*Leading_Coefficient(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc1(j);
        idx := idx + 1;
        ydg(idx) := pct(i) + pw2(j);
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
      end loop;
      for j in cff'range loop -- all mixed second derivative terms
        for k in j+1..cff'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + pw1(j) + pw1(k);
          shared := pcf(i)*Second_Mixed_Derivative
                             (pdg(i).all,lc0,j,k,vrblvl-1);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc1(j)*lc1(k); -- /create(2.0);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw2(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc2(j)*lc2(k); -- /create(2.0);
        end loop;
      end loop;
      for j in cff'range loop -- all pure third derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 3.0*pw1(j);
        shared := pcf(i)*Third_Derivative(pdg(i).all,lc0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc1(j)*lc1(j)/create(6.0); -- Taylor series
        idx := idx + 1;
        ydg(idx) := pct(i) + 3.0*pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*lc2(j)*lc2(j)/create(6.0); -- Taylor series
      end loop;
      for j in cff'range loop -- all semi mixed third derivative terms
        for k in j+1..cff'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + 2.0*pw1(j) + pw1(k);
          shared := pcf(i)*Third_Semi_Mixed_Derivative
                             (pdg(i).all,lc0,j,k,vrblvl-1);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc1(j)*lc1(k)/create(2.0);
          idx := idx + 1;
          ydg(idx) := pct(i) + 2.0*pw2(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc2(j)*lc2(k)/create(2.0);
          idx := idx + 1; -- flip role of j and k
          ydg(idx) := pct(i) + 2.0*pw1(k) + pw1(j);
          shared := pcf(i)*Third_Semi_Mixed_Derivative
                             (pdg(i).all,lc0,k,j,vrblvl-1);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc1(k)*lc1(j)/create(2.0);
          idx := idx + 1; -- flip role of j and k
          ydg(idx) := pct(i) + 2.0*pw2(k) + pw2(j);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc2(k)*lc2(j)/create(2.0);
        end loop;
      end loop;
      for j in cff'range loop -- all fully mixed third derivative terms
        for k in j+1..cff'last loop
          for L in k+1..cff'last loop
            idx := idx + 1;
            ydg(idx) := pct(i) + pw1(j) + pw1(k) + pw1(L);
            shared := pcf(i)*Third_Fully_Mixed_Derivative
                               (pdg(i).all,lc0,j,k,L,vrblvl-1);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*lc1(j)*lc1(k)*lc1(L);
            idx := idx + 1;
            ydg(idx) := pct(i) + pw2(j) + pw2(k) + pw2(L);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*lc2(j)*lc2(k)*lc2(L);
          end loop;
        end loop;
      end loop;
    end loop;
    if vrblvl > 0 then
      put("idx : "); put(idx,1); new_line;
      put("ycf'last : "); put(ycf'last,1); new_line;
    end if;
    Double_Real_Powered_Series.Sort(ydg,ycf); 
    Double_Real_Powered_Series.Normalize(ycf,ydg);
  end Third_Derivative_Second_Order;

-- ON A POLYNOMIAL HOMOTOPY :

  procedure First_Derivative_First_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    nbr,idx,size : integer32;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_derivative_first_order 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,1,1,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        First_Derivative_First_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cf0,cf1,pw1,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the first order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
             put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        pw2(i) := ydg(idx);
        cf2(i) := ycf(idx);
      end;
    end loop;
  end First_Derivative_First_Order;

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
    nbr,idx,size : integer32;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_derivative_first_order 3 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,1,1,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
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
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
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
      size := Size_Evaluation(dim,2,1,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Second_Derivative_First_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cf0,cf1,pw1,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the second order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
            put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        pw2(i) := ydg(idx);
        cf2(i) := ycf(idx);
      end;
    end loop;
  end Second_Derivative_First_Order;

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
      put_line("second_derivative_first_order 3 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,2,1,nbr);
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
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    nbr,size,idx : integer32;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("third_derivative_first_order 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,3,1,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Third_Derivative_First_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cf0,cf1,pw1,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the third order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
            put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        cf2(i) := ycf(idx);
        pw2(i) := ydg(idx);
      end;
    end loop;
  end Third_Derivative_First_Order;

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

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("third_derivative_first_order 3 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,3,1,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Third_Derivative_First_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cff,pwr,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the third order evaluation of the polynomial ");
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
