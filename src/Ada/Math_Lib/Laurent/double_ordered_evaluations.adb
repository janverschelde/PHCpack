with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_IO;       use Standard_Integer_Vectors_IO;
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
      elsif ndf = 4 then
        res := 1 + dim + dim*(dim+1)/2 
                 + dim + dim*(dim-1) + dim*(dim-1)*(dim-2)/6
                 + dim + 3*dim*(dim-1)/2 + dim*(dim-1)*(dim-2)/2
                 + dim*(dim-1)*(dim-2)*(dim-3)/24;
        res := res*nbr;
      elsif ndf = 5 then
        res := 1 + dim + dim*(dim+1)/2 
                 + dim + dim*(dim-1) + dim*(dim-1)*(dim-2)/6
                 + dim + 3*dim*(dim-1)/2 + dim*(dim-1)*(dim-2)/2
                       + dim*(dim-1)*(dim-2)*(dim-3)/24
                 + dim + 2*dim*(dim-1) + dim*(dim-1)*(dim-2)
                       + dim*(dim-1)*(dim-2)*(dim-3)/6
                       + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120;
        res := res*nbr;
      elsif ndf = 6 then
        res := 1 + dim + dim*(dim+1)/2 
                 + dim + dim*(dim-1) + dim*(dim-1)*(dim-2)/6
                 + dim + 3*dim*(dim-1)/2 + dim*(dim-1)*(dim-2)/2
                       + dim*(dim-1)*(dim-2)*(dim-3)/24
                 + dim + 2*dim*(dim-1) + dim*(dim-1)*(dim-2)
                       + dim*(dim-1)*(dim-2)*(dim-3)/6
                       + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120
                 + dim + 5*dim*(dim-1)/2 + 10*dim*(dim-1)*(dim-2)/6
                       + 10*dim*(dim-1)*(dim-2)*(dim-3)/24
                       + 5*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120
                       + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)/720;
        res := res*nbr;
      elsif ndf = 7 then
        res := 1 + dim + dim*(dim+1)/2 
                 + dim + dim*(dim-1) + dim*(dim-1)*(dim-2)/6
                 + dim + 3*dim*(dim-1)/2 + dim*(dim-1)*(dim-2)/2
                       + dim*(dim-1)*(dim-2)*(dim-3)/24
                 + dim + 2*dim*(dim-1) + dim*(dim-1)*(dim-2)
                       + dim*(dim-1)*(dim-2)*(dim-3)/6
                       + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120
                 + dim + 5*dim*(dim-1)/2 + 10*dim*(dim-1)*(dim-2)/6
                       + 10*dim*(dim-1)*(dim-2)*(dim-3)/24
                       + 5*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120
                       + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)/720
                 + dim + 3*dim*(dim-1) + 15*dim*(dim-1)*(dim-2)/6
                       + 20*dim*(dim-1)*(dim-2)*(dim-3)/24
                       + 15*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120
                       + 6*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)/720
                 + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)*(dim-6)/5040;
        res := res*nbr;
      elsif ndf = 8 then
        res := 1 + dim + dim*(dim+1)/2 
                 + dim + dim*(dim-1) + dim*(dim-1)*(dim-2)/6
                 + dim + 3*dim*(dim-1)/2 + dim*(dim-1)*(dim-2)/2
                       + dim*(dim-1)*(dim-2)*(dim-3)/24
                 + dim + 2*dim*(dim-1) + dim*(dim-1)*(dim-2)
                       + dim*(dim-1)*(dim-2)*(dim-3)/6
                       + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120
                 + dim + 5*dim*(dim-1)/2 + 10*dim*(dim-1)*(dim-2)/6
                       + 10*dim*(dim-1)*(dim-2)*(dim-3)/24
                       + 5*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120
                       + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)/720
                 + dim + 3*dim*(dim-1) + 15*dim*(dim-1)*(dim-2)/6
                       + 20*dim*(dim-1)*(dim-2)*(dim-3)/24
                       + 15*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120
                       + 6*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)/720
                 + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)*(dim-6)/5040
          + dim + 7*dim*(dim-1)/2 + 21*dim*(dim-1)*(dim-2)/6
                + 35*dim*(dim-1)*(dim-2)*(dim-3)/24
                + 35*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)/120
                + 21*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)/720
          + 7*dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)*(dim-6)/5040
          + dim*(dim-1)*(dim-2)*(dim-3)*(dim-4)*(dim-5)*(dim-6)*(dim-7)/40320;
        res := res*nbr;
      end if;
    elsif ord = 2 then
      if ndf = 1 then -- first derivative second order
        res := (2*dim+1)*nbr;
      elsif ndf = 2 then
        res := (1 + 2*dim + 3*dim + 2*dim*(dim-1))*nbr;
      elsif ndf = 3 then
        res := 1 + 2*dim + 3*dim + 2*dim*(dim-1)
             + 4*dim + 6*dim*(dim-1) + 8*dim*(dim-1)*(dim-2)/6;
        res := res*nbr;
      end if;
    end if;
    return res;
  end Size_Evaluation;

-- ON ONE POLYNOMIAL :

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
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;
    shared : Complex_Number;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_derivative_second_order 0 ...");
    end if;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,0,vrblvl-1);
      for j in cf0'range loop -- all first order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j);
        shared := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j);
        idx := idx + 1;              -- second order terms
        ydg(idx) := pct(i) + pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf2(j);
      end loop;
    end loop;
    if vrblvl > 0 then
      put("idx : "); put(idx,1); new_line;
      put("ycf'last : "); put(ycf'last,1); new_line;
    end if;
    Double_Real_Powered_Series.Sort(ydg,ycf); 
    Double_Real_Powered_Series.Normalize(ycf,ydg);
  end First_Derivative_Second_Order;

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
    First_Derivative_Second_Order
      (pcf,pct,pdg,lc0,lc1,lc2,pw1,pw2,ycf,ydg,vrblvl);
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
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;
    shared : Complex_Number;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("second_derivative_second_order 0 ...");
    end if;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,0,vrblvl-1);
      for j in cf0'range loop -- all first order terms
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j);
        shared := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j);
        idx := idx + 1;
        ydg(idx) := pct(i) + pw2(j); -- second order term
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf2(j);
      end loop;
      for j in cf0'range loop -- all pure second derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*pw1(j);
        shared := pcf(i)*Second_Derivative(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j)*cf1(j)/create(2.0); -- Taylor series
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf2(j)*cf2(j)/create(2.0); -- Taylor series
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j) + pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j)*cf2(j); -- double compensates /2
      end loop;
      for j in cf0'range loop -- all mixed second derivative terms
        for k in j+1..cf0'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + pw1(j) + pw1(k);
          shared
            := pcf(i)*Second_Mixed_Derivative(pdg(i).all,cf0,j,k,vrblvl-1);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(j)*cf1(k);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw2(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf2(j)*cf2(k);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw1(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(j)*cf2(k);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw2(j) + pw1(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf2(j)*cf1(k);
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
          ycf(idx) := ycf(idx)*lc1(j)*lc1(k);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw2(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc2(j)*lc2(k);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw1(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc1(j)*lc2(k);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw2(j) + pw1(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*lc1(j)*lc1(k);
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
        ycf(idx) := ycf(idx)*cf1(j)*cf1(j)*cf1(j)/create(6.0); -- Taylor series
      end loop;
      for j in cf0'range loop -- all semi mixed third derivative terms
        for k in j+1..cf0'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + 2.0*pw1(j) + pw1(k);
          ycf(idx) := pcf(i)*Third_Semi_Mixed_Derivative
                               (pdg(i).all,cf0,j,k,vrblvl-1);
          ycf(idx) := ycf(idx)*cf1(j)*cf1(j)*cf1(k)/create(2.0);
          idx := idx + 1; -- flip role of j and k
          ydg(idx) := pct(i) + 2.0*pw1(k) + pw1(j);
          ycf(idx) := pcf(i)*Third_Semi_Mixed_Derivative
                               (pdg(i).all,cf0,k,j,vrblvl-1);
          ycf(idx) := ycf(idx)*cf1(k)*cf1(k)*cf1(j)/create(2.0);
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
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;
    shared : Complex_Number;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("third_derivative_second_order 0 ...");
    end if;
    for i in pcf'range loop -- run over all monomials
      idx := idx + 1;
      ydg(idx) := pct(i);     -- first the constant term
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,0,vrblvl-1);
      for j in cf0'range loop -- all first derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j);
        shared := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j);
        idx := idx + 1;
        ydg(idx) := pct(i) + pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf2(j);
      end loop;
      for j in cf0'range loop -- all pure second derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*pw1(j);
        shared := pcf(i)*Second_Derivative(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j)*cf1(j)/create(2.0); -- Taylor series
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf2(j)*cf2(j)/create(2.0); -- Taylor series
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j) + pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j)*cf2(j);
      end loop;
      for j in cf0'range loop -- all mixed second derivative terms
        for k in j+1..cf0'last loop
          idx := idx + 1;
          ydg(idx) := pct(i) + pw1(j) + pw1(k);
          shared := pcf(i)*Second_Mixed_Derivative
                             (pdg(i).all,cf0,j,k,vrblvl-1);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(j)*cf1(k); 
          idx := idx + 1;
          ydg(idx) := pct(i) + pw2(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf2(j)*cf2(k);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw1(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(j)*cf2(k);
          idx := idx + 1;
          ydg(idx) := pct(i) + pw2(j) + pw1(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf2(j)*cf1(k);
        end loop;
      end loop;
      for j in cf0'range loop -- all pure third derivative terms
        idx := idx + 1;
        ydg(idx) := pct(i) + 3.0*pw1(j);
        shared := pcf(i)*Third_Derivative(pdg(i).all,cf0,j,vrblvl-1);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j)*cf1(j)*cf1(j)/create(6.0);
        idx := idx + 1;
        ydg(idx) := pct(i) + 3.0*pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf2(j)*cf2(j)*cf2(j)/create(6.0);
        idx := idx + 1;
        ydg(idx) := pct(i) + 2.0*pw1(j) + pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j)*cf1(j)*cf2(j)/create(2.0);
        idx := idx + 1;
        ydg(idx) := pct(i) + pw1(j) + 2.0*pw2(j);
        ycf(idx) := shared;
        ycf(idx) := ycf(idx)*cf1(j)*cf2(j)*cf2(j)/create(2.0);
      end loop;
      for j in cf0'range loop -- all semi mixed third derivative terms
        for k in j+1..cf0'last loop
         -- (2,1) for components (j,k)
          shared := pcf(i)*Third_Semi_Mixed_Derivative
                             (pdg(i).all,cf0,j,k,vrblvl-1);
          idx := idx + 1; -- 2*alpha_j + alpha_k
          ydg(idx) := pct(i) + 2.0*pw1(j) + pw1(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(j)*cf1(j)*cf1(k)/create(2.0);
          idx := idx + 1; -- 2*alpha_j + beta_k
          ydg(idx) := pct(i) + 2.0*pw1(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(j)*cf1(j)*cf2(k)/create(2.0);
          idx := idx + 1; -- 2*beta_j + alpha_k
          ydg(idx) := pct(i) + 2.0*pw2(j) + pw1(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf2(j)*cf2(j)*cf1(k)/create(2.0);
          idx := idx + 1; -- 2*beta_j + beta_k
          ydg(idx) := pct(i) + 2.0*pw2(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf2(j)*cf2(j)*cf2(k)/create(2.0);
          idx := idx + 1; -- alpha_j + beta_j + alpha_k
          ydg(idx) := pct(i) + pw1(j) + pw2(j) + pw1(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(j)*cf2(j)*cf1(k); -- /create(2.0);
          idx := idx + 1; -- alpha_j + beta_j + beta_k
          ydg(idx) := pct(i) + pw1(j) + pw2(j) + pw2(k);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(j)*cf2(j)*cf2(k); -- /create(2.0);
         -- (2,1) for components (k,j) flipping role of j and k
          shared := pcf(i)*Third_Semi_Mixed_Derivative
                             (pdg(i).all,cf0,k,j,vrblvl-1);
          idx := idx + 1; -- 2*alpha_k + alpha_j
          ydg(idx) := pct(i) + 2.0*pw1(k) + pw1(j);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(k)*cf1(k)*cf1(j)/create(2.0);
          idx := idx + 1; -- 2*alpha_k + beta_j
          ydg(idx) := pct(i) + 2.0*pw1(k) + pw2(j);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(k)*cf1(k)*cf2(j)/create(2.0);
          idx := idx + 1;  -- 2*beta_k + alpha_j
          ydg(idx) := pct(i) + 2.0*pw2(k) + pw1(j);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf2(k)*cf2(k)*cf1(j)/create(2.0);
          idx := idx + 1;  -- 2*beta_k + beta_j
          ydg(idx) := pct(i) + 2.0*pw2(k) + pw2(j);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf2(k)*cf2(k)*cf2(j)/create(2.0);
          idx := idx + 1;  -- alpha_k + beta_k + alpha_j
          ydg(idx) := pct(i) + pw1(k) + pw2(k) + pw1(j);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(k)*cf2(k)*cf1(j); -- /create(2.0);
          idx := idx + 1;  -- alpha_k + beta_k + beta_j
          ydg(idx) := pct(i) + pw1(k) + pw2(k) + pw2(j);
          ycf(idx) := shared;
          ycf(idx) := ycf(idx)*cf1(k)*cf2(k)*cf2(j); -- /create(2.0);
        end loop;
      end loop;
      for j in cf0'range loop -- all fully mixed third derivative terms
        for k in j+1..cf0'last loop
          for L in k+1..cf0'last loop -- 8 = 2^3 combinations
            idx := idx + 1;
            ydg(idx) := pct(i) + pw1(j) + pw1(k) + pw1(L);
            shared := pcf(i)*Third_Fully_Mixed_Derivative
                               (pdg(i).all,cf0,j,k,L,vrblvl-1);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*cf1(j)*cf1(k)*cf1(L);
            idx := idx + 1;
            ydg(idx) := pct(i) + pw1(j) + pw1(k) + pw2(L);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*cf1(j)*cf1(k)*cf2(L);
            idx := idx + 1;
            ydg(idx) := pct(i) + pw1(j) + pw2(k) + pw1(L);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*cf1(j)*cf2(k)*cf1(L);
            idx := idx + 1;
            ydg(idx) := pct(i) + pw1(j) + pw2(k) + pw2(L);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*cf1(j)*cf2(k)*cf2(L);
            idx := idx + 1;
            ydg(idx) := pct(i) + pw2(j) + pw1(k) + pw1(L);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*cf2(j)*cf1(k)*cf1(L);
            idx := idx + 1;
            ydg(idx) := pct(i) + pw2(j) + pw1(k) + pw2(L);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*cf2(j)*cf1(k)*cf2(L);
            idx := idx + 1;
            ydg(idx) := pct(i) + pw2(j) + pw2(k) + pw1(L);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*cf2(j)*cf2(k)*cf1(L);
            idx := idx + 1;
            ydg(idx) := pct(i) + pw2(j) + pw2(k) + pw2(L);
            ycf(idx) := shared;
            ycf(idx) := ycf(idx)*cf2(j)*cf2(k)*cf2(L);
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
    Third_Derivative_Second_Order
      (pcf,pct,pdg,lc0,lc1,lc2,pw1,pw2,ycf,ydg,vrblvl);
  end Third_Derivative_Second_Order;

-- INDEXED DERIVATIVES ON ONE POLYNOMIAL :

  procedure Fixed_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                difsum : in integer32; idxnxt : in out integer32;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := cf0'last;
    lstidx : integer32; -- last index
    cnt : integer32 := 0;
    deg : Standard_Integer_Vectors.Vector(1..dim);
    pdgidx : integer32; -- index of the monomial

    procedure Next ( difidx : in Standard_Integer_Vectors.Vector;
                     continue : out boolean ) is

      val : Complex_Number;

    begin
      cnt := cnt + 1;
      if vrblvl > 0 then
        put(cnt,3); put(" :"); put(difidx); new_line;
      end if;
      deg := pdg(pdgidx).all;
      val := Double_Leading_Evaluations.Indexed_Derivative(deg,cf0,difidx);
      val := pcf(pdgidx)*val;
      for k in difidx'range loop   -- multiply with cf1(k)**difidx(k)
        for j in 1..difidx(k) loop
          val := val*cf1(k);
          val := val/double_float(j); -- divide by difidx(k)!
        end loop;
      end loop;
      ycf(idxnxt) := val;
      ydg(idxnxt) := pct(pdgidx);
      for k in difidx'range loop
        if difidx(k) > 0 
         then ydg(idxnxt) := ydg(idxnxt) + double_float(difidx(k))*pw1(k);
        end if;
      end loop;
      idxnxt := idxnxt + 1;
      continue := true;
    end Next;

    procedure Differentiation_Indices is
      new Double_Leading_Evaluations.Enumerate_Indices(Next);

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("fixed_derivative_first_order 0 ...");
      put("difsum : "); put(difsum,1);
      put(", idxnxt : "); put(idxnxt,1); new_line;
    end if;
    pdgidx := pdg'first-1;
    while pdgidx < pdg'last loop
      pdgidx := pdgidx + 1;
      Differentiation_Indices(dim,difsum);
    end loop;
    if vrblvl > 0
     then put("idxnxt : "); put(idxnxt,1); new_line;
    end if;
    lstidx := idxnxt-1;
    Double_Real_Powered_Series.Sort(ydg(ydg'first..lstidx),
                                    ycf(ycf'first..lstidx)); 
    Double_Real_Powered_Series.Normalize(ycf(ycf'first..lstidx),
                                         ydg(ydg'first..lstidx));
  end Fixed_Derivative_First_Order;

  function Positive_Count
             ( v : Standard_Integer_Vectors.Vector ) return integer32 is

  -- DESCRIPTION :
  --   Returns the number of positive values in v.

    res : integer32 := 0;

  begin
    for i in v'range loop
      if v(i) > 0 
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Positive_Count;

  function Positive_Index
             ( v : Standard_Integer_Vectors.Vector ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index of the first positive value in v.

  begin
    for i in v'range loop
      if v(i) > 0
       then return i;
      end if;
    end loop;
    return 0;
  end Positive_Index;

  procedure Fixed_Derivative_Second_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                difsum : in integer32; idxnxt : in out integer32;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := cf0'last;
    lstidx : integer32; -- last index
    cnt : integer32 := 0;
    deg : Standard_Integer_Vectors.Vector(1..dim);
    pdgidx : integer32; -- index of the monomial

    procedure Accumulate
                ( k : in integer32;
                  difidx : in Standard_Integer_Vectors.Vector;
                  ycfval : in Complex_Number; ydgval : in double_float ) is

    -- DESCRIPTION :
    --   Accumulates the values for ycf and ydg at the current index k
    --   of the differentiation indices.

      newycf : Complex_Number := ycfval;
      newydg : double_float := ydgval;
      difval : integer32;
 
    begin
      if vrblvl > 0 then
        put("accumulating at k = "); put(k,1); put_line(" ...");
      end if;
      if k > difidx'last then
        if vrblvl > 0 then
          put("exponent "); put(idxnxt,1); put(" :"); put(ydgval); new_line;
          put("coefficient : "); put(ycfval); new_line;
        end if;
        ycf(idxnxt) := ycfval;
        ydg(idxnxt) := ydgval;
        idxnxt := idxnxt + 1;
      elsif difidx(k) = 0 then -- skip component k
        Accumulate(k+1,difidx,ycfval,ydgval);
      else -- make the combinations
        difval := difidx(k);
        newycf := ycfval;
        for j in 1..difval loop -- multiply with cf1(k)**difval
          newycf := newycf*cf1(k)/double_float(j);
        end loop;
        newydg := ydgval + double_float(difval)*pw1(k);
        Accumulate(k+1,difidx,newycf,newydg);
        for j in 1..difval-1 loop
          newycf := ycfval;
          for kk in 1..difval-j loop -- multiply with cf1(k)
            newycf := newycf*cf1(k)/double_float(kk);
          end loop;
          for kk in 1..j loop -- multiply with cf2(k)
            newycf := newycf*cf2(k)/double_float(kk);
          end loop;
          newydg := ydgval + double_float(difval-j)*pw1(k)
                           + double_float(j)*pw2(k);
          Accumulate(k+1,difidx,newycf,newydg);
        end loop;
        newycf := ycfval;
        for j in 1..difval loop -- multiply with cf2(k)**difval
          newycf := newycf*cf2(k)/double_float(j);
        end loop;
        newydg := ydgval + double_float(difval)*pw2(k);
        Accumulate(k+1,difidx,newycf,newydg);
      end if;
    end Accumulate;

    procedure Next ( difidx : in Standard_Integer_Vectors.Vector;
                     continue : out boolean ) is

      val : Complex_Number;
      poscnt : constant integer32 := Positive_Count(difidx);
      posidx,posval : integer32;

    begin
      cnt := cnt + 1;
      if vrblvl > 0 then
        put(cnt,3); put(" :"); put(difidx); new_line;
      end if;
      deg := pdg(pdgidx).all;
      val := Double_Leading_Evaluations.Indexed_Derivative(deg,cf0,difidx);
      val := pcf(pdgidx)*val;
      ycf(idxnxt) := val;
      if poscnt = 1 then -- no mixed derivative
        posidx := Positive_Index(difidx);
        posval := difidx(posidx);
        for j in 1..posval loop -- multiply with cf1(posidx)**posval
          ycf(idxnxt) := ycf(idxnxt)*cf1(posidx)/double_float(j);
        end loop;
        ydg(idxnxt) := pct(pdgidx) + double_float(posval)*pw1(posidx);
        if vrblvl > 0 then
          put("exponent "); put(idxnxt,1); put(" :");
          put(ydg(idxnxt)); new_line;
          put("coefficient : "); put(ycf(idxnxt)); new_line;
        end if;
        for j in 1..posval-1 loop
          ycf(idxnxt+j) := val;
          for k in 1..posval-j loop -- multiply with cf1(posidx)
            ycf(idxnxt+j) := ycf(idxnxt+j)*cf1(posidx)/double_float(k);
          end loop;
          for k in 1..j loop -- multiply with cf2(posidx)
            ycf(idxnxt+j) := ycf(idxnxt+j)*cf2(posidx)/double_float(k);
          end loop;
          ydg(idxnxt+j) := pct(pdgidx)
                         + double_float(posval-j)*pw1(posidx)
                         + double_float(j)*pw2(posidx);
          if vrblvl > 0 then
            put("exponent "); put(idxnxt+j,1); put(" :");
            put(ydg(idxnxt+j)); new_line;
            put("coefficient : "); put(ycf(idxnxt+j)); new_line;
          end if;
        end loop;
        idxnxt := idxnxt + posval;
        ycf(idxnxt) := val;
        for j in 1..posval loop -- multiply with cf2(posidx)**posval
          ycf(idxnxt) := ycf(idxnxt)*cf2(posidx)/double_float(j);
        end loop;
        ydg(idxnxt) := pct(pdgidx) + double_float(posval)*pw2(posidx);
        if vrblvl > 0 then
          put("exponent "); put(idxnxt,1); put(" :");
          put(ydg(idxnxt)); new_line;
          put("coefficient : "); put(ycf(idxnxt)); new_line;
        end if;
        idxnxt := idxnxt + 1;
      else
        Accumulate(difidx'first,difidx,val,pct(pdgidx));
      end if;
      continue := true;
    end Next;

    procedure Differentiation_Indices is
      new Double_Leading_Evaluations.Enumerate_Indices(Next);

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("fixed_derivative_second_order 0 ...");
      put("difsum : "); put(difsum,1);
      put(", idxnxt : "); put(idxnxt,1); new_line;
    end if;
    pdgidx := pdg'first-1;
    while pdgidx < pdg'last loop
      pdgidx := pdgidx + 1;
      Differentiation_Indices(dim,difsum);
    end loop;
    if vrblvl > 0
     then put("idxnxt : "); put(idxnxt,1); new_line;
    end if;
    lstidx := idxnxt-1;
    Double_Real_Powered_Series.Sort(ydg(ydg'first..lstidx),
                                    ycf(ycf'first..lstidx)); 
    Double_Real_Powered_Series.Normalize(ycf(ycf'first..lstidx),
                                         ydg(ydg'first..lstidx));
  end Fixed_Derivative_Second_Order;

  procedure First_Order_Evaluation
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                difmax : in integer32;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    idx : integer32 := ycf'first-1;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_order_evaluation 1 ...");
      put("difmax : "); put(difmax,1); new_line;
    end if;
    for i in pcf'range loop -- first compute all zero-th order terms
      idx := idx + 1;
      ydg(idx) := pct(i);
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,0,vrblvl-1);
    end loop;
    Double_Real_Powered_Series.Sort(ydg(ydg'first..idx),
                                    ycf(ycf'first..idx)); 
    Double_Real_Powered_Series.Normalize(ycf(ycf'first..idx),
                                         ydg(ydg'first..idx));
    idx := idx + 1;
    for difsum in 1..difmax loop
      Fixed_Derivative_First_Order
        (pcf,pct,pdg,cf0,cf1,pw1,difsum,idx,ycf,ydg,vrblvl-1);
    end loop;
    if vrblvl > 0 then
      put("ycf'last : "); put(ycf'last,1);
      put(", idx : "); put(idx,1); new_line;
    end if;
  end First_Order_Evaluation;

  procedure Second_Order_Evaluation
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                difmax : in integer32;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    idx : integer32 := ycf'first-1;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("second_order_evaluation 1 ...");
      put("difmax : "); put(difmax,1); new_line;
    end if;
    for i in pcf'range loop -- first compute all zero-th order terms
      idx := idx + 1;
      ydg(idx) := pct(i);
      ycf(idx) := pcf(i)*Leading_Coefficient(pdg(i).all,cf0,0,vrblvl-1);
    end loop;
    Double_Real_Powered_Series.Sort(ydg(ydg'first..idx),
                                    ycf(ycf'first..idx)); 
    Double_Real_Powered_Series.Normalize(ycf(ycf'first..idx),
                                         ydg(ydg'first..idx));
    idx := idx + 1;
    for difsum in 1..difmax loop
      Fixed_Derivative_Second_Order
        (pcf,pct,pdg,cf0,cf1,cf2,pw1,pw2,difsum,idx,ycf,ydg,vrblvl-1);
    end loop;
    if vrblvl > 0 then
      put("ycf'last : "); put(ycf'last,1);
      put(", idx : "); put(idx,1); new_line;
    end if;
  end Second_Order_Evaluation;

-- ON A LAURENT HOMOTOPY :

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
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                cf3 : out Standard_Complex_Vectors.Vector;
                pw3 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    nbr,idx,size : integer32;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_derivative_second_order 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,1,2,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        First_Derivative_Second_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cf0,cf1,cf2,pw1,pw2,
           ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the first derivative second order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
            put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        pw3(i) := ydg(idx);
        cf3(i) := ycf(idx);
      end;
    end loop;
  end First_Derivative_Second_Order;

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
    nbr,idx,size : integer32;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_derivative_second_order 3 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,1,2,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        First_Derivative_Second_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cff,pwr,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the first derivative second order evaluation of polynomial ");
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
          put("the second derivative first order evaluation of polynomial ");
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
          put("the second derivative first order evaluation of polynomial ");
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
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                cf3 : out Standard_Complex_Vectors.Vector;
                pw3 : out Standard_Floating_Vectors.Vector;
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
      size := Size_Evaluation(dim,2,2,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Second_Derivative_Second_Order
          (hcf(i).all,hct(i).all,hdg(i).all,
           cf0,cf1,cf2,pw1,pw2,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the second derivative second order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
            put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        pw3(i) := ydg(idx);
        cf3(i) := ycf(idx);
      end;
    end loop;
  end Second_Derivative_Second_Order;

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
      put_line("second_derivative_second_order 3 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,2,2,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Second_Derivative_Second_Order
          (hcf(i).all,hct(i).all,hdg(i).all,cff,pwr,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the second derivative second order evaluation of polynomial ");
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
          put("third derivative first order evaluation of polynomial ");
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
          put("third derivative first order evaluation of the polynomial ");
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

  procedure Third_Derivative_Second_Order
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                cf3 : out Standard_Complex_Vectors.Vector;
                pw3 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    nbr,size,idx : integer32;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("third_derivative_second_order 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,3,2,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Third_Derivative_Second_Order
          (hcf(i).all,hct(i).all,hdg(i).all,
           cf0,cf1,cf2,pw1,pw2,ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("the third derivative second order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
            put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        pw3(i) := ydg(idx);
        cf3(i) := ycf(idx);
      end;
    end loop;
  end Third_Derivative_Second_Order;

-- INDEXED DERIVATIVES ON A LAURENT HOMOTOPY :

  procedure First_Order_Evaluation
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                difmax : in integer32;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    nbr,size,idx : integer32;

  begin
    if vrblvl > 0 then
      put("-> in Double_Ordered_Evaluations.");
      put_line("first_order_evaluation 2 ...");
    end if;
    for i in hcf'range loop
      nbr := hcf(i)'last;
      size := Size_Evaluation(dim,difmax,1,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        First_Order_Evaluation
          (hcf(i).all,hct(i).all,hdg(i).all,cf0,cf1,pw1,difmax,
           ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("derivative "); put(difmax,1);
          put("first order evaluation of polynomial ");
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
  end First_Order_Evaluation;

  procedure Second_Order_Evaluation
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cf1 : in Standard_Complex_Vectors.Vector;
                cf2 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                pw2 : in Standard_Floating_Vectors.Vector;
                difmax : in integer32;
                cf3 : out Standard_Complex_Vectors.Vector;
                pw3 : out Standard_Floating_Vectors.Vector;
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
      size := Size_Evaluation(dim,difmax,2,nbr);
      declare
        ycf : Standard_Complex_Vectors.Vector(1..size);
        ydg : Standard_Floating_Vectors.Vector(1..size);
      begin
        Second_Order_Evaluation
          (hcf(i).all,hct(i).all,hdg(i).all,cf0,cf1,cf2,pw1,pw2,difmax,
           ycf,ydg,vrblvl-1);
        if vrblvl > 0 then
          put("derivative "); put(difmax,1);
          put("second order evaluation of polynomial ");
          put(i,1); put_line(" :");
          for i in ycf'range loop
            put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
          end loop;
        end if;
        idx := Double_Real_Powered_Series.Positive_Minimum_Index(ycf,ydg);
        pw3(i) := ydg(idx);
        cf3(i) := ycf(idx);
      end;
    end loop;
  end Second_Order_Evaluation;

end Double_Ordered_Evaluations;
