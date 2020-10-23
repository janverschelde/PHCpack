with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Deca_Double_Numbers;                 use Deca_Double_Numbers;

package body DecaDobl_Hessian_Updaters is

  procedure Speel1 ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec ) is

    m1 : integer32;
    powfac : deca_double; -- multiplier factor of two powers

  begin
    m1 := xps(fac(1)); -- the monomial is c*x**m1, m1 >= 2.
    powfac := Deca_Double_Numbers.Create(m1*(m1-1));
    if m1 = 2 then
      H(idx(1),idx(1)) := H(idx(1),idx(1)) + c*powfac;
    elsif m1 = 3 then
      H(idx(1),idx(1)) := H(idx(1),idx(1)) + c*powfac*x(fac(1));
    else -- m > 3, if m = 4, then x**2 at pwt(fac(1))(1)
      H(idx(1),idx(1)) := H(idx(1),idx(1)) + c*powfac*(pwt(fac(1))(m1-3));
    end if;
  end Speel1;

  procedure Speel2 ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec ) is

    m1,m2 : integer32;
    powfac : deca_double; -- multiplier factor of two powers
    acc : Complex_Number;

  begin
    m1 := xps(fac(1));
    powfac := Deca_Double_Numbers.Create(m1*(m1-1));
    if fac'last = 1 then -- the other variable has no higher power
      if fac(1) = idx(1) -- do not forget to multiply with the other var
       then acc := c*powfac*x(idx(2));
       else acc := c*powfac*x(idx(1));
      end if;
    else -- the other variable appears with a higher power
      m2 := xps(fac(2));
      acc := c*powfac*pwt(fac(2))(m2-1);
    end if;
    if m1 = 2 then
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc;
    elsif m1 = 3 then
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*x(fac(1));
    else
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*(pwt(fac(1))(m1-3));
    end if;
    if fac'last = 1 then
      powfac := Deca_Double_Numbers.Create(m1);
      if m1 = 2 then
        H(idx(1),idx(2)) := H(idx(1),idx(2)) + c*powfac*x(fac(1));
      else
        H(idx(1),idx(2)) := H(idx(1),idx(2)) + c*powfac*(pwt(fac(1))(m1-2));
      end if;
    else
      powfac := Deca_Double_Numbers.Create(m2*(m2-1));
      acc := c*powfac*(pwt(fac(1))(m1-1));
      if m2 = 2 then
        H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc;
      elsif m2 = 3 then
        H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*x(fac(2));
      else
        H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*(pwt(fac(2))(m2-3));
      end if;
      powfac := Deca_Double_Numbers.Create(m1*m2);
      acc := c*powfac;
      if m1 = 2
       then acc := acc*x(fac(1));
       else acc := acc*(pwt(fac(1))(m1-2));
      end if;
      if m2 = 2
       then acc := acc*x(fac(2));
       else acc := acc*(pwt(fac(2))(m2-2));
      end if;
      H(idx(1),idx(2)) := H(idx(1),idx(2)) + acc;
    end if;
  end Speel2;

  procedure Speel3 ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec ) is

    m1,m2 : integer32;
    powfac : deca_double; -- multiplier factor of two powers
    acc : Complex_Number;
    offdiagfac : Complex_Number; -- common off diagonal factor
    ondiagfac : Complex_Number;  -- common on diagonal factor

  begin
    offdiagfac := c; ondiagfac := c;  -- off/on diagonal factors
    for i in fac'range loop
      m1 := xps(fac(i));
      if m1 = 2 then
        offdiagfac := offdiagfac*x(fac(i));
      elsif m1 = 3 then
        offdiagfac := offdiagfac*(pwt(fac(i))(1));
        ondiagfac := ondiagfac*x(fac(i));
      else
        offdiagfac := offdiagfac*(pwt(fac(i))(m1-2));
        ondiagfac := ondiagfac*(pwt(fac(i))(m1-3));
      end if;
    end loop;
   -- compute the off diagonal elements of the Hessian
    powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(2)));
    H(idx(1),idx(2)) := H(idx(1),idx(2)) + offdiagfac*powfac*x(idx(3));
    powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(3)));
    H(idx(1),idx(3)) := H(idx(1),idx(3)) + offdiagfac*powfac*x(idx(2));
    powfac := Deca_Double_Numbers.Create(xps(idx(2))*xps(idx(3)));
    H(idx(2),idx(3)) := H(idx(2),idx(3)) + offdiagfac*powfac*x(idx(1));
   -- ten cases for the on diagonal element of the Hessian
    if fac'last = 3 then -- all variables raised to higher power
      m1 := xps(fac(1)); powfac := Deca_Double_Numbers.Create(m1*(m1-1));
      acc := ondiagfac*(pwt(fac(2))(1))*(pwt(fac(3))(1));
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
      m1 := xps(fac(2)); powfac := Deca_Double_Numbers.Create(m1*(m1-1));
      acc := ondiagfac*(pwt(fac(1))(1))*(pwt(fac(3))(1));
      H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
      m1 := xps(fac(3)); powfac := Deca_Double_Numbers.Create(m1*(m1-1));
      acc := ondiagfac*(pwt(fac(1))(1))*(pwt(fac(2))(1));
      H(fac(3),fac(3)) := H(fac(3),fac(3)) + acc*powfac;
    elsif fac'last = 1 then -- one variable raised to higher power
      m1 := xps(fac(1)); powfac := Deca_Double_Numbers.Create(m1*(m1-1));
      if fac(1) = idx(1) then
        acc := ondiagfac*x(idx(2))*x(idx(3));
      elsif fac(1) = idx(2) then
        acc := ondiagfac*x(idx(1))*x(idx(3));
      else -- fac(1) = idx(3)
        acc := ondiagfac*x(idx(1))*x(idx(2));
      end if;
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
    else -- fac'last = 2, two variables raised to higher power
      m1 := xps(fac(1)); powfac := Deca_Double_Numbers.Create(m1*(m1-1));
      m2 := xps(fac(2));
      if fac(1) = idx(1) then
        if fac(2) = idx(2) then -- idx(3) has power 1
          acc := ondiagfac*(pwt(fac(2))(1))*x(idx(3));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := Deca_Double_Numbers.Create(m2*(m2-1));
          acc := ondiagfac*(pwt(fac(1))(1))*x(idx(3));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        else -- idx(2) has power 1
          acc := ondiagfac*x(idx(2))*(pwt(fac(2))(1));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := Deca_Double_Numbers.Create(m2*(m2-1));
          acc := ondiagfac*x(idx(2))*(pwt(fac(1))(1));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        end if;
      elsif fac(1) = idx(2) then
        if fac(2) = idx(1) then -- idx(3) has power 1
          acc := ondiagfac*(pwt(fac(2))(1))*x(idx(3));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := Deca_Double_Numbers.Create(m2*(m2-1));
          acc := ondiagfac*(pwt(fac(1))(1))*x(idx(3));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        else -- idx(1) has power 1
          acc := ondiagfac*x(idx(1))*(pwt(fac(2))(1));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := Deca_Double_Numbers.Create(m2*(m2-1));
          acc := ondiagfac*x(idx(1))*(pwt(fac(1))(1));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        end if;
      else --  fac(1) = idx(3)
        if fac(2) = idx(1) then -- idx(2) has power 1
          acc := ondiagfac*(pwt(fac(2))(1))*x(idx(2));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := Deca_Double_Numbers.Create(m2*(m2-1));
          acc := ondiagfac*(pwt(fac(1))(1))*x(idx(2));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        else -- idx(1) has power 1
          acc := ondiagfac*x(idx(1))*(pwt(fac(2))(1));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := Deca_Double_Numbers.Create(m2*(m2-1));
          acc := ondiagfac*x(idx(1))*(pwt(fac(1))(1));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        end if;
      end if;
    end if;
  end Speel3;

  procedure Speel4 ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec ) is

    m1 : integer32;
    powfac : deca_double; -- multiplier factor of two powers
    acc : Complex_Number;
    offdiagfac : Complex_Number; -- common off diagonal factor
    ondiagfac : Complex_Number;  -- common on diagonal factor
    fwd1 : Complex_Number; -- first forward product
    bck1 : Complex_Number; -- first backward prodcut

  begin
    offdiagfac := c; ondiagfac := c;  -- off/on diagonal factors
    for i in fac'range loop
      m1 := xps(fac(i));
      if m1 = 2 then
        offdiagfac := offdiagfac*x(fac(i));
      elsif m1 = 3 then
        offdiagfac := offdiagfac*(pwt(fac(i))(1));
        ondiagfac := ondiagfac*x(fac(i));
      else
        offdiagfac := offdiagfac*(pwt(fac(i))(m1-2));
        ondiagfac := ondiagfac*(pwt(fac(i))(m1-3));
      end if;
    end loop;
   -- the off diagonal elements use forward and backward products
    fwd1 := x(idx(1))*x(idx(2));
    bck1 := x(idx(4))*x(idx(3));
   -- the last element is a copy of fwd1, with a multiplier factor
    powfac := Deca_Double_Numbers.Create(xps(idx(3))*xps(idx(4)));
    H(idx(3),idx(4)) := H(idx(3),idx(4)) + offdiagfac*powfac*fwd1;
   -- the first element is a copy of bck(1), with a multiplier factor
    powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(2)));
    H(idx(1),idx(2)) := H(idx(1),idx(2)) + offdiagfac*powfac*bck1;
   -- the other off diagonal elements
    acc := offdiagfac*x(idx(2));
    powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(3)));
    H(idx(1),idx(3)) := H(idx(1),idx(3)) + acc*powfac*x(idx(4));
    powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(4)));
    H(idx(1),idx(4)) := H(idx(1),idx(4)) + acc*powfac*x(idx(3));
    acc := offdiagfac*x(idx(1));
    powfac := Deca_Double_Numbers.Create(xps(idx(2))*xps(idx(3)));
    H(idx(2),idx(3)) := H(idx(2),idx(3)) + acc*powfac*x(idx(4));
    powfac := Deca_Double_Numbers.Create(xps(idx(2))*xps(idx(4)));
    H(idx(2),idx(4)) := H(idx(2),idx(4)) + acc*powfac*x(idx(3));
   -- compute the diagonal elements
    for k in fac'range loop
      m1 := xps(fac(k)); powfac := Deca_Double_Numbers.Create(m1*(m1-1));
      acc := powfac*ondiagfac; -- acc is the cofactor
      for i in idx'range loop
        if idx(i) /= fac(k) then -- skip the current factor
          if xps(idx(i)) = 1 
           then acc := acc*x(idx(i));
           else acc := acc*(pwt(idx(i))(1));
          end if;
        end if;
      end loop;
      H(fac(k),fac(k)) := H(fac(k),fac(k)) + acc;
    end loop;
  end Speel4;

  procedure SpeelN ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     fwd : in DecaDobl_Complex_Vectors.Link_to_Vector;
                     bck : in DecaDobl_Complex_Vectors.Link_to_Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec ) is

    sz : constant integer32 := idx'last;
    m1 : integer32;
    powfac : deca_double; -- multiplier factor of two powers
    acc : Complex_Number;
    offdiagfac : Complex_Number; -- common off diagonal factor
    ondiagfac : Complex_Number;  -- common on diagonal factor

  begin
    offdiagfac := c; ondiagfac := c;  -- off/on diagonal factors
    for i in fac'range loop
      m1 := xps(fac(i));
      if m1 = 2 then
        offdiagfac := offdiagfac*x(fac(i));
      elsif m1 = 3 then
        offdiagfac := offdiagfac*(pwt(fac(i))(1));
        ondiagfac := ondiagfac*x(fac(i));
      else
        offdiagfac := offdiagfac*(pwt(fac(i))(m1-2));
        ondiagfac := ondiagfac*(pwt(fac(i))(m1-3));
      end if;
    end loop;
   -- the off diagonal elements use forward and backward products
   -- last element is copy of fwd(sz-3), multiplied with c
    powfac := Deca_Double_Numbers.Create(xps(idx(sz-1))*xps(idx(sz)));
    H(idx(sz-1),idx(sz)) := H(idx(sz-1),idx(sz)) + offdiagfac*powfac*fwd(sz-3);
   -- first element is copy of bck(sz-3), multiplied with c
    powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(2)));
    H(idx(1),idx(2)) := H(idx(1),idx(2)) + offdiagfac*powfac*bck(sz-3);
   -- first row is special, starts with x(idx(2)) after diagonal
    acc := offdiagfac*x(idx(2));
    powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(3)));
    H(idx(1),idx(3)) := H(idx(1),idx(3)) + acc*powfac*bck(sz-4);
    for k in 4..sz-2 loop
      acc := acc*x(idx(k-1));
      powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(k)));
      H(idx(1),idx(k)) := H(idx(1),idx(k)) + acc*powfac*bck(sz-k-1);
    end loop;
    acc := acc*x(idx(sz-2));
    powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(sz-1)));
    H(idx(1),idx(sz-1)) := H(idx(1),idx(sz-1)) + acc*powfac*x(idx(sz));
    powfac := Deca_Double_Numbers.Create(xps(idx(1))*xps(idx(sz)));
    H(idx(1),idx(sz)) := H(idx(1),idx(sz)) + acc*powfac*x(idx(sz-1));
   -- second row is special, starts with x(idx(1)) after diagonal
    acc := offdiagfac*x(idx(1));
    powfac := Deca_Double_Numbers.Create(xps(idx(2))*xps(idx(3)));
    H(idx(2),idx(3)) := H(idx(2),idx(3)) + acc*powfac*bck(sz-4);
    for k in 4..sz-2 loop
      acc := acc*x(idx(k-1));
      powfac := Deca_Double_Numbers.Create(xps(idx(2))*xps(idx(k)));
      H(idx(2),idx(k)) := H(idx(2),idx(k)) + acc*powfac*bck(sz-k-1);
    end loop;
    acc := acc*x(idx(sz-2));
    powfac := Deca_Double_Numbers.Create(xps(idx(2))*xps(idx(sz-1)));
    H(idx(2),idx(sz-1)) := H(idx(2),idx(sz-1)) + acc*powfac*x(idx(sz));
    powfac := Deca_Double_Numbers.Create(xps(idx(2))*xps(idx(sz)));
    H(idx(2),idx(sz)) := H(idx(2),idx(sz)) + acc*powfac*x(idx(sz-1));
   -- the row with index sz-2 has a general formula
    acc := offdiagfac*fwd(sz-4);
    powfac := Deca_Double_Numbers.Create(xps(idx(sz-2))*xps(idx(sz-1)));
    H(idx(sz-2),idx(sz-1)) := H(idx(sz-2),idx(sz-1)) + acc*powfac*x(idx(sz));
    powfac := Deca_Double_Numbers.Create(xps(idx(sz-2))*xps(idx(sz)));
    H(idx(sz-2),idx(sz)) := H(idx(sz-2),idx(sz)) + acc*powfac*x(idx(sz-1));
    for rw in 3..sz-3 loop  -- row rw starts with fwd(rw-2)
      acc := offdiagfac*fwd(rw-2);
      powfac := Deca_Double_Numbers.Create(xps(idx(rw))*xps(idx(rw+1)));
      H(idx(rw),idx(rw+1)) := H(idx(rw),idx(rw+1)) + acc*powfac*bck(sz-rw-2);
      for k in rw+2..sz-2 loop
        acc := acc*x(idx(k-1));
        powfac := Deca_Double_Numbers.Create(xps(idx(rw))*xps(idx(k)));
        H(idx(rw),idx(k)) := H(idx(rw),idx(k)) + acc*powfac*bck(sz-k-1);
      end loop;
      acc := acc*x(idx(sz-2));
      powfac := Deca_Double_Numbers.Create(xps(idx(rw))*xps(idx(sz-1)));
      H(idx(rw),idx(sz-1)) := H(idx(rw),idx(sz-1)) + acc*powfac*x(idx(sz));
      powfac := Deca_Double_Numbers.Create(xps(idx(rw))*xps(idx(sz)));
      H(idx(rw),idx(sz)) := H(idx(rw),idx(sz)) + acc*powfac*x(idx(sz-1));
    end loop;
   -- compute the diagonal elements
    for k in fac'range loop
      m1 := xps(fac(k)); powfac := Deca_Double_Numbers.Create(m1*(m1-1));
      acc := powfac*ondiagfac; -- acc is the cofactor
      for i in idx'range loop
        if idx(i) /= fac(k) then -- skip the current factor
          if xps(idx(i)) = 1 
           then acc := acc*x(idx(i));
           else acc := acc*(pwt(idx(i))(1));
          end if;
        end if;
      end loop;
      H(fac(k),fac(k)) := H(fac(k),fac(k)) + acc;
    end loop;
   -- the above loop for the diagonal elements applies a loop
   -- for the cofactor, a similar triple loop with forward, backward,
   -- and cross porducts is possible for all fac'last cofactors
  end SpeelN;

end DecaDobl_Hessian_Updaters;
