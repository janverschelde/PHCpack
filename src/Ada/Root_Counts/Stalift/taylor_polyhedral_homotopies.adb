with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;

package body Taylor_Polyhedral_Homotopies is

  function Offset ( mic : Mixed_Cell; idx : integer32 )
                  return double_float is

    res : double_float := 0.0;
    lpt : constant Standard_Floating_Vectors.Link_to_Vector
        := Head_Of(mic.pts(idx));

  begin
    for k in lpt'range loop
      res := res + lpt(k)*mic.nor(k);
    end loop;
    return res;
  end Offset;

  procedure Make_Powers
              ( file : in file_type;
                cfq : in Standard_Complex_VecVecs.VecVec;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; 
                pow : out Standard_Floating_VecVecs.VecVec ) is

    tmp : Lists_of_Floating_Vectors.List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;
    pwr : double_float;
    idx : integer32 := 0; -- index of polynomial
    monidx : integer32;   -- index of a monomial coefficient
    cff : Standard_Complex_Numbers.Complex_Number;
    tol_zero : constant double_float := 1.0e-12;

  begin
    for k in lif'range loop
      for i in 1..mix(k) loop
        tmp := lif(k);
        idx := idx + 1;
        declare
          len : constant integer32 := integer32(Length_Of(lif(k)));
        begin
          pow(idx) := new Standard_Floating_Vectors.Vector(1..len);
          monidx := 0;
          while not Is_Null(tmp) loop
            monidx := monidx + 1;
            cff := cfq(idx)(monidx);
            put(file,cff);
            lpt := Head_Of(tmp);
            pwr := lpt(lpt'last) - Offset(mic,k);
            for i in cfq'range loop
              put(file," "); put(file,natural32(lpt(i)),1);
              pwr := pwr + lpt(i)*mic.nor(i);
            end loop;
            put(file," : "); put(file,pwr); new_line(file);
            if pwr < tol_zero
             then pwr := 0.0;
            end if;
            pow(idx)(monidx) := pwr;
            tmp := Tail_Of(tmp);
          end loop;
        end;
      end loop;
    end loop;
  end Make_Powers;

  function Smallest_Nonzero_Power
             ( pow : Standard_Floating_VecVecs.VecVec )
             return double_float is

    res : double_float := -1.0;
    pwr : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pow'range loop
      pwr := pow(i);
      for k in pwr'range loop
        if pwr(k) > 0.0 then
          if res = -1.0 then
            res := pwr(k);
          elsif pwr(k) < res then
            res := pwr(k);
          end if;
        end if;
      end loop;
    end loop;
    return res;
  end Smallest_Nonzero_Power;

  procedure Scale_Powers
              ( pow : in out Standard_Floating_VecVecs.VecVec;
                val : in double_float ) is

    pwr : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pow'range loop
      pwr := pow(i);
      for k in pwr'range loop
        if pwr(k) > 0.0
         then pwr(k) := pwr(k)/val;
        end if;
      end loop;
    end loop;
  end Scale_Powers;

  procedure Make_Homotopy
              ( file : in file_type;
                deg : in integer32;
                pnt : in double_float;
                cfq : in Standard_Complex_VecVecs.VecVec;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                pow : in Standard_Floating_VecVecs.VecVec;
                thm : out Taylor_Homotopy ) is

    tmp : Lists_of_Floating_Vectors.List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;
    pwr : double_float;
    idx : integer32 := 0; -- index of polynomial
    monidx : integer32;   -- index of a monomial coefficient
    cff : Standard_Complex_Numbers.Complex_Number;
    exp : Standard_Integer_Vectors.Vector(cfq'range);
    tol_zero : constant double_float := 1.0e-14;

  begin
    for k in lif'range loop
      for i in 1..mix(k) loop
        tmp := lif(k);
        idx := idx + 1;
        declare
          len : constant integer32 := integer32(Length_Of(lif(k)));
          tmv : Taylor_Monomial_Vector(1..len);
          cnt : integer32 := 0;
        begin
          monidx := 0;
          while not Is_Null(tmp) loop
            monidx := monidx + 1;
            cff := cfq(idx)(monidx);
            lpt := Head_Of(tmp);
            for i in exp'range loop
              exp(i) := integer32(lpt(i));
            end loop;
            pwr := pow(idx)(monidx);
            declare
              tm : constant Link_to_Taylor_Monomial
                 := Make(deg,pwr,pnt,cff,exp);
            begin
              if pwr = 0.0 then -- initial monomial
                put_line(file,"added initial monomial");
                cnt := cnt + 1;
                tmv(cnt) := tm;
              else
                put(file,"lead coefficient of series :");
                put(file,tm.cff(0));
                if abs(tm.cff(0)) < tol_zero then
                  put_line(file," ignore");
                else
                  put_line(file," added");
                  cnt := cnt + 1;
                  tmv(cnt) := tm;
                end if;
              end if;
            end;
            tmp := Tail_Of(tmp);
          end loop;
          thm(idx) := new Taylor_Monomial_Vector'(tmv(1..cnt));
        end;
      end loop;
    end loop;
  end Make_Homotopy;

  procedure Make_Homotopy
              ( file : in file_type; deg : in integer32;
                pnt : in double_float;
                cfq : in Standard_Complex_VecVecs.VecVec;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; thm : out Taylor_Homotopy ) is

    powers : Standard_Floating_VecVecs.VecVec(1..cfq'last);
    minpwr : double_float := 1.0;

  begin
    Make_Powers(file,cfq,mix,lif,mic,powers);
    put_line(file,"The powers of t in the homotopy :");
    put(file,powers);
    minpwr := Smallest_Nonzero_Power(powers);
    put(file,"Smallest nonzero power :"); put(file,minpwr); new_line(file);
    Scale_Powers(powers,minpwr);
    put_line(file,"The powers of t in the homotopy after scaling :");
    put(file,powers);
    Make_Homotopy(file,deg,pnt,cfq,mix,lif,powers,thm);
    Standard_Floating_VecVecs.Clear(powers);
  end Make_Homotopy;

end Taylor_Polyhedral_Homotopies;
