with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body VarbPrec_Polynomial_Evaluations is

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Polynomials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float ) is

    use Standard_Complex_Numbers,Standard_Complex_Polynomials;

    value : Complex_Number := create(0.0);
    absum : double_float := 0.0;

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        for j in 1..t.dg(i) loop
          val := val*z(i);
        end loop;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Standard_Complex_Polynomials.Poly;
               z : Standard_Complex_Vectors.Vector ) return double_float is

    res,absfz,denrco : double_float;

  begin
    Inverse_Condition_Number(f,z,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Polynomials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double ) is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Polynomials;

    zero : constant double_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : double_double := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        for j in 1..t.dg(i) loop
          val := val*z(i);
        end loop;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Polynomials.Poly;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double is

    res,numrco,absfz : double_double;

  begin
    Inverse_Condition_Number(f,z,numrco,absfz,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Polynomials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double ) is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Polynomials;

    zero : constant quad_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : quad_double := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        for j in 1..t.dg(i) loop
          val := val*z(i);
        end loop;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Polynomials.Poly;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double is

    res,absfz,denrco : quad_double;

  begin
    Inverse_Condition_Number(f,z,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Polynomials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number ) is

    use Multprec_Complex_Numbers,Multprec_Complex_Polynomials;

    res : Floating_Number;
    zero : Floating_Number := create(0.0);
    value : Complex_Number := create(zero);
    absum : Floating_Number := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number;
      avl : Floating_Number;

    begin
      Copy(t.cf,val);
      for i in t.dg'range loop
        for j in 1..t.dg(i) loop
          Mul(val,z(i));
        end loop;
      end loop;
      Add(value,val);
      avl := AbsVal(val);
      Add(absum,avl);
      Clear(val); Clear(avl);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    absfz := AbsVal(value); Clear(value);
    denrco := absum;
    rco := absfz/absum; Clear(zero);
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Multprec_Complex_Polynomials.Poly;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number is

    res,absfz,denrco : Floating_Number;

  begin
    Inverse_Condition_Number(f,z,absfz,denrco,res);
    Clear(denrco); Clear(absfz);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Laurentials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float ) is

    use Standard_Complex_Numbers,Standard_Complex_Laurentials;

    value : Complex_Number := create(0.0);
    absum : double_float := 0.0;

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          for j in 1..t.dg(i) loop
            val := val*z(i);
          end loop;
        elsif t.dg(i) < 0 then
          for j in 1..(-t.dg(i)) loop
            val := val/z(i);
          end loop;
        end if;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Standard_Complex_Laurentials.Poly;
               z : Standard_Complex_Vectors.Vector ) return double_float is

    res,absfz,denrco : double_float;

  begin
    Inverse_Condition_Number(f,z,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Laurentials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double ) is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Laurentials;

    zero : constant double_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : double_double := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          for j in 1..t.dg(i) loop
            val := val*z(i);
          end loop;
        elsif t.dg(i) < 0 then
          for j in 1..(-t.dg(i)) loop
            val := val/z(i);
          end loop;
        end if;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Laurentials.Poly;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double is

    res,absfz,denrco : double_double;

  begin
    Inverse_Condition_Number(f,z,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Laurentials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double ) is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Laurentials;

    zero : constant quad_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : quad_double := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          for j in 1..t.dg(i) loop
            val := val*z(i);
          end loop;
        elsif t.dg(i) < 0 then
          for j in 1..(-t.dg(i)) loop
            val := val/z(i);
          end loop;
        end if;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Laurentials.Poly;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double is

    res,absfz,denrco : quad_double;

  begin
    Inverse_Condition_Number(f,z,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Laurentials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number ) is

    use Multprec_Complex_Numbers,Multprec_Complex_Laurentials;

    res : Floating_Number;
    zero : Floating_Number := create(0.0);
    value : Complex_Number := create(zero);
    absum : Floating_Number := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number;
      tmp : Floating_Number;

    begin
      Copy(t.cf,val);
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          for j in 1..t.dg(i) loop
            Mul(val,z(i));
          end loop;
        elsif t.dg(i) < 0 then
          for j in 1..(-t.dg(i)) loop
            Div(val,z(i));
          end loop;
        end if;
      end loop;
      Add(value,val);
      tmp := AbsVal(val);
      Add(absum,tmp);
      Clear(val); Clear(tmp);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    absfz := AbsVal(value); Clear(value);
    denrco := absum;
    rco := absfz/denrco; Clear(zero);
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Multprec_Complex_Laurentials.Poly;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number is

    res,absfz,denrco : Floating_Number;

  begin
    Inverse_Condition_Number(f,z,absfz,denrco,res);
    Clear(denrco); Clear(absfz);
    return res;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Standard_Complex_Poly_Systems.Poly_Sys;
               z : Standard_Complex_Vectors.Vector ) return double_float is

    res : double_float := Inverse_Condition_Number(f(f'first),z);
    rco : double_float;

  begin
    for i in f'first+1..f'last loop
      exit when res + 1.0 = 1.0;
      rco := Inverse_Condition_Number(f(i),z);
      if rco < res
       then res := rco;
      end if;
    end loop;
    return res;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double is

    res : double_double := Inverse_Condition_Number(f(f'first),z);
    one : constant double_double := create(1.0);
    rco : double_double;

  begin
    for i in f'first+1..f'last loop
      exit when res + one = one;
      rco := Inverse_Condition_Number(f(i),z);
      if rco < res
       then res := rco;
      end if;
    end loop;
    return res;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double is

    res : quad_double := Inverse_Condition_Number(f(f'first),z);
    one : constant quad_double := create(1.0);
    rco : quad_double;

  begin
    for i in f'first+1..f'last loop
      exit when res + one = one;
      rco := Inverse_Condition_Number(f(i),z);
      if rco < res
       then res := rco;
      end if;
    end loop;
    return res;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Multprec_Complex_Poly_Systems.Poly_Sys;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number is

    res : Floating_Number := Inverse_Condition_Number(f(f'first),z);
    rco : Floating_Number;

  begin
    for i in f'first+1..f'last loop
      rco := Inverse_Condition_Number(f(i),z);
      if rco < res
       then Copy(rco,res);
      end if;
      Clear(rco);
    end loop;
    return res;
  end Inverse_Condition_Number;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Polynomials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               rco : out double_float;
               fz : out Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Numbers,Standard_Complex_Polynomials;

    value : Complex_Number := create(0.0);
    absum : double_float := 0.0;

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        for j in 1..t.dg(i) loop
          val := val*z(i);
        end loop;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    rco := AbsVal(value)/absum;
    fz := value;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Polynomials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               rco : out double_double;
               fz : out DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Polynomials;

    zero : constant double_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : double_double := zero;

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        for j in 1..t.dg(i) loop
          val := val*z(i);
        end loop;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    rco := AbsVal(value)/absum;
    fz := value;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Polynomials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               rco : out quad_double;
               fz : out QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Polynomials;

    zero : constant quad_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : quad_double := zero;

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        for j in 1..t.dg(i) loop
          val := val*z(i);
        end loop;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    rco := AbsVal(value)/absum;
    fz := value;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Polynomials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               rco : out Floating_Number;
               fz : out Multprec_Complex_Numbers.Complex_Number ) is

    use Multprec_Complex_Numbers,Multprec_Complex_Polynomials;

    res : Floating_Number;
    zero : Floating_Number := create(0.0);
    value : Complex_Number := create(zero);
    absum : Floating_Number := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number;
      avl : Floating_Number;

    begin
      Copy(t.cf,val);
      for i in t.dg'range loop
        for j in 1..t.dg(i) loop
          Mul(val,z(i));
        end loop;
      end loop;
      Add(value,val);
      avl := AbsVal(val);
      Add(absum,avl);
      Clear(val); Clear(avl);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    res := AbsVal(value);
    Div(res,absum);
    Clear(absum); Clear(zero);
    rco := res;
    fz := value;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Laurentials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               rco : out double_float;
               fz : out Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Numbers,Standard_Complex_Laurentials;

    value : Complex_Number := create(0.0);
    absum : double_float := 0.0;

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          for j in 1..t.dg(i) loop
            val := val*z(i);
          end loop;
        elsif t.dg(i) < 0 then
          for j in 1..(-t.dg(i)) loop
            val := val/z(i);
          end loop;
        end if;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    rco := AbsVal(value)/absum;
    fz := value;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Laurentials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               rco : out double_double;
               fz : out DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Laurentials;

    zero : constant double_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : double_double := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          for j in 1..t.dg(i) loop
            val := val*z(i);
          end loop;
        elsif t.dg(i) < 0 then
          for j in 1..(-t.dg(i)) loop
            val := val/z(i);
          end loop;
        end if;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    rco := AbsVal(value)/absum;
    fz := value;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Laurentials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               rco : out quad_double;
               fz : out QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Laurentials;

    zero : constant quad_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : quad_double := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number := t.cf;

    begin
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          for j in 1..t.dg(i) loop
            val := val*z(i);
          end loop;
        elsif t.dg(i) < 0 then
          for j in 1..(-t.dg(i)) loop
            val := val/z(i);
          end loop;
        end if;
      end loop;
      value := value + val;
      absum := absum + AbsVal(val);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    rco := AbsVal(value)/absum;
    fz := value;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Laurentials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               rco : out Floating_Number;
               fz : out Multprec_Complex_Numbers.Complex_Number ) is

    use Multprec_Complex_Numbers,Multprec_Complex_Laurentials;

    res : Floating_Number;
    zero : Floating_Number := create(0.0);
    value : Complex_Number := create(zero);
    absum : Floating_Number := create(0.0);

    procedure Evaluate_Term ( t : in Term; continue : out boolean ) is

      val : Complex_Number;
      tmp : Floating_Number;

    begin
      Copy(t.cf,val);
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          for j in 1..t.dg(i) loop
            Mul(val,z(i));
          end loop;
        elsif t.dg(i) < 0 then
          for j in 1..(-t.dg(i)) loop
            Div(val,z(i));
          end loop;
        end if;
      end loop;
      Add(value,val);
      tmp := AbsVal(val);
      Add(absum,tmp);
      Clear(val); Clear(tmp);
      continue := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Visiting_Iterator(Evaluate_Term);

  begin
    Evaluate_Terms(f);
    res := AbsVal(value);
    Div(res,absum);
    Clear(absum); Clear(zero);
    rco := res;
    fz := value;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
               z : in Standard_Complex_Vectors.Vector;
               rco : out double_float;
               fz : out Standard_Complex_Vectors.Vector ) is

    val : double_float;

  begin
    Evaluate_with_Inverse_Condition(f(f'first),z,rco,fz(fz'first));
    for i in f'first+1..f'last loop
      Evaluate_with_Inverse_Condition(f(i),z,val,fz(i));
      if val < rco
       then rco := val;
      end if;
    end loop;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               rco : out double_double;
               fz : out DoblDobl_Complex_Vectors.Vector ) is

    val : double_double;

  begin
    Evaluate_with_Inverse_Condition(f(f'first),z,rco,fz(fz'first));
    for i in f'first+1..f'last loop
      Evaluate_with_Inverse_Condition(f(i),z,val,fz(i));
      if val < rco
       then rco := val;
      end if;
    end loop;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               rco : out quad_double;
               fz : out QuadDobl_Complex_Vectors.Vector ) is

    val : quad_double;

  begin
    Evaluate_with_Inverse_Condition(f(f'first),z,rco,fz(fz'first));
    for i in f'first+1..f'last loop
      Evaluate_with_Inverse_Condition(f(i),z,val,fz(i));
      if val < rco
       then rco := val;
      end if;
    end loop;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               rco : out Floating_Number;
               fz : out Multprec_Complex_Vectors.Vector ) is

    val : Floating_Number;

  begin
    Evaluate_with_Inverse_Condition(f(f'first),z,rco,fz(fz'first));
    for i in f'first+1..f'last loop
      Evaluate_with_Inverse_Condition(f(i),z,val,fz(i));
      if val < rco
       then Copy(val,rco);
      end if;
    end loop;
    Clear(val);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Laur_Systems.Laur_Sys;
               z : in Standard_Complex_Vectors.Vector;
               rco : out double_float;
               fz : out Standard_Complex_Vectors.Vector ) is

    val : double_float;

  begin
    Evaluate_with_Inverse_Condition(f(f'first),z,rco,fz(fz'first));
    for i in f'first+1..f'last loop
      Evaluate_with_Inverse_Condition(f(i),z,val,fz(i));
      if val < rco
       then rco := val;
      end if;
    end loop;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               rco : out double_double;
               fz : out DoblDobl_Complex_Vectors.Vector ) is

    val : double_double;

  begin
    Evaluate_with_Inverse_Condition(f(f'first),z,rco,fz(fz'first));
    for i in f'first+1..f'last loop
      Evaluate_with_Inverse_Condition(f(i),z,val,fz(i));
      if val < rco
       then rco := val;
      end if;
    end loop;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               rco : out quad_double;
               fz : out QuadDobl_Complex_Vectors.Vector ) is

    val : quad_double;

  begin
    Evaluate_with_Inverse_Condition(f(f'first),z,rco,fz(fz'first));
    for i in f'first+1..f'last loop
      Evaluate_with_Inverse_Condition(f(i),z,val,fz(i));
      if val < rco
       then rco := val;
      end if;
    end loop;
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Laur_Systems.Laur_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               rco : out Floating_Number;
               fz : out Multprec_Complex_Vectors.Vector ) is

    val : Floating_Number;

  begin
    Evaluate_with_Inverse_Condition(f(f'first),z,rco,fz(fz'first));
    for i in f'first+1..f'last loop
      Evaluate_with_Inverse_Condition(f(i),z,val,fz(i));
      if val < rco
       then Copy(val,rco);
      end if;
    end loop;
    Clear(val);
  end Evaluate_with_Inverse_Condition;

end VarbPrec_Polynomial_Evaluations;
