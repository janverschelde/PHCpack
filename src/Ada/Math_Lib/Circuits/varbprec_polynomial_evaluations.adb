with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body VarbPrec_Polynomial_Evaluations is

-- PART I : ordinary polynomials

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Polynomials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               fz : out Standard_Complex_Numbers.Complex_Number;
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
    fz := value;
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Standard_Complex_Polynomials.Poly;
               z : Standard_Complex_Vectors.Vector ) return double_float is

    res,absfz,denrco : double_float;
    fz : Standard_Complex_Numbers.Complex_Number;

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Polynomials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               fz : out DoblDobl_Complex_Numbers.Complex_Number;
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
    fz := value;
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Polynomials.Poly;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double is

    res,numrco,absfz : double_double;
    fz : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    Inverse_Condition_Number(f,z,fz,numrco,absfz,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Polynomials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               fz : out QuadDobl_Complex_Numbers.Complex_Number;
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
    fz := value;
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Polynomials.Poly;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double is

    res,absfz,denrco : quad_double;
    fz : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Polynomials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               fz : out Multprec_Complex_Numbers.Complex_Number;
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
    fz := value;
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/absum; Clear(zero);
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Multprec_Complex_Polynomials.Poly;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number is

    res,absfz,denrco : Floating_Number;
    fz : Multprec_Complex_Numbers.Complex_Number;

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    Multprec_Complex_Numbers.Clear(fz);
    Clear(absfz); Clear(denrco);
    return res;
  end Inverse_Condition_Number;

-- PART II : Laurent polynomials

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Laurentials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               fz : out Standard_Complex_Numbers.Complex_Number;
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
    fz := value;
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Standard_Complex_Laurentials.Poly;
               z : Standard_Complex_Vectors.Vector ) return double_float is

    res,absfz,denrco : double_float;
    fz : Standard_Complex_Numbers.Complex_Number;

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Laurentials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               fz : out DoblDobl_Complex_Numbers.Complex_Number;
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
    fz := value;
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Laurentials.Poly;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double is

    res,absfz,denrco : double_double;
    fz : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Laurentials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               fz : out QuadDobl_Complex_Numbers.Complex_Number;
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
    fz := value;
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Laurentials.Poly;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double is

    res,absfz,denrco : quad_double;
    fz : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Laurentials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               fz : out Multprec_Complex_Numbers.Complex_Number;
               absfz,denrco,rco : out Floating_Number ) is

    use Multprec_Complex_Numbers,Multprec_Complex_Laurentials;

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
    fz := value;
    absfz := AbsVal(value);
    denrco := absum;
    rco := absfz/denrco; Clear(zero);
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Multprec_Complex_Laurentials.Poly;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number is

    res,absfz,denrco : Floating_Number;
    fz : Multprec_Complex_Numbers.Complex_Number;

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    Multprec_Complex_Numbers.Clear(fz);
    Clear(denrco); Clear(absfz);
    return res;
  end Inverse_Condition_Number;

-- PART III : ordinary polynomial systems

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
               z : in Standard_Complex_Vectors.Vector;
               fz : out Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float ) is

    wrk,wrk_absfz,wrk_denrco : double_float;

  begin
    Inverse_Condition_Number(f(f'first),z,fz(fz'first),absfz,denrco,rco);
    for i in f'first+1..f'last loop
      exit when rco + 1.0 = 1.0;
      Inverse_Condition_Number(f(i),z,fz(i),wrk_absfz,wrk_denrco,wrk);
      if wrk < rco
       then rco := wrk; absfz := wrk_absfz; denrco := wrk_denrco;
      end if;
    end loop;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Standard_Complex_Poly_Systems.Poly_Sys;
               z : Standard_Complex_Vectors.Vector ) return double_float is

    res,absfz,denrco : double_float;
    fz : Standard_Complex_Vectors.Vector(f'range);

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               fz : out DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double ) is

    one : constant double_double := create(1.0);
    wrk,wrk_absfz,wrk_denrco : double_double;

  begin
    Inverse_Condition_Number(f(f'first),z,fz(fz'first),absfz,denrco,rco);
    for i in f'first+1..f'last loop
      exit when rco + one = one;
      Inverse_Condition_Number(f(i),z,fz(i),wrk_absfz,wrk_denrco,wrk);
      if wrk < rco
       then rco := wrk; absfz := wrk_absfz; denrco := wrk_denrco;
      end if;
    end loop;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double is

    res,absfz,denrco : double_double;
    fz : DoblDobl_Complex_Vectors.Vector(f'range);

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               fz : out QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double ) is

    one : constant quad_double := create(1.0);
    wrk,wrk_absfz,wrk_denrco : quad_double;

  begin
    Inverse_Condition_Number(f(f'first),z,fz(fz'first),absfz,denrco,rco);
    for i in f'first+1..f'last loop
      exit when rco + one = one;
      Inverse_Condition_Number(f(i),z,fz(i),wrk_absfz,wrk_denrco,wrk);
      if wrk < rco
       then rco := wrk; absfz := wrk_absfz; denrco := wrk_denrco;
      end if;
    end loop;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double is

    res,absfz,denrco : quad_double;
    fz : QuadDobl_Complex_Vectors.Vector(f'range);

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               fz : out Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number ) is

    wrk,wrk_absfz,wrk_denrco : Floating_Number;

  begin
    Inverse_Condition_Number(f(f'first),z,fz(fz'first),absfz,denrco,rco);
    for i in f'first+1..f'last loop
      Inverse_Condition_Number(f(i),z,fz(i),wrk_absfz,wrk_denrco,wrk);
      if wrk < rco 
       then Copy(wrk,rco); Copy(wrk_absfz,absfz); Copy(wrk_denrco,denrco);
      end if;
      Clear(wrk); Clear(wrk_absfz); Clear(wrk_denrco);
    end loop;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Multprec_Complex_Poly_Systems.Poly_Sys;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number is

    res,absfz,denrco : Floating_Number;
    fz : Multprec_Complex_Vectors.Vector(f'range);

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    Multprec_Complex_Vectors.Clear(fz);
    Clear(absfz); Clear(denrco);
    return res;
  end Inverse_Condition_Number;

-- PART IV : Laurent polynomial systems

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Laur_Systems.Laur_Sys;
               z : in Standard_Complex_Vectors.Vector;
               fz : out Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float ) is

    wrk,wrk_absfz,wrk_denrco : double_float;

  begin
    Inverse_Condition_Number(f(f'first),z,fz(fz'first),absfz,denrco,rco);
    for i in f'first+1..f'last loop
      exit when rco + 1.0 = 1.0;
      Inverse_Condition_Number(f(i),z,fz(i),wrk_absfz,wrk_denrco,wrk);
      if wrk < rco
       then rco := wrk; absfz := wrk_absfz; denrco := wrk_denrco;
      end if;
    end loop;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               fz : out DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double ) is

    one : constant double_double := create(1.0);
    wrk,wrk_absfz,wrk_denrco : double_double;

  begin
    Inverse_Condition_Number(f(f'first),z,fz(fz'first),absfz,denrco,rco);
    for i in f'first+1..f'last loop
      exit when rco + one = one;
      Inverse_Condition_Number(f(i),z,fz(i),wrk_absfz,wrk_denrco,wrk);
      if wrk < rco
       then rco := wrk; absfz := wrk_absfz; denrco := wrk_denrco;
      end if;
    end loop;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               fz : out QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double ) is

    one : constant quad_double := create(1.0);
    wrk,wrk_absfz,wrk_denrco : quad_double;

  begin
    Inverse_Condition_Number(f(f'first),z,fz(fz'first),absfz,denrco,rco);
    for i in f'first+1..f'last loop
      exit when rco + one = one;
      Inverse_Condition_Number(f(i),z,fz(i),wrk_absfz,wrk_denrco,wrk);
      if wrk < rco
       then rco := wrk; absfz := wrk_absfz; denrco := wrk_denrco;
      end if;
    end loop;
  end Inverse_Condition_Number;

  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Laur_Systems.Laur_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               fz : out Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number ) is

    wrk,wrk_absfz,wrk_denrco : Floating_Number;

  begin
    Inverse_Condition_Number(f(f'first),z,fz(fz'first),absfz,denrco,rco);
    for i in f'first+1..f'last loop
      Inverse_Condition_Number(f(i),z,fz(i),wrk_absfz,wrk_denrco,wrk);
      if wrk < rco 
       then Copy(wrk,rco); Copy(wrk_absfz,absfz); Copy(wrk_denrco,denrco);
      end if;
      Clear(wrk); Clear(wrk_absfz); Clear(wrk_denrco);
    end loop;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Standard_Complex_Laur_Systems.Laur_Sys;
               z : Standard_Complex_Vectors.Vector ) return double_float is

    res,absfz,denrco : double_float;
    fz : Standard_Complex_Vectors.Vector(f'range);

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double is

    res,absfz,denrco : double_double;
    fz : DoblDobl_Complex_Vectors.Vector(f'range);

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double is

    res,absfz,denrco : quad_double;
    fz : QuadDobl_Complex_Vectors.Vector(f'range);

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    return res;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Multprec_Complex_Laur_Systems.Laur_Sys;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number is

    res,absfz,denrco : Floating_Number;
    fz : Multprec_Complex_Vectors.Vector(f'range);

  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,res);
    Multprec_Complex_Vectors.Clear(fz);
    Clear(absfz); Clear(denrco);
    return res;
  end Inverse_Condition_Number;

-- PART V : the wrappers

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Polynomials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float;
               fz : out Standard_Complex_Numbers.Complex_Number ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Polynomials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double;
               fz : out DoblDobl_Complex_Numbers.Complex_Number ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Polynomials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double;
               fz : out QuadDobl_Complex_Numbers.Complex_Number ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Polynomials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number;
               fz : out Multprec_Complex_Numbers.Complex_Number ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Laurentials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float;
               fz : out Standard_Complex_Numbers.Complex_Number ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Laurentials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double;
               fz : out DoblDobl_Complex_Numbers.Complex_Number ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Laurentials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double;
               fz : out QuadDobl_Complex_Numbers.Complex_Number ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Laurentials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number;
               fz : out Multprec_Complex_Numbers.Complex_Number ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float;
               fz : out Standard_Complex_Vectors.Vector ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double;
               fz : out DoblDobl_Complex_Vectors.Vector ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double;
               fz : out QuadDobl_Complex_Vectors.Vector ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number;
               fz : out Multprec_Complex_Vectors.Vector ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Laur_Systems.Laur_Sys;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float;
               fz : out Standard_Complex_Vectors.Vector ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double;
               fz : out DoblDobl_Complex_Vectors.Vector ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double;
               fz : out QuadDobl_Complex_Vectors.Vector ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Laur_Systems.Laur_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number;
               fz : out Multprec_Complex_Vectors.Vector ) is
  begin
    Inverse_Condition_Number(f,z,fz,absfz,denrco,rco);
  end Evaluate_with_Inverse_Condition;

end VarbPrec_Polynomial_Evaluations;
