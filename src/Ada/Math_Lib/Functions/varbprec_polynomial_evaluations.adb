with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;

package body VarbPrec_Polynomial_Evaluations is

  function Inverse_Condition_Number
             ( f : Standard_Complex_Polynomials.Poly;
               z : Standard_Complex_Vectors.Vector ) return double_float is

    use Standard_Complex_Numbers;

    value : Complex_Number := create(0.0);
    absum : double_float := 0.0;

    use Standard_Complex_Polynomials;

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
    return AbsVal(value)/absum;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Polynomials.Poly;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double is

    use DoblDobl_Complex_Numbers;

    zero : constant double_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : double_double := create(0.0);

    use DoblDobl_Complex_Polynomials;

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
    return AbsVal(value)/absum;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Polynomials.Poly;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double is

    use QuadDobl_Complex_Numbers;

    zero : constant quad_double := create(0.0);
    value : Complex_Number := create(zero);
    absum : quad_double := create(0.0);

    use QuadDobl_Complex_Polynomials;

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
    return AbsVal(value)/absum;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number
             ( f : Multprec_Complex_Polynomials.Poly;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number is

    use Multprec_Complex_Numbers;

    res : Floating_Number;
    zero : Floating_Number := create(0.0);
    value : Complex_Number := create(zero);
    absum : Floating_Number := create(0.0);

    use Multprec_Complex_Polynomials;

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
    Clear(value); Clear(absum); Clear(zero);
    return res;
  end Inverse_Condition_Number;

end VarbPrec_Polynomial_Evaluations;
