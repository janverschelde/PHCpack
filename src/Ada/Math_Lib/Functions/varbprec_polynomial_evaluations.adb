with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body VarbPrec_Polynomial_Evaluations is

  function Inverse_Condition_Number
             ( f : Standard_Complex_Polynomials.Poly;
               z : Standard_Complex_Vectors.Vector ) return double_float is

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

end VarbPrec_Polynomial_Evaluations;
