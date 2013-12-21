with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Random_Numbers;
with Multprec_Natural_Coefficients;     use Multprec_Natural_Coefficients;

package body Extended_Random_Numbers is

  function Extended_Random
             ( n : Natural_Number; size : natural32 ) return Natural_Number is

    res : Natural_Number;
    the_base : constant natural32 := Multprec_Natural_Coefficients.Base;
    cff : constant Array_of_Naturals := Coefficients(n);
    arn : Array_of_Naturals(0..size);
    shift : natural32;

  begin
    if size >= cff'last then
      shift := size - cff'last;
      for i in cff'range loop
        exit when (i > size);
        arn(i+shift) := cff(i);
      end loop;
    else
      shift := cff'last - size;
      for i in arn'range loop
        arn(i) := cff(i+shift);
      end loop;
      shift := 0;
    end if;
    if shift > 0 then
      for i in 0..shift-1 loop
        arn(i) := natural32(Standard_Random_Numbers.Random
                    (integer32(0),integer32(the_base-1)));
      end loop;
    end if;
    res := Create(arn);
    return res;
  end Extended_Random;

  function Extended_Random
             ( i : Integer_Number; size : natural32 ) return Integer_Number is

    res : Integer_Number;
    nat : constant Natural_Number := Unsigned(i);
    res_nat : constant Natural_Number := Extended_Random(nat,size);

  begin
    res := Shallow_Create(res_nat);
    if Negative(i)
     then Min(res);
    end if;
    return res;
  end Extended_Random;

  function Extended_Random
             ( f : Floating_Number; size : natural32 )
             return Floating_Number is

    res : Floating_Number;
    frac : constant Integer_Number := Fraction(f);
    res_frac,res_expo : Integer_Number;
    szf : constant natural32 := Size_Fraction(f);
    df : integer32;

  begin
    res_frac := Extended_Random(frac,size);
    Copy(Exponent(f),res_expo);
    if size > szf then
      df := integer32(size) - integer32(szf);
      Sub(res_expo,df*integer32(Multprec_Natural_Coefficients.Exponent));
    elsif size < szf then
      df := integer32(szf) - integer32(size);
      Add(res_expo,df*integer32(Multprec_Natural_Coefficients.Exponent));
    end if;
    res := Create(res_frac,res_expo);
    return res;
  end Extended_Random;

  function Extended_Random
             ( c : Complex_Number; size : natural32 ) return Complex_Number is

    res : Complex_Number := c;
    re : Floating_Number := Real_Part(c);
    im : Floating_Number := Imag_Part(c);
    extre : Floating_Number := Extended_Random(re,size);
    extim : Floating_Number := Extended_Random(im,size);

  begin
    res := Create(extre,extim);
    Clear(re); Clear(im);
    Clear(extre); Clear(extim);
    return res;
  end Extended_Random;

end Extended_Random_Numbers;
