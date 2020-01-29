with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;

procedure ts_hesscnv is

-- DESCRIPTION :
--   Computing the value of the Hessian matrix of a convolution circuit.

  function Diff ( x : Standard_Complex_Vectors.Vector;
                  e : Standard_Integer_Vectors.Vector; i,j : integer32 )  
                return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the value of the second derivative of the monomial
  --   with exponents in e, with respect to i and j, evaluated at x.

  -- REQUIRED : i is in e'range and j is in e'range, 
  --   i and j can be the same if the second derivative
  --   with the same i is requested.

    use Standard_Complex_Numbers;

    res : Complex_Number := Create(0.0);
    fac : double_float;

  begin
    if i = j then
      if e(i) >= 2 then
        fac := double_float(e(i))*double_float(e(i)-1);
        res := Create(fac);
        for k in 1..e(i)-2 loop
          res := res*x(i);
        end loop;
        for k in e'range loop
          if k /= i then
            for j in 1..e(k) loop
              res := res*x(k);
            end loop;
          end if;
        end loop;
      end if;
    elsif e(i) < 1 or e(j) < 0 then
      null;
    end if;
    return res;
  end Diff;

  procedure Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Uses the dimension, degree, number of terms, and largest power
  --   to generate a random circuit for testing.

    use Standard_Speelpenning_Convolutions;

    c : constant Circuit(nbr,dim,dim-1,dim-2)
      := Standard_Random_Convolution_Circuit(dim,deg,nbr,pwr);
    d : Standard_Complex_Numbers.Complex_Number;
    x : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    idx1,idx2 : integer32 := 0;

  begin
    put_line("The exponents : ");
    for i in c.xps'range loop
      put(c.xps(i).all); new_line;
    end loop;
    put("Give the first index : "); get(idx1);
    put("Give the second index : "); get(idx2);
    d := Diff(x,c.xps(1).all,idx1,idx2);
    put("The value of the second derivative w.r.t. ");
    put(idx1,1); put(" and "); put(idx2,1); put_line(" :");
    put(d); new_line;
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension, degree, number of terms,
  --   and largest power to generate a random circuit.

    dim,deg,nbr,pwr : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    put("Give the number of terms : "); get(nbr);
    put("Give the largest power : "); get(pwr);
    Test(dim,deg,nbr,pwr);
  end Main;

begin
  Main;
end ts_hesscnv;
