with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Double_Laurent_Series;             use Double_Laurent_Series;

package body Test_Double_Laurent_Series is

  procedure Test_Multiply_Inverse_Divide ( deg : in integer32 ) is

    ale,ble,prodle,invble,quotle : integer32 := 0;
    a : constant Standard_Complex_Vectors.Vector(0..deg)
      := Standard_Random_Vectors.Random_Vector(0,deg);
    b : constant Standard_Complex_Vectors.Vector(0..deg)
      := Standard_Random_Vectors.Random_Vector(0,deg);
    prod : Standard_Complex_Vectors.Vector(0..deg);
    invb : Standard_Complex_Vectors.Vector(0..deg);
    quot : Standard_Complex_Vectors.Vector(0..deg);

  begin
    new_line;
    put("Give the leading exponent of a : "); get(ale);
    put("Give the leading exponent of b : "); get(ble);
    new_line;
    put_line("A random series a :"); Write(ale,a);
    new_line;
    put_line("A random series b :"); Write(ble,b);
    Multiply(deg,ale,ble,a,b,prodle,prod);
    new_line;
    put_line("The product of a with b :"); Write(prodle,prod);
    Divide(deg,prodle,ble,prod,b,quotle,quot,invb);
    new_line;
    put_line("The quotient of a*b with b :"); Write(quotle,quot);
    invble := -ble;
    Multiply(deg,ble,invble,b,invb,prodle,prod);
    put_line("The product of b with 1/b :"); Write(prodle,prod);
  end Test_Multiply_Inverse_Divide;

  procedure Test_Add_and_Subtract ( deg : in integer32 ) is

    ale,ble,sumle,difle : integer32 := 0;
    a : constant Standard_Complex_Vectors.Vector(0..deg)
      := Standard_Random_Vectors.Random_Vector(0,deg);
    b : constant Standard_Complex_Vectors.Vector(0..deg)
      := Standard_Random_Vectors.Random_Vector(0,deg);
    sum : Standard_Complex_Vectors.Vector(0..deg);
    dif : Standard_Complex_Vectors.Vector(0..deg);

  begin
    new_line;
    put("Give the leading exponent of a : "); get(ale);
    put("Give the leading exponent of b : "); get(ble);
    new_line;
    put_line("A random series a :"); Write(ale,a);
    new_line;
    put_line("A random series b :"); Write(ble,b);
    Add(deg,ale,ble,a,b,sumle,sum);
    new_line;
    put_line("The sum of a and b :"); Write(sumle,sum);
    Subtract(deg,sumle,ble,sum,b,difle,dif);
    new_line;
    put_line("The result of (a + b) - b :"); Write(difle,dif);
    Subtract(deg,ale,ale,a,a,difle,dif);
    new_line;
    put_line("The result of a - a :"); Write(difle,dif);
  end Test_Add_and_Subtract;

  procedure Main is

     d : integer32 := 0;

  begin
    new_line;
    put("Give the truncation degree : "); get(d);
    Test_Multiply_Inverse_Divide(d);
    Test_Add_and_Subtract(d);
  end Main;

end Test_Double_Laurent_Series;
