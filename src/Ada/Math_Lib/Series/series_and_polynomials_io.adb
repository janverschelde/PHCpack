with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;
with Series_and_Polynomials;                use Series_and_Polynomials;

package body Series_and_Polynomials_io is

  procedure get ( s : out Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( file : in file_type; s : out Series ) is

    p : Standard_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    Standard_Complex_Polynomials_io.get(p);
    s := Polynomial_to_Series(p);
    Standard_Complex_Polynomials.Clear(p);
  end get;

  procedure put ( s : in Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( file : in file_type; s : in Series ) is

    p : Standard_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    Standard_Complex_Polynomials_io.put(file,p);
    Standard_Complex_Polynomials.Clear(p);
  end put;

  procedure get ( p : out Standard_Series_Polynomials.Poly ) is
  begin
    get(standard_input,p);
  end get;

  procedure get ( file : in file_type;
                  p : out Standard_Series_Polynomials.Poly ) is

    q : Standard_Complex_Polynomials.Poly;

  begin
    Standard_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q);
    Standard_Complex_Polynomials.Clear(q);
  end get;

  procedure put ( p : in Standard_Series_Polynomials.Poly ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type;
                  p : in Standard_Series_Polynomials.Poly ) is

    q : Standard_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p);

  begin
    Standard_Complex_Polynomials_io.put(file,q);
    Standard_Complex_Polynomials.Clear(q);
  end put;

end Series_and_Polynomials_io;
