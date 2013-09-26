with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;          use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Complex_Polynomial_Matrices;       use Complex_Polynomial_Matrices;
with Complex_Polynomial_Matrices_io;    use Complex_Polynomial_Matrices_io;

procedure ts_cpm is

  procedure Main is

    lpm : Link_to_Polynomial_Matrix;
  
  begin
    Interactive_get(lpm);
    put_line("Your polynomial matrix : ");
    put(lpm.all);
    declare
      d : constant Standard_Integer_Vectors.Vector := Degrees(lpm.all);
      k : constant integer32 := Sum(d) + d'length;
      c : constant Standard_Complex_Vectors.Vector := Coefficients(k,lpm.all);
      n : constant integer32 := lpm'length(1);
      m : constant integer32 := lpm'length(2);
    begin
      put("The degrees : "); put(d); new_line;
      put("The total number of coefficients : "); put(k,1); new_line;
      put_line("The coefficients in the polynomial matrix : ");
      put_line(c);
      put_line("The recreated polynomial matrix : ");
      put(Create(n,m,d,c));
    end;
  end Main;

begin
  new_line;
  put_line("Representations of matrix of complex univariate polynomials.");
  new_line;
  Main;
end ts_cpm;
