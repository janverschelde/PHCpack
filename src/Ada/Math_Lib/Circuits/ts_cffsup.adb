with text_io;                             use text_io;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Standard_Complex_VecVecs;
with Standard_Natural_Vectors_io;         use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with Coefficient_Supported_Polynomials;   use Coefficient_Supported_Polynomials;

procedure ts_cffsup is

-- DESCRIPTION :
--   Tests the operations in coefficient_supported_polynomials.

  procedure Extract_Supports_and_Coefficients
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    c : Standard_Complex_VecVecs.VecVec(p'range);
    e : Standard_Natural_VecVecs.Array_of_VecVecs(p'range);

  begin
    Coefficients_and_Supports(p,c,e);
    for i in p'range loop
      put("coefficients and supports of polynomial "); put(i,1);
      put_line(" : ");
      for j in c(i)'range loop
        put(c(i)(j)); put(" "); put(e(i)(j).all); new_line;
      end loop;
    end loop;
  end Extract_Supports_and_Coefficients;

  procedure Main is

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("reading a polynomial system ..."); get(p);
    put_line("-> your system : "); put(p.all);
    Extract_Supports_and_Coefficients(p.all);
  end Main;

begin
  Main;
end ts_cffsup;
