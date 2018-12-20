with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;

package body QuadDobl_Complex_Monomials_io is

  procedure put ( file : in file_type; m : in Monomial ) is
  begin
    put(file,"coefficient : "); put(file,m.cff); new_line(file);
    put(file,"dimension : "); put(file,m.dim,1); new_line(file);
    put(file,"number of variables : "); put(file,m.nvr,1); new_line(file);
    put(file,"number of variables with exponent > 1 : ");
    put(file,m.n_base,1); new_line(file);
    put(file,"degree : "); put(file,Degree(m),1); new_line(file);
    put(file,"largest exponent : ");
    put(file,Largest_Exponent(m),1); new_line(file);
    if m.nvr > 0 then
      put(file,"positions : "); put(file,m.pos); new_line(file);
      put(file,"exponents : "); put(file,m.exp); new_line(file);
      if m.n_base > 0 then
        put_line(file,"for variables with exponents > 1 :");
        put(file,"positions : "); put(file,m.pos_base); new_line(file);
        put(file,"exponents : "); put(file,m.exp_base); new_line(file);
        put(file,"exponents minus 2 : ");
        put(file,m.exp_tbl_base); new_line(file);
      end if;
    end if;
  end put;

  procedure put ( m : in Monomial ) is
  begin
    put(standard_output,m);
  end put;

end QuadDobl_Complex_Monomials_io;
