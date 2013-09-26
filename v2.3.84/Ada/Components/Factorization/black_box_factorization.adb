with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Drivers_to_Factor_Polynomials;     use Drivers_to_Factor_Polynomials;

procedure Black_Box_Factorization
            ( infilename : in string; file : in file_type; p : in Poly ) is

  n : constant natural32 := Number_of_Unknowns(p);
  fail : boolean;
  factors : Link_to_Poly_Sys;

begin
  put(file,"1 ");
  put(file,n,1);
  new_line(file);
  put(file,p);
  new_line(file);
  new_line(file);
  Driver_to_Factor(file,false,false,n,p,fail,factors);
  if factors /= null 
   then Write_Factors(infilename,factors.all);
  end if;
end Black_Box_Factorization;
