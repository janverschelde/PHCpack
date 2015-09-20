with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with Drivers_to_Factor_Polynomials;     use Drivers_to_Factor_Polynomials;

package body Black_Box_Factorization is

  procedure Standard_Black_Box_Factorization
              ( infilename : in string; file : in file_type;
                p : in Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Poly_Systems;

    n : constant natural32
      := Standard_Complex_Polynomials.Number_of_Unknowns(p);
    fail : boolean;
    factors : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

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
  end Standard_Black_Box_Factorization;

  procedure DoblDobl_Black_Box_Factorization
              ( infilename : in string; file : in file_type;
                p : in DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Poly_Systems;

    n : constant natural32
      := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p);
    fail : boolean;
    factors : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

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
  end DoblDobl_Black_Box_Factorization;

  procedure QuadDobl_Black_Box_Factorization
              ( infilename : in string; file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Poly_Systems;

    n : constant natural32
      := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p);
    fail : boolean;
    factors : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

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
  end QuadDobl_Black_Box_Factorization;

end Black_Box_Factorization;
