with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Natural_Vectors;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Linear_Product_System;

package body Standard_Linear_Product_System_io is

  procedure get ( n : out natural32 ) is
  begin
    Standard_Linear_Product_System_io.get(Standard_Input,n);
  end get;

  procedure get ( file : in file_type; n : out natural32 ) is

    nn : natural32;

  begin
    get(file,nn); n := nn;
    Standard_Linear_Product_System.Init(nn);
    declare
      h : Vector(0..integer32(nn));
      pp : Poly;
      d : degrees := new Standard_Natural_Vectors.Vector'
                           (1..integer32(nn) => 0);
      stop : boolean;
    begin
      get(file,nn,pp); Clear(pp);
      for i in 1..nn loop
        stop := false;
        while not stop loop
          get(file,nn,pp);
          stop := (pp = Null_Poly);
          exit when stop;
          h(0) := Coeff(pp,d);
          for j in 1..integer32(nn) loop
            d(j) := 1;
            h(j) := Coeff(pp,d);
            d(j) := 0;
          end loop;
          Standard_Linear_Product_System.Add_Hyperplane(i,h);
        end loop;
      end loop;
      Clear(d);
    end;
  end get;

  procedure put ( n,fore,after,exp : in natural32 ) is
  begin
    Standard_Linear_Product_System_io.put(Standard_Output,n,fore,after,exp);
  end put;

  procedure put ( file : in file_type; n,fore,after,exp : in natural32 ) is

    h : Vector(0..integer32(n));

    procedure Write_Number ( file : in file_type; x : in Complex_Number ) is
    begin
      if IMAG_PART(x) + 1.0 = 1.0 then
        put(file,REAL_PART(x),fore,after,exp);
      else
        put(file,'(');
        put(file,REAL_PART(x),fore,after,exp);
        put(file,'+');
        put(file,IMAG_PART(x),fore,after,exp);
        put(file,')');
      end if;
    end Write_Number;

  begin
    for i in 1..n loop
      put(file,"The hyperplanes for the "); put(file,i,1); 
      put_line(file,"th equation :");
      for j in 1..Standard_Linear_Product_System.Number_Of_Hyperplanes(i) loop
        h := Standard_Linear_Product_System.Get_Hyperplane(i,j);
        put(file,' ');
        for k in 1..n loop
  	  Write_Number(file,h(integer32(k)));
	  put(file,'*');
	  Symbol_Table_io.put(file,Symbol_Table.Get(k));
	  put(file," + ");
        end loop;
        Write_Number(file,h(0));
        new_line(file);
      end loop;
    end loop;
  end put;

end Standard_Linear_Product_System_io;
