with integer_io;                         use integer_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Bracket_Monomials_io;               use Bracket_Monomials_io;

package body QuadDobl_Bracket_Polynomials_io is

  procedure Write_Number ( file : in file_type; c : in Complex_Number;
                           ic : out integer ) is

  -- DESCRIPTION :
  --   Writes the complex coefficient in a concise way.

    f : quad_double;
    i : integer := 0;
    zero : constant quad_double := create(0.0);

  begin
    if IMAG_PART(c) /= zero then
      put(file,"+(");
      put(file,REAL_PART(c)); 
      if IMAG_PART(c) > zero
       then put(file," +"); put(file,IMAG_PART(c));
       else put(file," -"); put(file,-IMAG_PART(c));
      end if;
      put(file,")");
    else
      f := REAL_PART(c);
      i := integer(hihi_part(f));
      if (f - Quad_Double_Numbers.Create(i)) /= zero then
        put(file,f);
      else
        if i > 0 then
          put(file,"+");
          if i /= 1
           then put(file,i,1);
          end if;
        elsif i = -1 then
          put(file,"-");
        else
          put(file,i,1);
        end if;
      end if;
    end if;
    ic := i;
  end Write_Number;

  procedure put ( t : in Bracket_Term ) is
  begin
    put(Standard_Output,t);
  end put;

  procedure put ( file : in file_type; t : in Bracket_Term ) is

    ic : integer;

  begin
    Write_Number(file,t.coeff,ic);
    if ic /= 1 and ic /= -1
     then put(file,"*");
    end if;
    put(file,t.monom);
  end put;

  procedure put ( p : in Bracket_Polynomial ) is
  begin
    put(Standard_Output,p);
  end put;

  procedure put ( file : in file_type; p : in Bracket_Polynomial ) is

    cnt : natural := 0;

    procedure Write_Term ( t : in Bracket_Term; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      if cnt > 4
       then cnt := 1;
            new_line(file);
      end if;
      put(file,t);
      continue := true;
    end Write_Term;

    procedure Write_Terms is new Enumerate_Terms(Write_Term);

  begin
    Write_Terms(p);
    put_line(file,";");
  end put;

  procedure put_line ( p : in Bracket_Polynomial ) is
  begin
    put_line(Standard_Output,p);
  end put_line;

  procedure put_line ( file : in file_type; p : in Bracket_Polynomial ) is

    procedure Write_Term ( t : in Bracket_Term; continue : out boolean ) is
    begin
      new_line(file);
      put(file,t);
      continue := true;
    end Write_Term;

    procedure Write_Terms is new Enumerate_Terms(Write_Term);

  begin
    Write_Terms(p);
    put_line(file,";");
  end put_line;

end QuadDobl_Bracket_Polynomials_io;
