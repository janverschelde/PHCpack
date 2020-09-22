with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with TripDobl_Complex_Numbers;           use TripDobl_Complex_Numbers;
with Symbol_Table,Symbol_Table_io;
with TripDobl_Polynomial_Convertors;     use TripDobl_Polynomial_Convertors;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laurentials_io;

package body TripDobl_Complex_Laurentials_io is

-- THE INPUT OPERATIONS :

  procedure get ( p : out Poly ) is
  begin
    get(standard_input,p);
  end get;

  procedure get ( file : in file_type; p : out Poly ) is

    mp : Multprec_Complex_Laurentials.Poly;

  begin
    Multprec_Complex_Laurentials_io.Set_Working_Precision(7);
    Multprec_Complex_Laurentials_io.get(file,mp);
    p := Multprec_Laurential_to_TripDobl_Complex(mp);
    Multprec_Complex_Laurentials.Clear(mp);
  end get;

-- AUXILIARIES FOR OUTPUT ROUTINES :

  function All_Zeroes ( d : Degrees ) return boolean is

  -- DESCRIPTION :
  --   Returns true if all d(i) = 0, for all i, otherwise false is returned.

  begin
    for i in d'range loop
      if d(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end All_Zeroes;

  procedure Write_Number ( file : in file_type; c : in Complex_Number ) is

  -- DESCRIPTION :
  --   Writes the triple double complex number to file
  --   in the format ( real_part(c) + imag_part(c)*I ).

    re : constant triple_double := REAL_PART(c);
    im : constant triple_double := IMAG_PART(c);

  begin
    put(file,"+ ("); put(file,re);
    if im >= 0.0
     then put(file,"+");
    end if;
    put(file,im);
    put(file,"*I)");
  end Write_Number;

  procedure Write_Factor ( file : in file_type; d,i : in integer32;
                           standard : in boolean; pow : in character ) is

  -- DESCRIPTION :
  --   Writes the factor corresponding with the ith unknown on file.

    sb : Symbol_Table.Symbol;

  begin
    if standard then
      put(file,'x');
      if i < 10
       then put(file,i,1);
       else put(file,i,2);
      end if;
    else 
      sb := Symbol_Table.get(natural32(i)); Symbol_Table_io.put(file,sb);
    end if;
    if d > 1 or d < 0 then
      if pow = '^'
       then put(file,'^');
       else put(file,"**");
      end if;
      if d < 10
       then put(file,d,1);
       else put(file,d,2);
      end if;
    end if;
  end Write_Factor;

-- THE OUTPUT OPERATIONS :

  procedure put ( p : in Poly ) is
  begin
    put_line(standard_output,p);
  end put;
 
  procedure put ( file : in file_type; p : in Poly ) is
  begin
    put_line(file,p);
  end put;

  procedure put_line ( p : in Poly ) is
  begin
    put_line(standard_output,p);
  end put_line;

  procedure put_line ( file : in file_type; p : in Poly ) is

    n : constant natural32 := Number_of_Unknowns(p);
    standard : constant boolean := ( Symbol_Table.Number < n );

    procedure Write_Term ( t : in Term; continue : out boolean ) is
 
    -- DESCRIPTION : 
    --   Writes a term to file.

    begin
      new_line(file);
      Write_Number(file,t.cf);
      if not All_Zeroes(t.dg) then
        for i in t.dg'range loop
          if t.dg(i) > 0 or t.dg(i) < 0 then
            put(file,'*');
            Write_Factor(file,t.dg(i),i,standard,'^');
          end if;
        end loop;
      end if;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
    put_line(file,";");
  end put_line;

end TripDobl_Complex_Laurentials_io;
