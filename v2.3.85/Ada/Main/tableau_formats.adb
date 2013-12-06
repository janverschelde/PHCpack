with integer_io;                         use integer_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Symbol_Table,Symbol_Table_io;       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Power_Lists;                        use Power_Lists;

package body Tableau_Formats is

-- AUXILIARIES :

  procedure Read_Symbols ( file : in file_type; n : in natural ) is

  -- DESCRIPTION :
  --   Reads n symbols from file.

  begin
    Symbol_Table.Init(n);
    for i in 1..n loop
      declare
        s : Symbol;
      begin
        Symbol_Table_io.get(file,s);
        Symbol_Table.Add(s);
      end;
    end loop;
  end Read_Symbols;

  procedure Write_Symbols ( file : in file_type ) is

  -- DESCRIPTION :
  --   Writes the content of the symbol table on file.

    n : constant natural := Symbol_Table.Number;

  begin
    put(file,n,1); new_line(file);
    for i in 1..n loop
      declare
        s : Symbol := Symbol_Table.Get(i);
      begin
        put(file," ");
        Symbol_Table_io.put(file,s);
      end;
    end loop;
    new_line(file);
  end Write_Symbols;

  procedure Read_Exponents ( file : in file_type;
                             exp : in out Array_of_Lists ) is

  -- DESCRIPTION :
  --   Reads exponent vectors from file.
  --   Each list of exponent vectors should be preceded by the number
  --   of elements in the list.

    n : constant natural := exp'length;
    m : natural;

  begin
    for i in exp'range loop
      get(file,m);
      get(file,n,m,exp(i));
    end loop;
  end Read_Exponents;

  procedure Write_Exponents ( file : in file_type; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Writes the supports of p on file.

    exp : Array_of_Lists(p'range) := Create(p);

  begin
    for i in exp'range loop
      put(file,Length_Of(exp(i)),1); new_line(file);
      put(file,exp(i));
    end loop;
    Deep_Clear(exp);
  end Write_Exponents;

  procedure Read_Coefficients
               ( file : in file_type;
                 exp : in Array_of_Lists; flt : in boolean;
                 cff : in out VecVec ) is

  -- DESCRIPTION :
  --   Reads the coefficient from file, given the supports.

    f : double_float;

  begin
    for i in cff'range loop
      cff(i) := new Standard_Complex_Vectors.Vector(1..Length_Of(exp(i)));
      for j in cff(i)'range loop
        if flt
         then get(file,f);         cff(i)(j) := Create(f);
         else get(file,cff(i)(j));
        end if;
      end loop;
    end loop;
  end Read_Coefficients;

  procedure Write_Coefficients ( file : in file_type; flt : in boolean;
                                 p : in Poly ) is

  -- DESCRIPTION :
  --   Writes the coefficients of the polynomial on file.
  --   If flt, then the imaginary parts will be omitted.

    procedure Write_Coefficient_of_Term ( t : in Term; cont : out boolean ) is
    begin
      if flt
       then put(file,REAL_PART(t.cf));
       else put(file,t.cf);
      end if;
      new_line(file);
      cont := true;
    end Write_Coefficient_of_Term;

    procedure Write_Coefficients_of_Terms is
      new Visiting_Iterator(Write_Coefficient_of_Term);

  begin
    Write_Coefficients_of_Terms(p);
  end Write_Coefficients;

  procedure Write_Coefficients ( file : in file_type; flt : in boolean;
                                 p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Writes the coefficients of the polynomial system on file.
  --   If flt, then the imaginary parts will be omitted.

  begin
    for i in p'range loop
      Write_Coefficients(file,flt,p(i));
    end loop;
  end Write_Coefficients;

  function Create ( exp : Standard_Integer_Vectors.Vector;
                    cff : Complex_Number ) return Term is

  -- DESCRIPTION :
  --   Creates a term from exponent vector and coefficient.

    res : Term;

  begin
    res.cf := cff;
    res.dg := new Standard_Natural_Vectors.Vector(exp'range);
    for i in exp'range loop
      res.dg(i) := exp(i);
    end loop;
    return res;
  end Create;

  function Create ( exp : List; cff : Standard_Complex_Vectors.Vector )
                  return Poly is

  -- DESCRIPTION :
  --   Creates a polynomial from the exponents and the coefficients.

    res : Poly;
    tmp : List := exp;
    lpt : Standard_Integer_Vectors.Link_to_Vector;

  begin
    for i in cff'range loop
      lpt := Head_Of(tmp);
      declare
        t : Term := Create(lpt.all,cff(i));
      begin
        Add(res,t);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  function Create ( exp : Array_of_Lists; cff : VecVec ) return Poly_Sys is

  -- DESCRIPTION :
  --   Creates the polynomial system from the supports and coefficients.

    res : Poly_Sys(exp'range);

  begin
    for i in res'range loop
      res(i) := Create(exp(i),cff(i).all);
    end loop;
    return res;
  end Create;

-- TARGET ROUTINES :

  procedure get ( file : in file_type; realcoeff : in boolean;
                  p : out Link_to_Poly_Sys ) is

    n : natural;

  begin
    get(file,n);
    new_line;
    put("Dimension : "); put(n,1); new_line;
    Read_Symbols(file,n);
    put("Symbols :"); 
    for i in 1..n loop
      put(" "); put(Symbol_Table.get(i));
    end loop;
    new_line;
    declare
      exp : Array_of_Lists(1..n);
      cff : VecVec(1..n);
      sys : Poly_Sys(1..n);
    begin
      Read_Exponents(file,exp);
      put_line("The exponents : ");
      put(exp);
      Read_Coefficients(file,exp,realcoeff,cff);
      put_line("The coefficients : ");
      for i in cff'range loop
        for j in cff(i)'range loop
          put(cff(i)(j)); new_line;
        end loop;
        new_line;
      end loop;
      sys := Create(exp,cff);
      put_line("The polynomial system : "); put(sys);
      p := new Poly_Sys'(sys);
    end;
  end get;

  procedure put ( file : in file_type; realcoeff : in boolean; 
                  p : in Poly_Sys ) is

  begin
    if not Symbol_Table.Empty
     then Write_Symbols(file);
    end if;
    Write_Exponents(file,p);
    Write_Coefficients(file,realcoeff,p);
  end put;

end Tableau_Formats;
