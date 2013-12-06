with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;

package body Standard_Tableau_Formats is

  procedure Write_Tableau
              ( file : in file_type;
                t : in Standard_Complex_Polynomials.Term ) is

    f : double_float;

  begin
    f := REAL_PART(t.cf); put(file,f);
    f := IMAG_PART(t.cf); put(file,f);
    for i in t.dg'range loop
      put(file," "); put(file,t.dg(i),1);
    end loop;
    new_line(file);
  end Write_Tableau;

  procedure Write_Tableau
              ( file : in file_type;
                t : in Standard_Complex_Laurentials.Term ) is

    f : double_float;

  begin
    f := REAL_PART(t.cf); put(file,f);
    f := IMAG_PART(t.cf); put(file,f);
    for i in t.dg'range loop
      put(file," "); put(file,t.dg(i),1);
    end loop;
    new_line(file);
  end Write_Tableau;

  procedure Write_Tableau
              ( file : in file_type;
                p : in Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Polynomials;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      Write_Tableau(file,t);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
  end Write_Tableau;

  procedure Write_Tableau
              ( file : in file_type;
                p : in Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      Write_Tableau(file,t);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
  end Write_Tableau;

  procedure Write_Tableau ( file : in file_type; p : in Poly_Sys ) is

    use Standard_Complex_Polynomials;

  begin
    put(file,p'last,1); new_line(file);
    for i in p'range loop
      put(file,Number_of_Terms(p(i)),1); new_line(file);
      Write_Tableau(file,p(i));
    end loop;
  end Write_Tableau;

  procedure Write_Tableau ( file : in file_type; p : in Laur_Sys ) is

    use Standard_Complex_Laurentials;

  begin
    put(file,p'last,1); new_line(file);
    for i in p'range loop
      put(file,Number_of_Terms(p(i)),1); new_line(file);
      Write_Tableau(file,p(i));
    end loop;
  end Write_Tableau;

  procedure Convert_Polynomial_into_Tableau_Format is

    lp : Link_to_Poly_Sys;
    outfile : file_type;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(outfile);
    Write_Tableau(outfile,lp.all);
  end Convert_Polynomial_into_Tableau_Format;

  procedure Convert_Laurent_into_Tableau_Format is

    lp : Link_to_Laur_Sys;
    outfile : file_type;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(outfile);
    Write_Tableau(outfile,lp.all);
  end Convert_Laurent_into_Tableau_Format;

  procedure Read_Tableau
              ( file : in file_type; n : in natural32;
                t : in out Standard_Complex_Polynomials.Term ) is

    f_re,f_im : double_float := 0.0;

  begin
    get(file,f_re); get(file,f_im);
    t.cf := Create(f_re,f_im);
    for i in 1..integer32(n) loop
      get(file,t.dg(i));
    end loop;
  end Read_Tableau;

  procedure Read_Tableau
              ( file : in file_type; n : in natural32;
                t : in out Standard_Complex_Laurentials.Term ) is

    f_re,f_im : double_float := 0.0;

  begin
    get(file,f_re); get(file,f_im);
    t.cf := Create(f_re,f_im);
    for i in 1..integer32(n) loop
      get(file,t.dg(i));
    end loop;
  end Read_Tableau;

  procedure Read_Tableau
              ( file : in file_type; n : in natural32;
                p : out Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Polynomials;

    t : Term;
    nt : natural32 := 0;

  begin
    t.cf := Create(0.0);
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    get(file,nt);
    for i in 1..nt loop
      Read_Tableau(file,n,t);
      Add(p,t);
    end loop;
    Clear(t);
  end Read_Tableau;

  procedure Read_Tableau
              ( file : in file_type; n : in natural32;
                p : out Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;

    t : Term;
    nt : natural32 := 0;

  begin
    t.cf := Create(0.0);
    t.dg := new Standard_Integer_Vectors.Vector'(1..integer32(n) => 0);
    get(file,nt);
    for i in 1..nt loop
      Read_Tableau(file,n,t);
      Add(p,t);
    end loop;
    Clear(t);
  end Read_Tableau;

  procedure Read_Tableau
              ( file : in file_type; n : in natural32; p : out Poly_Sys ) is
  begin
    for i in p'range loop
      Read_Tableau(file,n,p(i));
    end loop;
  end Read_Tableau;

  procedure Read_Tableau
              ( file : in file_type; n : in natural32; p : out Laur_Sys ) is
  begin
    for i in p'range loop
      Read_Tableau(file,n,p(i));
    end loop;
  end Read_Tableau;

  procedure Read_Tableau
              ( file : in file_type; p : out Link_to_Poly_Sys ) is

    n : natural32 := 0;

  begin
    get(file,n);
    declare
      s : Poly_Sys(1..integer32(n));
    begin
      Read_Tableau(file,n,s);
      p := new Poly_Sys'(s);
    end;
  end Read_Tableau;

  procedure Read_Tableau
              ( file : in file_type; p : out Link_to_Laur_Sys ) is

    n : natural32 := 0;

  begin
    get(file,n);
    declare
      s : Laur_Sys(1..integer32(n));
    begin
      Read_Tableau(file,n,s);
      p := new Laur_Sys'(s);
    end;
  end Read_Tableau;

  procedure Extract_Coefficients_and_Exponents
              ( p : in Standard_Complex_Polynomials.Poly;
                c : out Standard_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec ) is

    use Standard_Complex_Polynomials;

    ind : integer32 := 0;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      c(ind) := t.cf;
      e(ind) := Standard_Natural_Vectors.Link_to_Vector(t.dg);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
  end Extract_Coefficients_and_Exponents;

  procedure Extract_Coefficients_and_Exponents_Copies
              ( p : in Standard_Complex_Polynomials.Poly;
                c : out Standard_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec ) is

    use Standard_Complex_Polynomials;

    ind : integer32 := 0;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      et : Standard_Natural_Vectors.Vector(t.dg'range);

    begin
      ind := ind + 1;
      c(ind) := t.cf;
      for i in t.dg'range loop
        et(i) := t.dg(i);
      end loop;
      e(ind) := new Standard_Natural_Vectors.Vector'(et);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
  end Extract_Coefficients_and_Exponents_Copies;

  procedure Extract_Coefficients_and_Exponents
              ( p : in Standard_Complex_Laurentials.Poly;
                c : out Standard_Complex_Vectors.Vector;
                e : out Standard_Integer_VecVecs.VecVec ) is

    use Standard_Complex_Laurentials;

    ind : integer32 := 0;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      c(ind) := t.cf;
      e(ind) := Standard_Integer_Vectors.Link_to_Vector(t.dg);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
  end Extract_Coefficients_and_Exponents;

  procedure Extract_Coefficients_and_Exponents_Copies
              ( p : in Standard_Complex_Laurentials.Poly;
                c : out Standard_Complex_Vectors.Vector;
                e : out Standard_Integer_VecVecs.VecVec ) is

    use Standard_Complex_Laurentials;

    ind : integer32 := 0;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      et : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      ind := ind + 1;
      c(ind) := t.cf;
      for i in et'range loop
        et(i) := t.dg(i);
      end loop;
      e(ind) := new Standard_Integer_Vectors.Vector'(et);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
  end Extract_Coefficients_and_Exponents_Copies;

  procedure Extract_Coefficients_and_Exponents
              ( p : in Poly_Sys;
                c : out Standard_Complex_VecVecs.VecVec;
                e : out Standard_Natural_VecVecs.Array_of_VecVecs ) is

    m : natural32;

  begin
    for i in p'range loop
      m := Standard_Complex_Polynomials.Number_of_Terms(p(i));
      c(i) := new Standard_Complex_Vectors.Vector(1..integer32(m)); 
      e(i) := new Standard_Natural_VecVecs.VecVec(1..integer32(m)); 
      Extract_Coefficients_and_Exponents(p(i),c(i).all,e(i).all);
    end loop;
  end Extract_Coefficients_and_Exponents;

  procedure Extract_Coefficients_and_Exponents_Copies
              ( p : in Poly_Sys;
                c : out Standard_Complex_VecVecs.VecVec;
                e : out Standard_Natural_VecVecs.Array_of_VecVecs ) is

    m : natural32;

  begin
    for i in p'range loop
      m := Standard_Complex_Polynomials.Number_of_Terms(p(i));
      c(i) := new Standard_Complex_Vectors.Vector(1..integer32(m)); 
      e(i) := new Standard_Natural_VecVecs.VecVec(1..integer32(m)); 
      Extract_Coefficients_and_Exponents_Copies(p(i),c(i).all,e(i).all);
    end loop;
  end Extract_Coefficients_and_Exponents_Copies;

  procedure Extract_Coefficients_and_Exponents
              ( p : in Laur_Sys;
                c : out Standard_Complex_VecVecs.VecVec;
                e : out Standard_Integer_VecVecs.Array_of_VecVecs ) is

    m : natural32;

  begin
    for i in p'range loop
      m := Standard_Complex_Laurentials.Number_of_Terms(p(i));
      c(i) := new Standard_Complex_Vectors.Vector(1..integer32(m)); 
      e(i) := new Standard_Integer_VecVecs.VecVec(1..integer32(m)); 
      Extract_Coefficients_and_Exponents(p(i),c(i).all,e(i).all);
    end loop;
  end Extract_Coefficients_and_Exponents;

  procedure Extract_Coefficients_and_Exponents_Copies
              ( p : in Laur_Sys;
                c : out Standard_Complex_VecVecs.VecVec;
                e : out Standard_Integer_VecVecs.Array_of_VecVecs ) is

    m : natural32;

  begin
    for i in p'range loop
      m := Standard_Complex_Laurentials.Number_of_Terms(p(i));
      c(i) := new Standard_Complex_Vectors.Vector(1..integer32(m)); 
      e(i) := new Standard_Integer_VecVecs.VecVec(1..integer32(m)); 
      Extract_Coefficients_and_Exponents_Copies(p(i),c(i).all,e(i).all);
    end loop;
  end Extract_Coefficients_and_Exponents_Copies;

  procedure Convert_Tableau_into_Symbolic_Format is

    infile,outfile : file_type;
    sys,lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading the name of the input file for the tableau...");
    Read_Name_and_Open_File(infile);
    new_line;
    put_line("Reading a polynomial system to initialize the symbol table...");
    get(sys);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(outfile);
    Read_Tableau(infile,lp);
    put_pair(outfile,lp.all);
    Clear(sys); Clear(lp);
  end Convert_Tableau_into_Symbolic_Format;

  procedure Convert_Polynomial_System_to_Tableau_Structures is

    p : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(p);
    declare
      c : Standard_Complex_VecVecs.VecVec(p'range);
      e : Standard_Natural_VecVecs.Array_of_VecVecs(p'range);
    begin
      Extract_Coefficients_and_Exponents(p.all,c,e);
      new_line;
      put_line("The tableau structure :");
      put(p'last,1); new_line;
      for i in c'range loop
        put(c(i)'last,1); new_line;
        for j in c(i)'range loop
          put(c(i)(j));
          put("  ");
          put(e(i)(j).all);
          new_line;
        end loop;
      end loop;
    end;
  end Convert_Polynomial_System_to_Tableau_Structures;

  procedure Convert_Laurent_System_to_Tableau_Structures is

    p : Link_to_Laur_Sys;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(p);
    declare
      c : Standard_Complex_VecVecs.VecVec(p'range);
      e : Standard_Integer_VecVecs.Array_of_VecVecs(p'range);
    begin
      Extract_Coefficients_and_Exponents(p.all,c,e);
      new_line;
      put_line("The tableau structure :");
      put(p'last,1); new_line;
      for i in c'range loop
        put(c(i)'last,1); new_line;
        for j in c(i)'range loop
          put(c(i)(j));
          put("  ");
          put(e(i)(j).all);
          new_line;
        end loop;
      end loop;
    end;
  end Convert_Laurent_System_to_Tableau_Structures;

  procedure Main_Interactive_Driver is

    ans : character;

  begin
    new_line;
    put_line("MENU to convert formats of polynomial systems :");
    put_line("  1. convert symbolic polynomial system into tableau format;");
    put_line("  2. convert symbolic Laurent system into tableau format;");
    put_line("  3. convert tableau system into symbolic format;");
    put_line("  4. convert polynomial system with data structures;");
    put_line("  5. convert Laurent system with data structures;");
    put("Type 1, 2, 3, 4, or 5 to select the desired conversion : ");
    Ask_Alternative(ans,"12345");
    if ans = '1' then
      Convert_Polynomial_into_Tableau_Format;
    elsif ans = '2' then
      Convert_Laurent_into_Tableau_Format;
    elsif ans = '3' then
      Convert_Tableau_into_Symbolic_Format;
    elsif ans = '4' then
      Convert_Polynomial_System_to_Tableau_Structures;
    else
      Convert_Laurent_System_to_Tableau_Structures;
    end if;
  end Main_Interactive_Driver;

end Standard_Tableau_Formats;
