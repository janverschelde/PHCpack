with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with String_Splitters;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Jaco_Matrices;     use Standard_Complex_Jaco_Matrices;

package body Standard_Evaluator_Packages is

  function Vector_Symbol ( i : integer32 ) return Symbol is

    res : Symbol;
    s : constant String := "x(" & Convert(i) & ")";
    cnt : natural := res'first;

  begin
    for i in s'range loop
      res(cnt) := s(integer(i));
      cnt := cnt+1;
      exit when cnt > res'last;
    end loop;
    for i in cnt..res'last loop
      res(i) := ' ';
    end loop;
    return res;
  end Vector_Symbol;

  procedure Replace_Symbols is

  -- DESCRIPTION :
  --   Replaces all symbols in the symbol table with vector entries:
  --   x(1), x(2), up to x(n).

    n : constant natural32 := Symbol_Table.Number;

  begin
    Symbol_Table.Clear;
    Symbol_Table.Init(n);
    for i in 1..integer32(n) loop
      Symbol_Table.Add(Vector_Symbol(i));
    end loop;
  end Replace_Symbols;

  procedure Create_Inline_Term_Evaluator
              ( file : in file_type; t : in Term ) is

  -- DESCRIPTION :
  --   Writes the code to evaluate one term on file.

    cff : boolean := false;

  begin
    new_line(file);
    if t.cf = Create(1.0) then
      if Sum(t.dg) = 0
       then put(file,"      + Create(1.0)");
       else put(file,"      + ");
      end if;
    elsif t.cf = Create(-1.0) then
      if Sum(t.dg) = 0
       then put(file,"      - Create(1.0)");
       else put(file,"      - ");
      end if;
    else
      put(file,"     + Create(");
      put(file,REAL_PART(t.cf));
      put(file,",");
      put(file,IMAG_PART(t.cf));
      put(file,")");
      cff := true;
    end if;
    for i in t.dg'range loop
      if t.dg(i) /= 0
       then if cff
             then put(file,"*");
             else cff := true;
            end if;
            put(file,"x(" & Convert(i) & ")");
            if t.dg(i) > 1
             then put(file,"**");
                  put(file,t.dg(i),1);
            end if;
      end if;
    end loop;
  end Create_Inline_Term_Evaluator;

  procedure Create_Inline_Polynomial_Evaluator
                   ( file : in file_type; p : in Poly ) is

    procedure Visit_Term ( t : in Term; continue : out boolean ) is
    begin
      Create_Inline_Term_Evaluator(file,t);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    put_line(file,";");
  end Create_Inline_Polynomial_Evaluator;

  procedure Create_Inline_System_Evaluator
               ( file : in file_type; funname : in String; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Writes the body of a function for an evaluator for p on file.

  begin
    put_line(file,"  function " & funname
                                & " ( x : Vector ) return Vector is");
    new_line(file);
    put_line(file,"    y : Vector(x'range);");
    new_line(file);
    put_line(file,"  begin");
    for i in p'range loop
      put(file,"    y(" & Convert(i) & ") := ");
      Create_Inline_Polynomial_Evaluator(file,p(i));
    end loop;
    put_line(file,"    return y;");
    put_line(file,"  end " & funname & ";");
  end Create_Inline_System_Evaluator;

  procedure Create_Inline_Jacobian_Evaluator
               ( file : in file_type; funname : in String; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Writes the body of a function to evaluate the Jacobian matrix of
  --   p on file.

    jm : Jaco_Mat(p'range,p'range) := Create(p);

  begin
    put_line(file,"  function " & funname 
                                & " ( x : Vector ) return Matrix is");
    new_line(file);
    put_line(file,"    y : Matrix(x'range,x'range);");
    new_line(file);
    put_line(file,"  begin");
    for i in p'range loop
      for j in p'range loop
        put(file,"    y(" & Convert(i) & "," & Convert(j) & ") := ");
        Create_Inline_Polynomial_Evaluator(file,jm(i,j));
      end loop;
    end loop;
    put_line(file,"    return y;");
    put_line(file,"  end " & funname & ";");
    Clear(jm);
  end Create_Inline_Jacobian_Evaluator;

  function Read_Package_Name return String is

  -- DESCRIPTION :
  --   Reads the package name from standard input and returns the string.

  begin
    put_line("Reading the name of the package.");
    declare
      s : constant String := String_Splitters.Read_String;
    begin
      return s;
    end;
  end Read_Package_Name;

  procedure Create ( packname : in String; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Creates a package that allows to evaluate the polynomial system p.

    specfile,bodyfile : file_type;

  begin
    Replace_Symbols;
    Create(specfile,out_file,packname & ".ads");
    put_line(specfile,"with Standard_Complex_Vectors;           "
                     & "use Standard_Complex_Vectors;");
    put_line(specfile,"with Standard_Complex_Matrices;          "
                     & "use Standard_Complex_Matrices;");
    new_line(specfile);
    put_line(specfile,"package " & packname & " is");
    new_line(specfile);
    put_line(specfile,"  function Eval_Sys ( x : Vector ) return Vector;");
    put_line(specfile,"  function Eval_Jaco ( x : Vector ) return Matrix;");
    new_line(specfile);
    put_line(specfile,"end " & packname & ";");
    Close(specfile);
    Create(bodyfile,out_file,packname & ".adb");
    put_line(bodyfile,"with Standard_Floating_Numbers;          "
                     & "use Standard_Floating_Numbers;");
    put_line(bodyfile,"with Standard_Complex_Numbers;           "
                     & "use Standard_Complex_Numbers;");
    new_line(bodyfile);
    put_line(bodyfile,"package body " & packname & " is");
    new_line(bodyfile);
    Create_Inline_System_Evaluator(bodyfile,"Eval_Sys",p);
    new_line(bodyfile);
    Create_Inline_Jacobian_Evaluator(bodyfile,"Eval_Jaco",p);
    new_line(bodyfile);
    put_line(bodyfile,"end " & packname & ";");
    Close(bodyfile);
  end Create;

  procedure Create ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Creates a package that allows to evaluate the polynomial system p.

    packname : constant String := Read_Package_Name;

  begin
    Create(packname,p);
  end Create;

end Standard_Evaluator_Packages;
