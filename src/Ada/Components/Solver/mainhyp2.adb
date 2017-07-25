with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Vectors;          use DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;           use DoblDobl_Random_Vectors;
with DoblDobl_Complex_Matrices;         use DoblDobl_Complex_Matrices;
with Symbol_Table;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laurentials_io;   use DoblDobl_Complex_Laurentials_io;
with DoblDobl_Laur_Poly_Convertors;
with DoblDobl_Complex_Poly_Functions;   use DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;        use DoblDobl_Complex_Solutions;
with Hypersurfaces_and_Filters;         use Hypersurfaces_and_Filters;
with Intrinsic_Witness_Sets_io;         use Intrinsic_Witness_Sets_io;

procedure mainhyp2 ( polyfile,logfile : in string ) is

  procedure Make_Witness_Set
              ( file : in file_type; name : in string;
                n : in natural32;
                p : in DoblDobl_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Generates a random offset and random direction vector
  --   and then computes all roots of a univariate problem
  --   to make a witness set for the hypersurface defined by p.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     name of the input file, used as prefix
  --            for the witness set;
  --   n        number of variables in the polynomial p;
  --   p        an ordinary polynomial p in n veriables.

    use DoblDobl_Complex_Polynomials;

    d : constant integer32 := Degree(p);
    f : Eval_Poly := Create(p);
    b,v : Vector(1..integer32(n));
    s : Solution_List;
    res : double_double;
    plane : Matrix(1..integer32(n),0..1);
    sys : Poly_Sys(1..1);

  begin
    b := Random_Vector(1,integer32(n));
    v := Random_Vector(1,integer32(n));
   -- put_line("Calling RP_Hypersurface_Witness_Set...");
    RP_Hypersurface_Witness_Set(file,n,natural32(d),f,b,v,s,res);
   -- put_line("...finished with RP_Hypersurface_Witness_Set.");
    Clear(f);   
    sys(1) := p;
    for i in b'range loop
      plane(i,0) := b(i);
      plane(i,1) := v(i);
    end loop;
    Write_Witness_Set(file,name,n,n-1,sys,s,plane);
  end Make_Witness_Set;

  procedure Read_Input_Polynomial
              ( name : out Link_to_String;
                n : out natural32;
                p : out DoblDobl_Complex_Laurentials.Poly ) is

  -- DESCRIPTION :
  --   Read the polynomial from polyfile, or 
  --   if polyfile is the empty string, then prompts the user
  --   for a file name to read the polynomial from.

  -- ON RETURN :
  --   name     name of the input file;
  --   n        number of variables in the polynomial;
  --   p        the polynomial read.

    use DoblDobl_Complex_Laurentials;

    infile : file_type;
    m : natural32 := 0;
    q : Poly;
    filename : Link_to_String;

  begin
    if polyfile = "" then
      new_line;
      put_line("Reading the name of the file with the polynomial.");
      Read_Name_and_Open_File(infile,filename);
    else
      Open_Input_File(infile,polyfile,filename);
    end if;
    get(infile,m);   -- skip the 1
    get(infile,m);   -- get the number of variables
    Symbol_Table.Init(m);
    get(infile,q);
    put("number of variables : "); put(m,1); new_line;
    put("the polynomial : "); put(q); new_line;
    n := m;
    p := q;
    name := filename;
  exception
    when others =>
      put_line("An error occurred with the polynomial on input.");
  end Read_Input_Polynomial;

  procedure Main is

  -- DESCRIPTION :
  --   Reads the input and checks whether the input polynomial
  --   is a genuine Laurent polynomial.

    infile,outfile : file_type;
    n : natural32 := 0;
    p : DoblDobl_Complex_Polynomials.Poly;
    q : DoblDobl_Complex_Laurentials.Poly;
    name : Link_to_String;

  begin
    Read_Input_Polynomial(name,n,q);
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(q) then
      put_line("Processing of Laurent polynomials not yet supported.");
      return;
    else
      p := DoblDobl_Laur_Poly_Convertors.Positive_Laurent_Polynomial(q);
    end if;
    put_line("The input polynomial p :"); put(p);
    if logfile = "" then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,logfile);
    end if;
    Make_Witness_Set(outfile,name.all,n,p);
  end Main;

begin
  Main;
end mainhyp2;
