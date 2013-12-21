with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Hypersurfaces_and_Filters;         use Hypersurfaces_and_Filters;
with Intrinsic_Witness_Sets_io;         use Intrinsic_Witness_Sets_io;

procedure mainhyp ( polyfile,logfile : in string ) is

  procedure Create_Witness_Set
              ( file : in file_type; name : in string;
                n : in natural32; p : in Poly ) is

    d : constant integer32 := Degree(p);
    f : Eval_Poly := Create(p);
    b,v : Vector(1..integer32(n));
    s : Solution_List;
    res : double_float;
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
  end Create_Witness_Set;

  procedure Main is

    infile,outfile : file_type;
    n : natural32 := 0;
    p : Poly;
    name : Link_to_String;

  begin
    if polyfile = "" then
      new_line;
      put_line("Reading the name of the file with the polynomial.");
      Read_Name_and_Open_File(infile,name);
    else
      Open_Input_File(infile,polyfile,name);
    end if;
    get(infile,n);   -- skip the 1
    get(infile,n,p); -- read the polynomial
   -- put("the polynomial : "); put(p); new_line;
    if logfile = "" then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,logfile);
    end if;
    Create_Witness_Set(outfile,name.all,n,p);
  end Main;

begin
  Main;
end mainhyp;
