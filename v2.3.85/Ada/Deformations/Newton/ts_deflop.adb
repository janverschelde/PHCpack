with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;     use Standard_Complex_Jaco_Matrices;
with Standard_Embed_Polynomials;         use Standard_Embed_Polynomials;

procedure ts_deflop is

-- DESCRIPTION :
--   Applies the deflation operator to a system given by user.

  procedure Add_Multiplier_Symbols ( n : in natural32 ) is

  -- DESCRIPTION :
  --   Adds n multiplier symbols of the form mu[k], for k in 1..n.

    sb : Symbol;

  begin
    Symbol_Table.Enlarge(n);
    sb := (sb'range => ' ');
    for k in 1..n loop
      declare
        sb : Symbol;
        kst : constant string := Convert(integer32(k));
      begin
        sb := (sb'range => ' ');
        sb(1) := 'm'; sb(2) := 'u'; sb(3) := '[';
        for i in kst'range loop
          sb(3+i) := kst(i);
        end loop;
        sb(3+kst'last+1) := ']';
        Symbol_Table.Add(sb);
      end;
    end loop;
  end Add_Multiplier_Symbols;

  function Hyperplane ( n : natural32 ) return Poly is

  -- DESCRIPTION :
  --   Returns the slicing hyperplane in the n multiplier variables.

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..2*integer32(n) => 0);
    t.cf := Random1;
    res := Create(t);
    for i in 1..integer32(n) loop
      t.cf := Random1;
      t.dg(integer32(n)+i) := 1;
      Add(res,t);
      t.dg(integer32(n)+i) := 0;
    end loop;
    Clear(t);
    return res;
  end Hyperplane;

  procedure Deflate ( file : in file_type; p : in Poly_Sys ) is

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    jm : Jaco_Mat(1..nq,1..nv) := Create(p);
    ejm : Jaco_Mat(1..nq,1..nv) := Add_Variables(jm,natural32(nv));
    hp : Poly := Hyperplane(natural32(nv));
    jmu : Poly;
    trm : Term;

  begin
    Add_Multiplier_Symbols(natural32(nv));
    put(file,nq+1,1); put(file," "); put(file,2*nv,1); new_line(file);
    put(file,hp); new_line(file);
    Clear(hp);
    trm.dg := new Standard_Natural_Vectors.Vector'(1..2*nv => 0);
    trm.cf := Create(1.0);
    for i in 1..nq loop
      jmu := Null_Poly;
      for j in 1..nv loop
        trm.dg(nv+j) := 1;
        Mul(ejm(i,j),trm);
        Add(jmu,ejm(i,j));
        trm.dg(nv+j) := 0;
      end loop;
      put(file,jmu); new_line(file);
      Clear(jmu);
    end loop;
    Clear(trm); Clear(jm);
  end Deflate;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and an output file.
  --   A corank one first-order deflation is written to the file.

    p : Link_to_Poly_Sys;
    file : file_type;

  begin
    new_line;
    put_line("Corank one first-order deflation operator...");
    new_line;
    get(p);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Deflate(file,p.all);
  end Main;

begin
  Main;
end ts_deflop;
