with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Embed_Polynomials;        use DoblDobl_Embed_Polynomials;
with DoblDobl_Complex_Solutions;        use DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Witness_Sets_Formats;              use Witness_Sets_Formats;
with DoblDobl_Plane_Representations;    use DoblDobl_Plane_Representations;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with DoblDobl_Hypersurface_Witsets;     use DoblDobl_Hypersurface_Witsets;

package body DoblDobl_Hypersurface_Witsets_io is

  procedure Write_Witness_Set
               ( n : in integer32; p : in Poly; b,v,t : in Vector ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Write_Witness_Set(file,n,p,b,v,t);
  end Write_Witness_Set;

  procedure Write_Witness_Set
               ( file : in file_type;
                 n : in integer32; p : in Poly; b,v,t : in Vector ) is

    dim : constant integer32 := 2*n-1;
    ep : constant Poly := Add_Embedding(p,natural32(n-1));
    sols : constant Solution_List := Embedded_Extrinsic_Solutions(n,b,v,t);
    eqs : constant DoblDobl_Complex_VecVecs.VecVec(1..n-1) := Equations1(b,v);
    hyp : Poly;
    sys : Poly_Sys(1..dim);
    et : Term;

  begin
    Witness_Sets_io.Add_Embed_Symbols(natural32(n-1));
    sys(1) := ep;
    et.cf := Create(integer(1));
    et.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    for i in eqs'range loop
      hyp := Hyperplane(eqs(i).all);
      sys(n+i) := Add_Variables(hyp,natural32(n-1));
      et.dg(n+i) := 1;
      Add(sys(n+i),et);
      sys(i+1) := Create(et);
      et.dg(n+i) := 0;
    end loop;
    Clear(et);
   -- put(file,natural32(dim),sys);
    put(file,sys);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(dim),sols);
  end Write_Witness_Set;

end DoblDobl_Hypersurface_Witsets_io;
