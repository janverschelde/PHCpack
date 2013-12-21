with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Matrix_Homotopies;

package body Matrix_Homotopies_io is

  procedure Write is
  begin
    Write(Standard_Output);
  end Write;

  procedure Write ( file : in file_type ) is
  begin
    for i in 1..Matrix_Homotopies.Cardinality loop
      declare
	start : constant Matrix := Matrix_Homotopies.Eval(i,Create(0.0));
	target : constant Matrix := Matrix_Homotopies.Eval(i,Create(1.0));
      begin
        put(file,"Matrix homotopy no. "); put(file,i,1);
        put_line(file," :");
        put_line(file,"Start matrix : "); put(file,start);
        put_line(file,"Target matrix : "); put(file,target);
      end;
    end loop;
  end Write;

end Matrix_Homotopies_io;
