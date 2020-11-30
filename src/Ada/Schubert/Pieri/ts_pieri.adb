with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Main_Pieri_Homotopies;
with Main_Quantum_Pieri;

procedure ts_pieri is

-- DESCRIPTION :
--   Simply calls the main Pieri homotopy procedures.

  procedure Main is

    n,d,q : natural32 := 0;
    file : file_type;

  begin
    new_line;
    put_line("Pieri Homotopies for linear subspace intersections.");
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give n = m+p, dimension of working space : "); get(n);
    put("Give p, the dimension of the solution planes : "); get(d);
    put("Give q, degree of map (q=0 : simple Pieri) : "); get(q);
    put(file,"Pieri Homotopies for n = "); put(file,n,1);
    put(file,"  p = "); put(file,d,1);
    put(file,"  and q = "); put(file,q,1); new_line(file);
    if q = 0
     then Main_Pieri_Homotopies.Main(file,n,d);
     else Main_Quantum_Pieri.Main(file,n,d,q);
    end if;
  end Main;

begin
  Main;
end ts_pieri;
