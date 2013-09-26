with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Matrices;
with Symbol_Table;                       use Symbol_Table;
with Symbol_Table_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Plane_Representations;              use Plane_Representations;
with Minor_Computations;                 use Minor_Computations;

procedure ts_topos is

-- DESCRIPTION :
--   This procedure checks whether matrices on input are totally positive.

  function Convert ( cmpmat : Standard_Complex_Matrices.Matrix )
                   return Standard_Floating_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a matrix with the real parts of cmpmat.

    res : Standard_Floating_Matrices.Matrix(cmpmat'range(1),cmpmat'range(2));
 
  begin
    for i in cmpmat'range(1) loop
      for j in cmpmat'range(2) loop
        res(i,j) := REAL_PART(cmpmat(i,j));
      end loop;
    end loop;
    return res;
  end Convert;

  procedure Read_Matrices ( infile,outfile : in file_type;
                            n,m,k : in natural32 ) is

    mat : Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in 1..k loop
      get(infile,mat);
      put(outfile,"matrix "); put(outfile,i,1); put_line(outfile," :");
      put(outfile,mat);
      Minors(outfile,n,m,mat);
    end loop;
  end Read_Matrices;

  procedure Test_Matrices is

    infile,outfile : file_type;
    n,m,k : natural32 := 0;

  begin
    new_line;
    put_line("Reading the name of the file with the input planes.");
    Read_Name_and_Open_File(infile);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line;
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    put("Give number of matrices : "); get(k);
    Read_Matrices(infile,outfile,n,m,k);
  end Test_Matrices;

  function Read_Localization_Pattern 
             return Standard_Natural_Matrices.Matrix is

    n,m : natural32 := 0;

  begin
    new_line;
    put_line("Reading the localization pattern.");
    new_line;
    put("Give number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      res : Standard_Natural_Matrices.Matrix(1..integer32(n),1..integer32(m));
    begin
      put("Give an "); put(n,1); put("-by-"); put(m,1);
      put_line(" integer matrix to represent the localization : ");
      get(res); skip_line;
      return res;
    end;
  end Read_Localization_Pattern;

  procedure Initialize_Symbol_Table ( n : in natural32 ) is
  begin
    put("Give "); put(n,1); put_line(" symbols to initialize symbol table :");
    Symbol_Table.Init(n);
    for i in 1..n loop
      declare
        sb : Symbol;
      begin
        Symbol_Table_io.get(sb);
        Symbol_Table.Add(sb);
      end;
    end loop;
  end Initialize_Symbol_Table; 

  procedure Write_Signs ( file : in file_type;
                          sigs : in Standard_Integer_Vectors.Vector ) is
  begin
    put(file,"Signs : ");
    for i in sigs'range loop
      if sigs(i) = 0 then
        put(file,"0");
      elsif sigs(i) > 0 then
        put(file,"+");
      else
        put(file,"-");
      end if;
    end loop;
    new_line(file);
  end Write_Signs;

  procedure Test_Solutions is

    file : file_type;
    sols,tmp : Solution_List;
    cnt : natural32 := 0;

  begin
    new_line;
    put_line("Testing total positivity for solution planes.");
    new_line;
    Read(sols);
    declare
      locpat : constant Standard_Natural_Matrices.Matrix
             := Read_Localization_Pattern;
    begin
      Initialize_Symbol_Table(natural32(Head_Of(sols).n)); skip_line;
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
      tmp := sols;
      while not Is_Null(tmp) loop
        cnt := cnt + 1;
        declare
          ls : Link_to_Solution := Head_Of(tmp);
          cmpmatpla : constant Standard_Complex_Matrices.Matrix
                    := Matrix_Rep(locpat,ls.v);
          fltmatpla : constant Standard_Floating_Matrices.Matrix
                    := Convert(cmpmatpla);
          sigmins : constant Standard_Integer_Vectors.Vector
                  := Sign_of_Minors
                       (fltmatpla'length(1),fltmatpla'length(2),fltmatpla);
        begin
          put(file,"SOLUTION "); put(file,cnt,1); put_line(file," :");
          put(file,ls.all); new_line(file);
          put_line(file,"MATRIX REPRESENTATION : ");
          put(file,fltmatpla,3);
          Write_Signs(file,sigmins);
          Minors(file,fltmatpla'length(1),fltmatpla'length(2),fltmatpla);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end;
  end Test_Solutions;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing whether matrices are totally positive.");
    new_line;
    put_line("Choose one of the following : ");
    put_line("  1. Total positivity on a given set of matrices.");
    put_line("  2. Total positivity on a list of solution planes.");
    put("Type 1 or 2 to select : "); get(ans); skip_line;
    if ans = '1'
     then Test_Matrices;
     else Test_Solutions;
    end if;
  end Main;

begin
  Main;
end ts_topos;
