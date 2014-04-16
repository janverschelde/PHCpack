with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Random_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
--with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Symbol_Table;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Planes_and_Polynomials;             use Planes_and_Polynomials;
with Standard_Plane_Representations;     use Standard_Plane_Representations;
with Standard_Intrinsic_Solutions;       use Standard_Intrinsic_Solutions;

-- for testing Recentered version:
--with Standard_Complex_Vectors_io;  use Standard_Complex_Vectors_io;
--with Standard_Complex_Poly_SysFun; use Standard_Complex_Poly_SysFun;

package body Intrinsic_Witness_Sets_io is

-- AUXILIARY PROCEDURES :

  function Append_wd ( s : string; d : natural32 ) return string is

  -- DESCRIPTION :
  --   Returns s_wd, where s is a file name and d a number.

    strd : constant string := convert(integer32(d));

  begin
    return s & "_w" & strd;
  end Append_wd;

  function Append_sd ( s : string; d : natural32 ) return string is

  -- DESCRIPTION :
  --   Returns s_sd, where s is a file name and d a number.

    strd : constant string := convert(integer32(d));

  begin
    return s & "_s" & strd;
  end Append_sd;

  function Embed_System ( f : in Poly_Sys; n,d : in natural32;
                          slices : in Matrix ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns an embedding of a polynomial system, with n original
  --   variables, to encode a solution set of dimension d.
  --   The matrix slices is the coefficient matrix of the linear equations,
  --   to be added to the system f.

    nd : constant integer32 := integer32(n+d);
    res : Poly_Sys(1..nd);
    cff : Vector(0..nd);

  begin
    for i in f'range loop                   -- embed polynomial system
      res(i) := Add_Embedding(f(i),d);
    end loop;
    if f'last < integer32(n) then           -- add dummy equations
      declare
        t : Term;
        ind : integer32 := integer32(n);
      begin
        t.cf := Create(1.0);
        t.dg := new Standard_Natural_Vectors.Vector'(1..nd => 0);
        for i in f'last+1..integer32(n) loop
          ind := ind + 1;
          t.dg(ind) := 1;
          res(i) := Create(t);
          t.dg(ind) := 0;
        end loop;
      end;
    end if;
    for i in 1..integer32(d) loop           -- add the slices
      for j in 0..integer32(n) loop
        cff(j) := slices(i,j);
      end loop;
      for j in 1..integer32(d) loop
        cff(integer32(n)+j) := Create(0.0); 
      end loop;
      cff(integer32(n)+i) := Create(1.0);
      res(integer32(n)+i) := Hyperplane(cff);
    end loop;
    return res;
  end Embed_System;

  function Count_Dummies ( p : Poly_Sys; n,d : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of dummy equations in the system p.
  --   A dummy equation is an equation which consists of one term,
  --   which is just one embedding variable zz = 0.

  -- ON ENTRY :
  --   p        the polynomial system, in n unknowns;
  --   n        number of unknowns in the system;
  --   d        number of embedding variables.

  -- REQUIRED : the embedding variables are the last d ones.

    res : natural32 := 0;

  begin
    for i in p'range loop
      if Number_of_Terms(p(i)) = 1 and then Degree(p(i)) = 1 then
        for j in n-d..n loop
          if Degree(p(i),integer32(j)) = 1
           then res := res + 1;
          end if;
        end loop;
      end if;
    end loop;
    return res;
  end Count_Dummies;

  function Intrinsic_Slices ( p : Poly_Sys; n,d : natural32 ) return Matrix is

  -- DESCRIPTION :
  --   Returns the intrinsic representation of the last d linear equations
  --   of the polynomial system in n unknowns.  The basis is orthogonal.

    sli : VecVec(1..integer32(d)) := Slices(p,d);
    equ : constant Matrix(1..integer32(d),0..integer32(n))
        := Equations_to_Matrix(sli,integer32(n));
    gen : constant Matrix(1..integer32(n),0..integer32(n-d))
        := Generators(equ);
    res : constant Matrix(1..integer32(n),0..integer32(n-d))
        := Orthogonalize(gen);

  begin
   -- put_line("The equations of the slices : "); put(equ,3);
   -- put_line("The generators : "); put(gen,3);
    Clear(sli);
    return res;
  end Intrinsic_Slices;

  procedure Write_Hyperplanes ( file : in file_type; c : in Matrix ) is

  -- DESCRIPTION :
  --   Converts the hyperplanes whose coefficients are given in c
  --   to a polynomial system and writes this system to file.

    hyp : VecVec(c'range(1)) := Equations_to_VecVec(c);
    sys : Poly_Sys(hyp'range);

  begin
    for i in hyp'range loop
      sys(i) := Hyperplane(hyp(i).all);
    end loop;
    put_line(file,sys);
    Clear(hyp);
    Clear(sys);
  end Write_Hyperplanes;

-- INPUT/OUTPUT of WITNESS STONES :

  procedure Read_Witness_Stone
              ( d,k : out natural32;
                s : out Solution_List; p : out Link_to_Matrix ) is

    infile : file_type;
    lp : Link_to_Poly_Sys;
    nv,ne : natural32 := 0;
    found : boolean;
    esols : Solution_List;

  begin
    new_line;
    put_line("Reading the file name for the witness stone ...");
    Read_Name_and_Open_File(infile);
    get(infile,ne);
    get(infile,nv);
    d := ne;
    Symbol_Table.Clear;
    Symbol_Table.Init(nv);
    if ne > 0 then
      lp := new Poly_Sys(1..integer32(ne));
      get(infile,lp.all);
      -- put_line("The hyperplanes : "); put_line(lp.all);
      p := new Matrix'(Intrinsic_Slices(lp.all,nv,ne));
    end if;
    Scan_and_Skip(infile,"SOLUTIONS",found);
    if found then
      get(infile,esols);
      -- put_line("The solutions : ");
      -- put(standard_output,Head_Of(esols).n,Length_Of(esols),esols);
      s := Project(esols,p.all);
      -- put_line("The solutions in intrinsic coordinates : ");
      -- put(standard_output,Head_Of(s).n,Length_Of(s),s);
      Clear(esols);
    end if;
    put("Read a witness set of dimension "); put(d,1);
    put(" in "); put(nv,1); put_line("-space.");
    k := 0;
    put("How many equations define this set ? "); get(k);
  end Read_Witness_Stone;

  procedure Read_Witness_Stone
              ( f : out Link_to_Poly_Sys; d,k : out natural32;
                s : out Solution_List; p : out Link_to_Matrix ) is

    esols : Solution_List;
   -- es2 : Solution_List;

  begin
    Standard_Read_Embedding(f,esols,d);
    k := natural32(f'last) - Count_Dummies(f.all,natural32(f'last),d) - d;
    put("Read a witness set of dimension "); put(d,1);
    put(" defined by "); put(k,1); put_line(" equations.");
    p := new Matrix'(Intrinsic_Slices(f.all,natural32(f'last)-d,d));
   -- put_line("The matrix of the plane in intrinsic coordinates : ");
   -- put(p.all,3);
    s := Project(esols,p.all);
   -- put_line("The solutions in intrinsic coordinates : ");
   -- put(standard_output,Head_Of(s).n,Length_Of(s),s);
   -- es2 := Expand(s,p.all);
   -- put_line("The solutions in extrinsic coordinates : ");
   -- put(standard_output,Head_Of(es2).n,Length_Of(es2),es2);
    Clear(esols);
  end Read_Witness_Stone;

  procedure Write_Witness_Stone
              ( file : in file_type; filename : in string; nv,d : in natural32;
                s : in Solution_List; p : in Matrix ) is

    witname : constant string := Append_sd(filename,d);
    witfile : file_type;
    esols : Solution_List := Expand(s,p);
    slices : constant Matrix(1..integer32(d),0..integer32(nv)) := Equations(p);

  begin
    put(file,"See the file " & witname & " for witness stone of dimension ");
    put(file,d,1); put_line(file,".");
    Create(witfile,out_file,witname);
    if d > 0
     then Write_Hyperplanes(witfile,slices);
     else put(witfile," 0 "); put(witfile,nv,1); new_line(witfile);
    end if;
    new_line(witfile);
    put(witfile,"TITLE : witness stone of dimension "); 
    put(witfile,d,1);
    put_line(witfile,", see " & filename & " for diagnostics."); 
    new_line(witfile);
    put_line(witfile,"THE SOLUTIONS :");
    new_line(witfile);
    put(witfile,Length_Of(esols),natural32(Head_Of(esols).n),esols);
    close(witfile);
    Clear(esols);
  end Write_Witness_Stone;

  procedure Write_Witness_Stone
              ( file : in file_type; filename : in string; nv,d : in natural32;
                f : in Poly_Sys; s : in Solution_List; p : in Matrix ) is

    witname : constant string := Append_sd(filename,d);
    witfile : file_type;
    wsols : Solution_List := Expand(s,p);
    esols : Solution_List := Add_Embedding(wsols,d);
    slices : constant Matrix(1..integer32(d),0..integer32(nv)) := Equations(p);
    sys : constant Poly_Sys(1..integer32(nv+d)) := Embed_System(f,nv,d,slices);

  begin
    put(file,"See the file " & witname & " for witness stone of dimension ");
    put(file,d,1); put_line(file,".");
    Create(witfile,out_file,witname);
    if Symbol_Table.Number < nv+d
     then Add_Embed_Symbols(d);
    end if;
    put(witfile,natural32(sys'last),sys);
    new_line(witfile);
    put(witfile,"TITLE : witness stone of dimension "); 
    put(witfile,d,1);
    put_line(witfile,", see " & filename & " for diagnostics."); 
    new_line(witfile);
    put_line(witfile,"THE SOLUTIONS :");
    new_line(witfile);
    put(witfile,Length_Of(esols),natural32(Head_Of(esols).n),esols);
    close(witfile);
    Clear(wsols); Clear(esols);
  end Write_Witness_Stone;

-- OUTPUT OF ONE WITNESS SET :

  procedure Write_Witness_Set
              ( file : in file_type; filename : in string; nv,d : in natural32;
                s : in Solution_List; p : in Matrix ) is

    witfile : file_type;
    witname : constant string := Append_wd(filename,d);
    esols : Solution_List := Expand(s,p);

  begin
    put(file,"See the file " & witname & " for solutions of dimension ");
    put(file,d,1); put(" in "); put(nv,1); put_line(file,"-space.");
    create(witfile,out_file,witname);
    put(witfile,Length_Of(esols),natural32(Head_Of(esols).n),esols);
    close(witfile);
    Clear(esols);
  end Write_Witness_Set;

  procedure Write_Witness_Set
              ( file : in file_type; nv,d : in natural32;
                f : in Poly_Sys; s : in Solution_List; p : in Matrix ) is

    wsols : Solution_List := Expand(s,p);
    esols : Solution_List := Add_Embedding(wsols,d);
    slices : constant Matrix(1..integer32(d),0..integer32(nv)) := Equations(p);
    sys : constant Poly_Sys(1..integer32(nv+d)) := Embed_System(f,nv,d,slices);

  begin
    if Symbol_Table.Number < nv+d
     then Add_Embed_Symbols(d);
    end if;
    put(file,natural32(sys'last),sys);
    new_line(file);
    put(file,"TITLE : witness set of dimension "); 
    put(file,d,1); new_line(file);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    new_line(file);
    put(file,Length_Of(esols),natural32(Head_Of(esols).n),esols);
    Clear(wsols); Clear(esols);
  end Write_Witness_Set;

  procedure Write_Recentered_Witness_Set
              ( file : in file_type; nv,d : in natural32;
                f : in Poly_Sys; s : in Solution_List; p : in Matrix ) is

    slices : constant Matrix(1..integer32(d),0..integer32(nv)) := Equations(p);
    sys : constant Poly_Sys(1..integer32(nv+d)) := Embed_System(f,nv,d,slices);

  begin
   -- put_line("Writing recentered witness sets ...");
   -- put("f'last = "); put(f'last,1); 
   -- put("  d = "); put(d,1); put("  nv = "); put(nv,1); new_line;
    if Symbol_Table.Number < nv+d
     then Add_Embed_Symbols(d);
    end if;
    put(file,natural32(sys'last),sys);
    new_line(file);
    put(file,"TITLE : witness set of dimension "); 
    put(file,d,1); new_line(file);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    new_line(file);
    put(file,Length_Of(s),natural32(Head_Of(s).n),s);
   -- declare
   --   tmp : Solution_List := s;
   -- begin
   --   for i in 1..Length_Of(s) loop
   --     put("evaluation of solution "); put(i,1); put_line(" : ");
   --     declare
   --       y : constant Vector := Eval(sys,Head_Of(tmp).v);
   --     begin
   --       put_line(y);
   --     end;
   --     tmp := Tail_Of(tmp);
   --   end loop;
   -- end;
  end Write_Recentered_Witness_Set;

  procedure Write_Witness_Set_to_File
              ( filename : in string; n,k : in natural32; f : in Poly_Sys;
                p : in Matrix; sols : in Solu_Info_Array ) is

    file : file_type;
    d : constant natural32 := n-k;
    strd : constant string := Characters_and_Numbers.convert(integer32(d));
    name : constant string := filename & "_w" & strd;
    s : constant Solution_List := Shallow_Create(sols);

  begin
    new_line;
    put("Writing the witness set to file "); put(name); put_line(" ...");
    put("See the file "); put(filename); put_line(" for diagnostics.");
    Create_Output_File(file,name);
    Write_Recentered_Witness_Set(file,n,d,f,s,p);
  end Write_Witness_Set_to_File;

  procedure Write_Witness_Set
              ( file : in file_type; filename : in string; nv,d : in natural32;
                f : in Poly_Sys; s : in Solution_List; p : in Matrix ) is

    witname : constant string := Append_wd(filename,d);
    witfile : file_type;
    wsols : Solution_List := Expand(s,p);
    esols : Solution_List := Add_Embedding(wsols,d);
    slices : constant Matrix(1..integer32(d),0..integer32(nv)) := Equations(p);
    sys : constant Poly_Sys(1..integer32(nv+d)) := Embed_System(f,nv,d,slices);

  begin
    declare
    begin
      Create(witfile,out_file,witname);
      put(file,"See the file " & witname & " for solutions of dimension ");
      put(file,d,1); put_line(file,".");
    exception
      when others =>
         put_line("Could not create the file " & witname & "...");
         declare
           pin : integer32 := Standard_Random_Numbers.Random(1000,9999);
           new_witname : constant string := Append_wd(witname,natural32(pin));
         begin
           Create(witfile,out_file,new_witname);
           put(file,"See the file " & new_witname
                  & " for solutions of dimension ");
           put(file,d,1); put_line(file,".");
           put_line("...created the file " & new_witname & " instead.");
         end;
    end;
    if Symbol_Table.Number < nv+d
     then Add_Embed_Symbols(d);
    end if;
    put(witfile,natural32(sys'last),sys);
    new_line(witfile);
    put(witfile,"TITLE : witness set of dimension "); 
    put(witfile,d,1);
    put_line(witfile,", see " & filename & " for diagnostics."); 
    new_line(witfile);
    put_line(witfile,"THE SOLUTIONS :");
    new_line(witfile);
    put(witfile,Length_Of(esols),natural32(Head_Of(esols).n),esols);
    close(witfile);
    Clear(wsols); Clear(esols);
  end Write_Witness_Set;

-- OUTPUT OF AN ARRAY OF WITNESS SETS :

  procedure Write_Witness_Sets
              ( file : in file_type; filename : in string; nv : in natural32;
                witset : in Array_of_Solution_Lists; planes : in VecMat ) is

    d : natural32;
    empty_sets : boolean := true;

  begin
    for i in witset'range loop
      d := nv - natural32(i);
      if not Is_Null(witset(i)) then
        Write_Witness_Set(file,filename,nv,d,witset(i),planes(i).all);
        empty_sets := false;
     -- else
     --   put(file,"There is no witness set of dimension "); put(file,d,1);
     --   put_line(file,".");
      end if;
    end loop;
    if empty_sets
     then put_line(file,"No witness sets found ...");
    end if;
  end Write_Witness_Sets;

  procedure Write_Witness_Sets
              ( file : in file_type; filename : in string; nv : in natural32;
                f : in Poly_Sys;
                witset : in Array_of_Solution_Lists; planes : in VecMat ) is

    d : natural32;
    empty_sets : boolean := true;

  begin
    for i in witset'range loop
      d := nv - natural32(i);
      if not Is_Null(witset(i)) then
        Write_Witness_Set(file,filename,nv,d,f,witset(i),planes(i).all);
        empty_sets := false;
     -- else
     --   put(file,"There is no witness set of dimension "); put(file,d,1);
     --   put_line(file,".");
      end if;
    end loop;
    if empty_sets
     then put_line(file,"No witness sets found ...");
    end if;
  end Write_Witness_Sets;

end Intrinsic_Witness_Sets_io;
