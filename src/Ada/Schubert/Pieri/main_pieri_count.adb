with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Localization_Posets;               use Localization_Posets;

package body Main_Pieri_Count is

  procedure Pieri_Count ( file : in file_type; m,p,q : in natural32 ) is

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    nq : constant natural32 := m*p + q*(m+p);
    nb : natural32;
    level_poset : Array_of_Nodes(0..integer32(nq));

  begin
    Q_Top_Bottom_Create(lnkroot,root.bottom(integer32(p)),natural32(m+p));
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    nb := natural32(level_poset(level_poset'last).roco);
    put(file,"Pieri root count : "); put(file,nb,1); new_line(file);
    Clear(level_poset);
  end Pieri_Count;

  procedure Main ( infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    infile,outfile : file_type;
    m,p,q : natural32 := 0;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in main_pieri_count.Main ...");
    end if;
    Open_Input_File(infile,infilename);
    Create_Output_File(outfile,outfilename);
    get(infile,m); put(outfile,"m = "); put(outfile,m,1);
    get(infile,p); put(outfile,"  p = "); put(outfile,p,1);
    get(infile,q); put(outfile,"  q = "); put(outfile,q,1);
    new_line(outfile);
    Close(infile);
    Pieri_Count(outfile,m,p,q);
    Close(outfile);
  end Main;

end Main_Pieri_Count;
