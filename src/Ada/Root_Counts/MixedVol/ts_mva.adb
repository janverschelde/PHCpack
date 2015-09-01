with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Integer_VecVecs;
with Lists_of_Floating_Vectors;         use Lists_of_Floating_Vectors;
--with Lists_of_Floating_Vectors_io;      use Lists_of_Floating_Vectors_io;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;    use Floating_Mixed_Subdivisions_io;
with Cell_Stack;                        use Cell_Stack;
with MixedVol_Algorithm;                use MixedVol_Algorithm;

procedure ts_mva is

-- DESCRIPTION :
--   Standalone routine to call the MixedVol algorithm.

  function Is_In ( f : List; v : Standard_Floating_Vectors.Vector;
                   tol : in double_float ) return boolean is

    tmp : List := f;
    w : Standard_Floating_Vectors.Link_to_Vector;
    eql : boolean := false;

  begin
    while not Is_Null(tmp) loop
      w := Head_Of(tmp);
      eql := true;
      for i in v'range loop
        if abs(v(i) - w(i)) > tol
         then eql := false; exit;
        end if;
      end loop;
      exit when eql;
      tmp := Tail_Of(tmp);
    end loop;
    return eql;
  end Is_In;

  procedure Swap ( v : in out Standard_Floating_Vectors.Link_to_Vector ) is

    t : constant double_float := v(v'first);

  begin
    v(v'first) := v(v'last);
    v(v'last) := t;
  end Swap;

  procedure Swap ( t : in out List ) is

    tmp : List := t;
    v : Standard_Floating_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      v := Head_Of(tmp);
      Swap(v);
      Set_Head(tmp,v);
      tmp := Tail_Of(tmp);
    end loop;
  end Swap;

  procedure Show_Candidate_Tropisms 
              ( sub : in Mixed_Subdivision; tol : in double_float ) is

  -- DESCRIPTION :
  --   Prints all inner normals with first component less than tol
  --   in absolute value.

    res,res_last : List;
    tmp : Mixed_Subdivision := sub;
    mic : Mixed_Cell;
    v : Standard_Floating_Vectors.Link_to_Vector;
    cnt : integer32 := 0;

  begin
    while not Is_Null(tmp) loop
      v := Head_Of(tmp).nor;
      if abs(v(v'first)) < tol then
        put(v,3); new_line;
        cnt := cnt + 1;      
        if not Is_In(res,v.all,tol)
         then Append(res,res_last,v.all);
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    Swap(res);
    put("computed "); put(cnt,1); put_line(" candidate tropisms");
    put("#different candidates : "); put(Length_Of(res),1);
    put_line(" list of candidate tropisms : "); put(res);
  end Show_Candidate_Tropisms;

  procedure Mixed_Volume_Computation
              ( file : in file_type; p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system and computes
  --   its mixed volume.  A mixed subdivision is written to file.

    n : constant integer32 := p'last;
    stlb : constant double_float := 0.0;
    ans : character;
    m,r,size,nb : integer32;
    mixvol : natural32;
    cnt,ind : Standard_Integer_Vectors.Vector(1..n);
    sup,mtype,perm,idx : Standard_Integer_Vectors.Link_to_Vector;
    vtx : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;
    cells : CellStack;
    sub : Mixed_Subdivision;
    ulift : boolean;
    cellcnt : natural32 := 0;

    procedure write_labels ( idx : Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      cellcnt := cellcnt + 1;
      put("Cell "); put(cellcnt,1); put(" is spanned by ");
      put(idx.all); new_line;
    end write_labels;

  begin
    Extract_Supports(n,p,m,ind,cnt,sup);
    put("Special lifting ? (y/n) ");
    Ask_Yes_or_No(ans);
    ulift := (ans = 'y');
    if ulift then
      uliftmv(n,m,ind,cnt,sup.all,stlb,r,mtype,perm,idx,vtx,lft,
              size,nb,cells,mixvol);
    else
      new_line;
      put("Test callback function ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        mv_with_callback(n,m,ind,cnt,sup.all,stlb,r,mtype,perm,idx,vtx,lft,
                         size,nb,cells,mixvol,false,write_labels'access);
      else
        mv(n,m,ind,cnt,sup.all,stlb,r,mtype,perm,idx,vtx,lft,
           size,nb,cells,mixvol);
      end if;
    end if;
    put("The mixed volume is "); put(mixvol,1); put_line(".");
    put("There are "); put(nb,1); put_line(" mixed cells.");
    put("Permutation of the supports : "); put(perm); new_line;
    put_line("Creating a regular mixed-cell configuration ...");
    if r < n
     then Create_Mixed_Cell_Configuration
            (n,r,size,nb,mtype,Vtx,lft,cells,sub);
     else Create_Mixed_Cell_Configuration
            (n,r,size,nb,mtype,perm,Vtx,lft,cells,sub);
    end if;
    if ulift then
      Show_Candidate_Tropisms(sub,1.0E-12);
    end if;
    new_line;
    put_line("See the output file for the mixed subdivision...");
    new_line;
    declare
      mix : Standard_Integer_Vectors.Vector(1..r);
    begin
      for i in 1..r loop
        mix(i) := mtype(i-1);
      end loop;
      put(file,natural32(n),mix,sub);
    end;
  end Mixed_Volume_Computation;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume.

    lp : Link_to_Laur_Sys;
    file : file_type;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    Mixed_Volume_Computation(file,lp.all);
  end Main;

begin
  Main;
end ts_mva;
