with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Prompt_for_Systems;
with Witness_Sets_io;                    use Witness_Sets_io;
with Continuation_Parameters;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Greeting_Banners;
with Write_Seed_Number;
with Cascade_Homotopies;                 use Cascade_Homotopies;
with Embeddings_and_Cascades;            use Embeddings_and_Cascades;

package body Drivers_to_Cascade_Filtering is

  procedure Standard_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Standard_Laur_Poly_Convertors;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    filename : Link_to_String;
    tolzero : constant double_float := 1.0E-12;
    tolsing : constant double_float := 1.0E-8;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      Standard_Read_Embedding(infile,lq,sols,dim);
      filename := new string'(inpname);
    end if;
    if lq = null then
       new_line;
       declare
         name : constant string := Read_String;
       begin
         Open_Input_File(infile,name);
         Standard_Read_Embedding(infile,lq,sols,dim);
         filename := new string'(name);
       end;
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    Continuation_Parameters.Tune(0);
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    put_line("See the input and output file for results ...");
    new_line;
    close(infile);
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Witness_Generate
        (filename.all,outfile,nt,lq.all,sols,dim,0,tolzero,tolsing);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate
        (filename.all,outfile,nt,lp.all,sols,dim,0,tolzero,tolsing);
    end if;
  end Standard_Witness_Generate;

  procedure DoblDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Laur_Poly_Convertors;
    use DoblDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    filename : Link_to_String;
    tolzero : constant double_float := 1.0E-12;
    tolsing : constant double_float := 1.0E-8;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      DoblDobl_Read_Embedding(infile,lq,sols,dim);
      filename := new string'(inpname);
    end if;
    if lq = null then
       new_line;
       declare
         name : constant string := Read_String;
       begin
         Open_Input_File(infile,name);
         DoblDobl_Read_Embedding(infile,lq,sols,dim);
         filename := new string'(name);
       end;
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    Continuation_Parameters.Tune(0);
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    put_line("See the input and output file for results ...");
    new_line;
    close(infile);
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Witness_Generate
        (filename.all,outfile,nt,lq.all,sols,dim,0,tolzero,tolsing);
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate
        (filename.all,outfile,nt,lp.all,sols,dim,0,tolzero,tolsing);
    end if;
  end DoblDobl_Witness_Generate;

  procedure QuadDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Laur_Poly_Convertors;
    use QuadDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    filename : Link_to_String;
    tolzero : constant double_float := 1.0E-12;
    tolsing : constant double_float := 1.0E-8;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      QuadDobl_Read_Embedding(infile,lq,sols,dim);
      filename := new string'(inpname);
    end if;
    if lq = null then
       new_line;
       declare
         name : constant string := Read_String;
       begin
         Open_Input_File(infile,name);
         QuadDobl_Read_Embedding(infile,lq,sols,dim);
         filename := new string'(name);
       end;
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    Continuation_Parameters.Tune(0);
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    put_line("See the input and output file for results ...");
    new_line;
    close(infile);
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Witness_Generate
        (filename.all,outfile,nt,lq.all,sols,dim,0,tolzero,tolsing);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate
        (filename.all,outfile,nt,lp.all,sols,dim,0,tolzero,tolsing);
    end if;
  end QuadDobl_Witness_Generate;

  procedure Driver_to_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Witness_Generate(nt,inpname,outname);
      when '1' => DoblDobl_Witness_Generate(nt,inpname,outname);
      when '2' => QuadDobl_Witness_Generate(nt,inpname,outname);
      when others => null;
    end case;
  end Driver_to_Witness_Generate;

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string;
                verbose : in integer32 := 0 ) is

    use Standard_Laur_Poly_Convertors;

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    outfilename : Link_to_String;

  begin
    if verbose > 0 then
      put("-> in drivers_to_cascade_filtering.");
      put_line("Standard_Embed_and_Cascade ...");
    end if;
    Prompt_for_Systems.Read_System(infile,inpname,lq,sysonfile);
    Create_Output_File(outfile,outname,outfilename);
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Standard_Embed_and_Cascade
        (outfile,outfilename.all,nt,lq.all,true,true,verbose-1);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Standard_Embed_and_Cascade
        (outfile,outfilename.all,nt,lp.all,true,true,verbose-1);
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Standard_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string;
                verbose : in integer32 := 0 ) is

    use DoblDobl_Laur_Poly_Convertors;

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    outfilename : Link_to_String;

  begin
    if verbose > 0 then
      put("-> in drivers_to_cascade_filtering.");
      put_line("DoblDobl_Embed_and_Cascade ...");
    end if;
    Prompt_for_Systems.Read_System(infile,inpname,lq,sysonfile);
    Create_Output_File(outfile,outname,outfilename);
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      DoblDobl_Embed_and_Cascade
        (outfile,outfilename.all,nt,lq.all,true,true,verbose-1);
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      DoblDobl_Embed_and_Cascade
        (outfile,outfilename.all,nt,lp.all,true,true,verbose-1);
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end DoblDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Laur_Poly_Convertors;

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    outfilename : Link_to_String;

  begin
    if verbose > 0 then
      put("-> in drivers_to_cascade_filtering.");
      put_line("QuadDobl_Embed_and_Cascade ...");
    end if;
    Prompt_for_Systems.Read_System(infile,inpname,lq,sysonfile);
    Create_Output_File(outfile,outname,outfilename);
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      QuadDobl_Embed_and_Cascade
        (outfile,outfilename.all,nt,lq.all,true,true,verbose-1);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      QuadDobl_Embed_and_Cascade
        (outfile,outfilename.all,nt,lp.all,true,true,verbose-1);
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end QuadDobl_Embed_and_Cascade;

  procedure Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string;
                verbose : in integer32 := 0 ) is

    prc : constant character := Prompt_for_Precision;

  begin
    if verbose > 0 then
      put_line("-> in drivers_to_cascade_filtering.Embed_and_Cascade ...");
    end if;
    case prc is
      when '0' => Standard_Embed_and_Cascade(nt,inpname,outname,verbose-1);
      when '1' => DoblDobl_Embed_and_Cascade(nt,inpname,outname,verbose-1);
      when '2' => QuadDobl_Embed_and_Cascade(nt,inpname,outname,verbose-1);
      when others => null;
    end case;
  end Embed_and_Cascade;

end Drivers_to_Cascade_Filtering;
