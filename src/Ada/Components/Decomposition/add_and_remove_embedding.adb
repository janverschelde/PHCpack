with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Prompt_for_Systems;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Square_and_Embed_Systems;           use Square_and_Embed_Systems;
with Greeting_Banners;
with Write_Seed_Number;

package body Add_and_Remove_Embedding is

  procedure Standard_Square_and_Embed ( iptname,optname : in string ) is

    use Standard_Laur_Poly_Convertors;

    lp,ep : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq,eq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    k : natural32;

  begin
    Prompt_for_Systems.Read_System(infile,iptname,lq,sysonfile);
    Create_Output_File(outfile,optname);
    new_line;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Interactive_Square_and_Embed(outfile,lq.all,eq,k);     
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Interactive_Square_and_Embed(outfile,lp.all,ep,k);
    end if;
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Standard_Square_and_Embed;

  procedure DoblDobl_Square_and_Embed ( iptname,optname : in string ) is

    use DoblDobl_Laur_Poly_Convertors;

    lp,ep : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq,eq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    k : natural32;

  begin
    Prompt_for_Systems.Read_System(infile,iptname,lq,sysonfile);
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Interactive_Square_and_Embed(outfile,lq.all,eq,k);     
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Create_Output_File(outfile,optname);
      new_line;
      Interactive_Square_and_Embed(outfile,lp.all,ep,k);
    end if;
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end DoblDobl_Square_and_Embed;

  procedure QuadDobl_Square_and_Embed ( iptname,optname : in string ) is

    use QuadDobl_Laur_Poly_Convertors;

    lp,ep : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq,eq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    k : natural32;

  begin
    Prompt_for_Systems.Read_System(infile,iptname,lq,sysonfile);
    Create_Output_File(outfile,optname);
    new_line;
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Interactive_Square_and_Embed(outfile,lq.all,eq,k);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Interactive_Square_and_Embed(outfile,lp.all,ep,k);
    end if;
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end QuadDobl_Square_and_Embed;

  procedure Driver_to_Square_and_Embed ( iptname,optname : in string ) is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Square_and_Embed(iptname,optname);
      when '1' => DoblDobl_Square_and_Embed(iptname,optname);
      when '2' => QuadDobl_Square_and_Embed(iptname,optname);
      when others => null;
    end case;
  end Driver_to_Square_and_Embed;

  procedure Standard_Remove_Embedding
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k,ns : in natural32 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    rp : constant Poly_Sys := Remove_Embedding(p,k,ns);
    nq : constant natural32 := natural32(rp'last);
    nv : constant natural32 := Number_of_Unknowns(rp(rp'first));
    rsols : Solution_List;

  begin
    if k + ns > 0
     then rsols := Remove_Embedding(sols,k+ns);
     else rsols := sols;
    end if;
    put(file,nq,1);
    put(file," ");
    put(file,nv,1);
    new_line(file);
    put(file,rp);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
  end Standard_Remove_Embedding;

  procedure Standard_Remove_Embedding
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k,ns : in natural32 ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    rp : constant Laur_Sys := Remove_Embedding(p,k,ns);
    nq : constant natural32 := natural32(rp'last);
    nv : constant natural32 := Number_of_Unknowns(rp(rp'first));
    rsols : Solution_List;

  begin
    if k + ns > 0
     then rsols := Remove_Embedding(sols,k+ns);
     else rsols := sols;
    end if;
    put(file,nq,1);
    put(file," ");
    put(file,nv,1);
    new_line(file);
    put(file,rp);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
  end Standard_Remove_Embedding;

  procedure DoblDobl_Remove_Embedding
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k,ns : in natural32 ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    rp : constant Poly_Sys := Remove_Embedding(p,k,ns);
    nq : constant natural32 := natural32(rp'last);
    nv : constant natural32 := Number_of_Unknowns(rp(rp'first));
    rsols : Solution_List;

  begin
    if k + ns > 0
     then rsols := Remove_Embedding(sols,k+ns);
     else rsols := sols;
    end if;
    put(file,nq,1);
    put(file," ");
    put(file,nv,1);
    new_line(file);
    put(file,rp);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
  end DoblDobl_Remove_Embedding;

  procedure DoblDobl_Remove_Embedding
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k,ns : in natural32 ) is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    rp : constant Laur_Sys := Remove_Embedding(p,k,ns);
    nq : constant natural32 := natural32(rp'last);
    nv : constant natural32 := Number_of_Unknowns(rp(rp'first));
    rsols : Solution_List;

  begin
    if k + ns > 0
     then rsols := Remove_Embedding(sols,k+ns);
     else rsols := sols;
    end if;
    put(file,nq,1);
    put(file," ");
    put(file,nv,1);
    new_line(file);
    put(file,rp);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
  end DoblDobl_Remove_Embedding;

  procedure QuadDobl_Remove_Embedding
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k,ns : in natural32 ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    rp : constant Poly_Sys := Remove_Embedding(p,k,ns);
    nq : constant natural32 := natural32(rp'last);
    nv : constant natural32 := Number_of_Unknowns(rp(rp'first));
    rsols : Solution_List;

  begin
    if k + ns > 0
     then rsols := Remove_Embedding(sols,k+ns);
     else rsols := sols;
    end if;
    put(file,nq,1);
    put(file," ");
    put(file,nv,1);
    new_line(file);
    put(file,rp);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
  end QuadDobl_Remove_Embedding;

  procedure QuadDobl_Remove_Embedding
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k,ns : in natural32 ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    rp : constant Laur_Sys := Remove_Embedding(p,k,ns);
    nq : constant natural32 := natural32(rp'last);
    nv : constant natural32 := Number_of_Unknowns(rp(rp'first));
    rsols : Solution_List;

  begin
    if k + ns > 0
     then rsols := Remove_Embedding(sols,k+ns);
     else rsols := sols;
    end if;
    put(file,nq,1);
    put(file," ");
    put(file,nv,1);
    new_line(file);
    put(file,rp);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(rsols),natural32(Head_Of(rsols).n),rsols);
  end QuadDobl_Remove_Embedding;

  procedure Standard_Remove_Embedding ( inpname,outname : in string ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Standard_Laur_Poly_Convertors;
    use Standard_Complex_Solutions;

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys := null;
    lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys := null;
    sols : Solution_List;
    infile,outfile : file_type;
    k,ns : natural32 := 0;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      Standard_Read_Embedding(infile,lq,sols,k,ns);
    end if;
    if lq = null then
      new_line;
      put_line("Removing an Embedding of a Polynomial System.");
      Standard_Read_Embedding(lq,sols,k,ns);
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    put("Number of embed variables in the system : ");
    put(k,1); new_line;
    put("Number of extra slack variables to remove : ");
    put(ns,1); new_line;
    new_line;
    put_line("See the output file for the system ...");
    new_line;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Standard_Remove_Embedding(outfile,lq.all,sols,k,ns);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Standard_Remove_Embedding(outfile,lp.all,sols,k,ns);
    end if;
    Close(outfile);
  end Standard_Remove_Embedding;

  procedure DoblDobl_Remove_Embedding ( inpname,outname : in string ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Laur_Poly_Convertors;
    use DoblDobl_Complex_Solutions;

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    infile,outfile : file_type;
    k,ns : natural32 := 0;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      DoblDobl_Read_Embedding(infile,lq,sols,k,ns);
    end if;
    if lq = null then
      new_line;
      put_line("Removing an Embedding of a Polynomial System.");
      DoblDobl_Read_Embedding(lq,sols,k,ns);
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    put("Number of embed variables in the system : ");
    put(k,1); new_line;
    put("Number of extra slack variables to remove : ");
    put(ns,1); new_line;
    new_line;
    put_line("See the output file for the system ...");
    new_line;
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      DoblDobl_Remove_Embedding(outfile,lq.all,sols,k,ns);
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      DoblDobl_Remove_Embedding(outfile,lp.all,sols,k,ns);
    end if;
    close(outfile);
  end DoblDobl_Remove_Embedding;

  procedure QuadDobl_Remove_Embedding ( inpname,outname : in string ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Laur_Poly_Convertors;
    use QuadDobl_Complex_Solutions;

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys := null;
    lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys := null;
    sols : Solution_List;
    infile,outfile : file_type;
    k,ns : natural32 := 0;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      QuadDobl_Read_Embedding(infile,lq,sols,k,ns);
    end if;
    if lq = null then
      new_line;
      put_line("Removing an Embedding of a Polynomial System.");
      QuadDobl_Read_Embedding(lq,sols,k,ns);
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    put("Number of embed variables in the system : ");
    put(k,1); new_line;
    put("Number of extra slack variables to remove : ");
    put(ns,1); new_line;
    new_line;
    put_line("See the output file for the system ...");
    new_line;
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      QuadDobl_Remove_Embedding(outfile,lq.all,sols,k,ns);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      QuadDobl_Remove_Embedding(outfile,lp.all,sols,k,ns);
    end if;
    close(outfile);
  end QuadDobl_Remove_Embedding;

  procedure Driver_to_Remove_Embedding ( inpname,outname : in string ) is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Remove_Embedding(inpname,outname);
      when '1' => DoblDobl_Remove_Embedding(inpname,outname);
      when '2' => QuadDobl_Remove_Embedding(inpname,outname);
      when others => null;
    end case;
  end Driver_to_Remove_Embedding;

end Add_and_Remove_Embedding;
