with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Path_Counts_Table;
with Embeddings_and_Cascades;            use Embeddings_and_Cascades;
with Standard_Witness_Solutions;
with DoblDobl_Witness_Solutions;
with QuadDobl_Witness_Solutions;
with Store_Witness_Solutions;            use Store_Witness_Solutions;
with Write_Witness_Solutions;            use Write_Witness_Solutions;

procedure ts_witsols is

-- DESCRIPTION :
--   Test on the package to store witness solutions.

  procedure Prompt_for_Options ( filter,factor : out boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the settings of the options
  --   filter and factor, returned on output.

    ans : character;

  begin
    new_line;
    put("Filter the witness supersets ? (y/n) ");
    Ask_Yes_or_No(ans);
    filter := (ans = 'y');
    if not filter then
      factor := false;
    else
      put("Factor the witness sets ? (y/n) ");
      Ask_Yes_or_No(ans);
      factor := (ans = 'y');
    end if;
  end Prompt_for_Options;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   and the prompts for the top dimension.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nt,nq,nv,topdim,lowdim : natural32 := 0;
    filter,factor : boolean;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    nq := natural32(lp'last);
    nv := Number_of_Unknowns(lp(lp'first));
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Prompt_for_Options(filter,factor);
    Standard_Witness_Solutions.Initialize(topdim);
    new_line;
    put_line("Computing ...");
    new_line;
    Standard_Solve_with_Callback
      (nt,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
    Standard_Write(topdim,lowdim);
    Write_Counts(filter,factor,pc,fc,idxfac);
  end Standard_Main;

  procedure Standard_Laurent_Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system,
  --   and the prompts for the top dimension.

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;

    lp : Link_to_Laur_Sys;
    nt,nq,nv,topdim,lowdim : natural32 := 0;
    filter,factor : boolean;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    nq := natural32(lp'last);
    nv := Number_of_Unknowns(lp(lp'first));
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Prompt_for_Options(filter,factor);
    new_line;
    put_line("Computing ...");
    new_line;
    Standard_Witness_Solutions.Initialize(topdim);
    Standard_Solve_with_Callback
      (nt,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
    Standard_Write(topdim,lowdim);
    Write_Counts(filter,factor,pc,fc,idxfac);
  end Standard_Laurent_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   and the prompts for the top dimension.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nt,nq,nv,topdim,lowdim : natural32 := 0;
    filter,factor : boolean;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    nq := natural32(lp'last);
    nv := Number_of_Unknowns(lp(lp'first));
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Prompt_for_Options(filter,factor);
    new_line;
    put_line("Computing ...");
    new_line;
    DoblDobl_Witness_Solutions.Initialize(topdim);
    DoblDobl_Solve_with_Callback
      (nt,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
    DoblDobl_Write(topdim,lowdim);
    Write_Counts(filter,factor,pc,fc,idxfac);
  end DoblDobl_Main;

  procedure DoblDobl_Laurent_Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system,
  --   and the prompts for the top dimension.

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;

    lp : Link_to_Laur_Sys;
    nt,nq,nv,topdim,lowdim : natural32 := 0;
    filter,factor : boolean;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    nq := natural32(lp'last);
    nv := Number_of_Unknowns(lp(lp'first));
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Prompt_for_Options(filter,factor);
    new_line;
    put_line("Computing ...");
    new_line;
    DoblDobl_Witness_Solutions.Initialize(topdim);
    DoblDobl_Solve_with_Callback
      (nt,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
    DoblDobl_Write(topdim,lowdim);
    Write_Counts(filter,factor,pc,fc,idxfac);
  end DoblDobl_Laurent_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   and the prompts for the top dimension.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nt,nq,nv,topdim,lowdim : natural32 := 0;
    filter,factor : boolean;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    nq := natural32(lp'last);
    nv := Number_of_Unknowns(lp(lp'first));
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Prompt_for_Options(filter,factor);
    new_line;
    put_line("Computing ...");
    new_line;
    QuadDobl_Witness_Solutions.Initialize(topdim);
    QuadDobl_Solve_with_Callback
      (nt,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
    QuadDobl_Write(topdim,lowdim);
    Write_Counts(filter,factor,pc,fc,idxfac);
  end QuadDobl_Main;

  procedure QuadDobl_Laurent_Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system,
  --   and the prompts for the top dimension.

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;

    lp : Link_to_Laur_Sys;
    nt,nq,nv,topdim,lowdim : natural32 := 0;
    filter,factor : boolean;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    nq := natural32(lp'last);
    nv := Number_of_Unknowns(lp(lp'first));
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Prompt_for_Options(filter,factor);
    new_line;
    put_line("Computing ...");
    new_line;
    QuadDobl_Witness_Solutions.Initialize(topdim);
    QuadDobl_Solve_with_Callback
      (nt,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
    QuadDobl_Write(topdim,lowdim);
    Write_Counts(filter,factor,pc,fc,idxfac);
  end QuadDobl_Laurent_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision
  --   and then calls the corresponding main test.

    prec,laur : character;

  begin
    new_line;
    put_line("MENU for the precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, 2 to select the precision : ");
    Ask_Alternative(prec,"012");
    new_line;
    put("Laurent polynomial system ? (y/n) ");
    Ask_Yes_or_No(laur);
    if laur = 'y' then
      case prec is
        when '0' => Standard_Laurent_Main;
        when '1' => DoblDobl_Laurent_Main;
        when '2' => QuadDobl_Laurent_Main;
        when others => null;
      end case;
    else
      case prec is
        when '0' => Standard_Main;
        when '1' => DoblDobl_Main;
        when '2' => QuadDobl_Main;
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_witsols;
