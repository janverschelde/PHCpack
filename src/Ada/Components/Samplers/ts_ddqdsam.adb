with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with DoblDobl_Complex_VecVecs;          use DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;          use QuadDobl_Complex_VecVecs;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Witness_Sets;
with Witness_Sets_io;                   use Witness_Sets_io;
with DoblDobl_Sampling_Machine;
with DoblDobl_Sample_Points;
with DoblDobl_Sample_Lists;
with DoblDobl_Sample_Grids;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sample_Points;
with QuadDobl_Sample_Lists;
with QuadDobl_Sample_Grids;

procedure ts_ddqdsam is

-- DESCRIPTION :
--   Interactive tests on the sampling machine in double double arithmetic.

  procedure Compute_DoblDobl_Samples
              ( sols : in DoblDobl_Complex_Solutions.Solution_List;
                n,d : integer32 ) is

  -- DESCRIPTION :
  --   Computes new samples on the solution set defined by a witness set,
  --   initialized in the sampling machine.

  -- REQUIRED :
  --   The double double sampling machine is initialized.

  -- ON INPUT :
  --   sols     solutions to start the sampling;
  --   n        ambient dimension of the embedded polynomial system;
  --   d        dimension of the solution set.

    hyp : DoblDobl_Complex_VecVecs.VecVec(1..d)
        := Witness_Sets.Random_Hyperplanes(natural32(d),natural32(n));
    newsol : DoblDobl_Complex_Solutions.Solution(n);

  begin
    DoblDobl_Sampling_Machine.Sample
      (standard_output,true,
       DoblDobl_Complex_Solutions.Head_Of(sols).all,hyp,newsol);
  end Compute_DoblDobl_Samples;

  procedure Compute_QuadDobl_Samples
              ( sols : in QuadDobl_Complex_Solutions.Solution_List;
                n,d : integer32 ) is

  -- DESCRIPTION :
  --   Computes new samples on the solution set defined by a witness set,
  --   initialized in the sampling machine.

  -- REQUIRED :
  --   The quad double sampling machine is initialized.

  -- ON INPUT :
  --   sols     solutions to start the sampling;
  --   n        ambient dimension of the embedded polynomial system;
  --   d        dimension of the solution set.

    hyp : QuadDobl_Complex_VecVecs.VecVec(1..d)
        := Witness_Sets.Random_Hyperplanes(natural32(d),natural32(n));
    newsol : QuadDobl_Complex_Solutions.Solution(n);

  begin
    QuadDobl_Sampling_Machine.Sample
      (standard_output,true,
       QuadDobl_Complex_Solutions.Head_Of(sols).all,hyp,newsol);
  end Compute_QuadDobl_Samples;

  procedure Test_DoblDobl_Sampling_Machine
              ( emb : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : natural32 ) is

  -- DESCRIPTION :
  --   Tests the sampling machine with double double precision arithmetic.

  -- ON ENTRY :
  --   emb      embedded polynomial system;
  --   sols     generic points on the solution set;
  --   dim      dimension of the solution set.

  begin
    new_line;
    put_line("Initializing the sampling machine ...");
    DoblDobl_Sampling_Machine.Initialize(emb);
    new_line;
    put_line("The embedded polynomial system : ");
    put(DoblDobl_Sampling_Machine.Embedded_System);
    new_line;
    put_line("Tuning the sampler ...");
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(standard_output,0);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    Compute_DoblDobl_Samples(sols,emb'last,integer32(dim));
  end Test_DoblDobl_Sampling_Machine;

  procedure Test_QuadDobl_Sampling_Machine
              ( emb : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : natural32 ) is

  -- DESCRIPTION :
  --   Tests the sampling machine with quad double precision arithmetic.

  -- ON ENTRY :
  --   emb      embedded polynomial system;
  --   sols     generic points on the solution set;
  --   dim      dimension of the solution set.

  begin
    new_line;
    put_line("Initializing the sampling machine ...");
    QuadDobl_Sampling_Machine.Initialize(emb);
    new_line;
    put_line("The embedded polynomial system : ");
    put(QuadDobl_Sampling_Machine.Embedded_System);
    new_line;
    put_line("Tuning the sampler ...");
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(standard_output,0);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    Compute_QuadDobl_Samples(sols,emb'last,integer32(dim));
  end Test_QuadDobl_Sampling_Machine;

  procedure Test_DoblDobl_Sampling is

  -- DESCRIPTION :
  --   Prompts the user for a witness set and takes samples.

    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    dim : natural32;

  begin
    DoblDobl_Read_Embedding(p,sols,dim);
    new_line;
    put("The dimension : "); put(dim,1); new_line;
    put("The degree : ");
    put(DoblDobl_Complex_Solutions.Length_Of(sols),1); new_line;
    put("Ambient dimension : "); put(p'last,1); new_line;
    Test_DoblDobl_Sampling_Machine(p.all,sols,dim);
  end Test_DoblDobl_Sampling;

  procedure Test_QuadDobl_Sampling is

  -- DESCRIPTION :
  --   Prompts the user for a witness set and takes samples.

    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    dim : natural32;

  begin
    QuadDobl_Read_Embedding(p,sols,dim);
    new_line;
    put("The dimension : "); put(dim,1); new_line;
    put("The degree : ");
    put(QuadDobl_Complex_Solutions.Length_Of(sols),1); new_line;
    put("Ambient dimension : "); put(p'last,1); new_line;
    Test_QuadDobl_Sampling_Machine(p.all,sols,dim);
  end Test_QuadDobl_Sampling;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Sampling a positive dimensional solution set ...");
    put_line("  1. in double double arithmetic; or");
    put_line("  2. in quad double arithmetic.");
    put("Type 1 or 2 to make your choice : "); Ask_Alternative(ans,"12");
    if ans = '1'
     then Test_DoblDobl_Sampling;
     else Test_QuadDobl_Sampling;
    end if;
  end Main;

begin
  Main;
end ts_ddqdsam;
