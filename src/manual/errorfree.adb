with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with TripDobl_Complex_Vectors;
with TripDobl_Random_Vectors;
with TripDobl_Complex_Vector_Norms;      use TripDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with PentDobl_Complex_Vectors;
with PentDobl_Random_Vectors;
with PentDobl_Complex_Vector_Norms;      use PentDobl_Complex_Vector_Norms;
with OctoDobl_Complex_Vectors;
with OctoDobl_Random_Vectors;
with OctoDobl_Complex_Vector_Norms;      use OctoDobl_Complex_Vector_Norms;
with DecaDobl_Complex_Vectors;
with DecaDobl_Random_Vectors;
with DecaDobl_Complex_Vector_Norms;      use DecaDobl_Complex_Vector_Norms;
with HexaDobl_Complex_Vectors;
with HexaDobl_Random_Vectors;
with HexaDobl_Complex_Vector_Norms;      use HexaDobl_Complex_Vector_Norms;

procedure errorfree is

-- DESCRIPTION :
--   If in any numerical computation with multiple doubles,
--   the first double equals the result without representation error, 
--   then the second double equals the multiple double precision
--   used to compute the result.

  procedure write ( first,second : in double_float ) is

  -- DESCRIPTION :
  --   Writes the two doubles, checking the sign of second.

  begin
    put(first);
    if second > 0.0
     then put(" +"); put(second);
     else put(" -"); put(-second);
    end if;
    new_line;
  end write;

  procedure DoblDobl_Two_Norm ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a complex vector of dimension dim 
  --   in double double precision and writes its 2-norm.

    rnv : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    nrm : constant double_double := Norm2(rnv);
    hi_norm : constant double_float := hi_part(nrm);
    lo_norm : constant double_float := lo_part(nrm);

  begin
    put("double double :"); write(hi_norm,lo_norm);
  end DoblDobl_Two_Norm;

  procedure TripDobl_Two_Norm ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a complex vector of dimension dim 
  --   in triple double precision and writes its 2-norm.

    rnv : constant TripDobl_Complex_Vectors.Vector(1..dim)
        := TripDobl_Random_Vectors.Random_Vector(1,dim);
    nrm : constant triple_double := Norm2(rnv);
    hi_norm : constant double_float := hi_part(nrm);
    mi_norm : constant double_float := mi_part(nrm);

  begin
    put("triple double :"); write(hi_norm,mi_norm);
  end TripDobl_Two_Norm;

  procedure QuadDobl_Two_Norm ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a complex vector of dimension dim 
  --   in quad double precision and writes its 2-norm.

    rnv : constant QuadDobl_Complex_Vectors.Vector(1..dim)
        := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    nrm : constant quad_double := Norm2(rnv);
    hihi_norm : constant double_float := hihi_part(nrm);
    lohi_norm : constant double_float := lohi_part(nrm);

  begin
    put("  quad double :"); write(hihi_norm,lohi_norm);
  end QuadDobl_Two_Norm;

  procedure PentDobl_Two_Norm ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a complex vector of dimension dim 
  --   in penta double precision and writes its 2-norm.

    rnv : constant PentDobl_Complex_Vectors.Vector(1..dim)
        := PentDobl_Random_Vectors.Random_Vector(1,dim);
    nrm : constant penta_double := Norm2(rnv);
    thumb_norm : constant double_float := thumb_part(nrm);
    index_norm : constant double_float := index_part(nrm);

  begin
    put(" penta double :"); write(thumb_norm,index_norm);
  end PentDobl_Two_Norm;

  procedure OctoDobl_Two_Norm ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a complex vector of dimension dim 
  --   in octo double precision and writes its 2-norm.

    rnv : constant OctoDobl_Complex_Vectors.Vector(1..dim)
        := OctoDobl_Random_Vectors.Random_Vector(1,dim);
    nrm : constant octo_double := Norm2(rnv);
    hihihi_norm : constant double_float := hihihi_part(nrm);
    lohihi_norm : constant double_float := lohihi_part(nrm);

  begin
    put("  octo double :"); write(hihihi_norm,lohihi_norm);
  end OctoDobl_Two_Norm;

  procedure DecaDobl_Two_Norm ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a complex vector of dimension dim 
  --   in deca double precision and writes its 2-norm.

    rnv : constant DecaDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Random_Vectors.Random_Vector(1,dim);
    nrm : constant deca_double := Norm2(rnv);
    rtb_norm : constant double_float := thumb_right(nrm);
    rix_norm : constant double_float := index_right(nrm);

  begin
    put("  deca double :"); write(rtb_norm,rix_norm);
  end DecaDobl_Two_Norm;

  procedure HexaDobl_Two_Norm ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a complex vector of dimension dim 
  --   in hexa double precision and writes its 2-norm.

    rnv : constant HexaDobl_Complex_Vectors.Vector(1..dim)
        := HexaDobl_Random_Vectors.Random_Vector(1,dim);
    nrm : constant hexa_double := Norm2(rnv);
    hihihihi_norm : constant double_float := hihihihi_part(nrm);
    lohihihi_norm : constant double_float := lohihihi_part(nrm);

  begin
    put("  hexa double :"); write(hihihihi_norm,lohihihi_norm);
  end HexaDobl_Two_Norm;

  procedure Main is

  begin
    new_line;
    put_line("Computing the 2-norm of a vector of dimension 64");
    put_line("of random complex numbers on the unit circle equals 8.");
    put_line("Observe the second double of the multiple double 2-norm.");
    new_line;
    DoblDobl_Two_Norm(64);
    TripDobl_Two_Norm(64);
    QuadDobl_Two_Norm(64);
    PentDobl_Two_Norm(64);
    OctoDobl_Two_Norm(64);
    DecaDobl_Two_Norm(64);
    HexaDobl_Two_Norm(64);
  end Main;

begin
  Main;
end errorfree;
