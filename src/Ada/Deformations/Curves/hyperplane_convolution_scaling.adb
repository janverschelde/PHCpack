with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_Polar;
with Hyperplane_Solution_Scaling;

package body Hyperplane_Convolution_Scaling is

-- 1-HOMOGENIZATION :

  procedure Adjust ( cff : in Standard_Complex_VecVecs.VecVec;
                     cst : in Standard_Complex_Vectors.Link_to_Vector;
                     sol : in Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    lnk : Standard_Complex_Vectors.Link_to_Vector;
    val : Complex_Number := cst(0);

  begin
    for k in sol'range loop
      lnk := cff(k);
      val := val + lnk(0)*sol(k);
    end loop;
    cst(0) := cst(0) - val;
  end Adjust;

  procedure Adjust ( cff : in DoblDobl_Complex_VecVecs.VecVec;
                     cst : in DoblDobl_Complex_Vectors.Link_to_Vector;
                     sol : in DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    val : Complex_Number := cst(0);

  begin
    for k in sol'range loop
      lnk := cff(k);
      val := val + lnk(0)*sol(k);
    end loop;
    cst(0) := cst(0) - val;
  end Adjust;

  procedure Adjust ( cff : in QuadDobl_Complex_VecVecs.VecVec;
                     cst : in QuadDobl_Complex_Vectors.Link_to_Vector;
                     sol : in QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    val : Complex_Number := cst(0);

  begin
    for k in sol'range loop
      lnk := cff(k);
      val := val + lnk(0)*sol(k);
    end loop;
    cst(0) := cst(0) - val;
  end Adjust;

  procedure Adjust_Last_Constant
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in Standard_Complex_Vectors.Vector ) is

    crc : constant Standard_Speelpenning_Convolutions.Link_to_Circuit
        := hom.crc(hom.crc'last);

  begin
    Adjust(crc.cff,crc.cst,sol);
  end Adjust_Last_Constant;

  procedure Adjust_Last_Constant
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in DoblDobl_Complex_Vectors.Vector ) is

    crc : constant DoblDobl_Speelpenning_Convolutions.Link_to_Circuit
        := hom.crc(hom.crc'last);

  begin
    Adjust(crc.cff,crc.cst,sol);
  end Adjust_Last_Constant;

  procedure Adjust_Last_Constant
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in QuadDobl_Complex_Vectors.Vector ) is

    crc : constant QuadDobl_Speelpenning_Convolutions.Link_to_Circuit
        := hom.crc(hom.crc'last);

  begin
    Adjust(crc.cff,crc.cst,sol);
  end Adjust_Last_Constant;

  procedure Adjust_Last_Radius
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System ) is

    homcrc : constant Standard_Speelpenning_Convolutions.Link_to_Circuit
           := hom.crc(hom.crc'last);
    abhcrc : constant Standard_Speelpenning_Convolutions.Link_to_Circuit
           := abh.crc(abh.crc'last);
    homlnk : constant Standard_Complex_Vectors.Link_to_Vector := homcrc.cst;
    abhlnk : constant Standard_Complex_Vectors.Link_to_Vector := abhcrc.cst;
    rad : double_float;

  begin
    rad := Standard_Complex_Numbers_Polar.Radius(homlnk(0));
    abhlnk(0) := Standard_Complex_Numbers.Create(rad);
  end Adjust_Last_Radius;

  procedure Adjust_Last_Radius
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System ) is

    homcrc : constant DoblDobl_Speelpenning_Convolutions.Link_to_Circuit
           := hom.crc(hom.crc'last);
    abhcrc : constant DoblDobl_Speelpenning_Convolutions.Link_to_Circuit
           := abh.crc(abh.crc'last);
    homlnk : constant DoblDobl_Complex_Vectors.Link_to_Vector := homcrc.cst;
    abhlnk : constant DoblDobl_Complex_Vectors.Link_to_Vector := abhcrc.cst;
    rad : double_double;

  begin
    rad := DoblDobl_Complex_Numbers_Polar.Radius(homlnk(0));
    abhlnk(0) := DoblDobl_Complex_Numbers.Create(rad);
  end Adjust_Last_Radius;

  procedure Adjust_Last_Radius
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System ) is

    homcrc : constant QuadDobl_Speelpenning_Convolutions.Link_to_Circuit
           := hom.crc(hom.crc'last);
    abhcrc : constant QuadDobl_Speelpenning_Convolutions.Link_to_Circuit
           := abh.crc(abh.crc'last);
    homlnk : constant QuadDobl_Complex_Vectors.Link_to_Vector := homcrc.cst;
    abhlnk : constant QuadDobl_Complex_Vectors.Link_to_Vector := abhcrc.cst;
    rad : quad_double;

  begin
    rad := QuadDobl_Complex_Numbers_Polar.Radius(homlnk(0));
    abhlnk(0) := QuadDobl_Complex_Numbers.Create(rad);
  end Adjust_Last_Radius;

  procedure Scale_and_Adjust
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector ) is

    crc : constant Standard_Speelpenning_Convolutions.Link_to_Circuit
        := hom.crc(hom.crc'last);

  begin
    Hyperplane_Solution_Scaling.Scale(sol);
    Adjust(crc.cff,crc.cst,sol);
  end Scale_and_Adjust;

  procedure Scale_and_Adjust
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector ) is

    crc : constant DoblDobl_Speelpenning_Convolutions.Link_to_Circuit
        := hom.crc(hom.crc'last);

  begin
    Hyperplane_Solution_Scaling.Scale(sol);
    Adjust(crc.cff,crc.cst,sol);
  end Scale_and_Adjust;

  procedure Scale_and_Adjust
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector ) is

    crc : constant QuadDobl_Speelpenning_Convolutions.Link_to_Circuit
        := hom.crc(hom.crc'last);

  begin
    Hyperplane_Solution_Scaling.Scale(sol);
    Adjust(crc.cff,crc.cst,sol);
  end Scale_and_Adjust;

-- MULTI-HOMOGENIZATION :

  procedure Adjust ( cff : in Standard_Complex_VecVecs.VecVec;
                     cst : in Standard_Complex_Vectors.Link_to_Vector;
                     sol : in Standard_Complex_Vectors.Vector;
                     idz : in Standard_Natural_Vectors.Link_to_Vector;
                     m,i : in integer32 ) is

    use Standard_Complex_Numbers;

    lnk : Standard_Complex_Vectors.Link_to_Vector;
    val : Complex_Number := cst(0);
    dim : constant integer32 := sol'last - m; -- original dimension
    idxcff : integer32 := cff'first-1;

  begin
    for k in sol'first..dim loop
      if integer32(idz(k)) = i then -- only use variables from the i-th set
        idxcff := idxcff + 1;
        lnk := cff(idxcff);
        val := val + lnk(0)*sol(k);
      end if;
    end loop;
    lnk := cff(idxcff+1);           -- the last coefficient is for z(i)
    val := val + lnk(0)*sol(dim+i);
    cst(0) := cst(0) - val;
  end Adjust;

  procedure Adjust ( cff : in DoblDobl_Complex_VecVecs.VecVec;
                     cst : in DoblDobl_Complex_Vectors.Link_to_Vector;
                     sol : in DoblDobl_Complex_Vectors.Vector;
                     idz : in Standard_Natural_Vectors.Link_to_Vector;
                     m,i : in integer32 ) is

    use DoblDobl_Complex_Numbers;

    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    val : Complex_Number := cst(0);
    dim : constant integer32 := sol'last - m; -- original dimension
    idxcff : integer32 := cff'first-1;

  begin
    for k in sol'first..dim loop
      if integer32(idz(k)) = i then -- only use variables from the i-th set
        idxcff := idxcff + 1;
        lnk := cff(idxcff);
        val := val + lnk(0)*sol(k);
      end if;
    end loop;
    lnk := cff(idxcff+1);           -- the last coefficient is for z(i)
    val := val + lnk(0)*sol(dim+i);
    cst(0) := cst(0) - val;
  end Adjust;

  procedure Adjust ( cff : in QuadDobl_Complex_VecVecs.VecVec;
                     cst : in QuadDobl_Complex_Vectors.Link_to_Vector;
                     sol : in QuadDobl_Complex_Vectors.Vector;
                     idz : in Standard_Natural_Vectors.Link_to_Vector;
                     m,i : in integer32 ) is

    use QuadDobl_Complex_Numbers;

    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    val : Complex_Number := cst(0);
    dim : constant integer32 := sol'last - m; -- original dimension
    idxcff : integer32 := cff'first-1;

  begin
    for k in sol'first..dim loop
      if integer32(idz(k)) = i then -- only use variables from the i-th set
        idxcff := idxcff + 1;
        lnk := cff(idxcff);
        val := val + lnk(0)*sol(k);
      end if;
    end loop;
    lnk := cff(idxcff+1);           -- the last coefficient is for z(i)
    val := val + lnk(0)*sol(dim+i);
    cst(0) := cst(0) - val;
  end Adjust;

  procedure Adjust_Last_Radii
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                m : in integer32 ) is

    homcrc,abhcrc : Standard_Speelpenning_Convolutions.Link_to_Circuit;
    homlnk,abhlnk : Standard_Complex_Vectors.Link_to_Vector;
    rad : double_float;

  begin
    for i in 1..m loop
      homcrc := hom.crc(hom.crc'last-m+i); homlnk := homcrc.cst;
      abhcrc := abh.crc(hom.crc'last-m+i); abhlnk := abhcrc.cst;
      rad := Standard_Complex_Numbers_Polar.Radius(homlnk(0));
      abhlnk(0) := Standard_Complex_Numbers.Create(rad);
    end loop;
  end Adjust_Last_Radii;

  procedure Adjust_Last_Radii
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                m : in integer32 ) is

    homcrc,abhcrc : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
    homlnk,abhlnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    rad : double_double;

  begin
    for i in 1..m loop
      homcrc := hom.crc(hom.crc'last-m+i); homlnk := homcrc.cst;
      abhcrc := abh.crc(hom.crc'last-m+i); abhlnk := abhcrc.cst;
      rad := DoblDobl_Complex_Numbers_Polar.Radius(homlnk(0));
      abhlnk(0) := DoblDobl_Complex_Numbers.Create(rad);
    end loop;
  end Adjust_Last_Radii;

  procedure Adjust_Last_Radii
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                m : in integer32 ) is

    homcrc,abhcrc : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
    homlnk,abhlnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    rad : quad_double;

  begin
    for i in 1..m loop
      homcrc := hom.crc(hom.crc'last-m+i); homlnk := homcrc.cst;
      abhcrc := abh.crc(hom.crc'last-m+i); abhlnk := abhcrc.cst;
      rad := QuadDobl_Complex_Numbers_Polar.Radius(homlnk(0));
      abhlnk(0) := QuadDobl_Complex_Numbers.Create(rad);
    end loop;
  end Adjust_Last_Radii;

  procedure Scale_and_Adjust
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                m : in integer32 ) is

    crc : Standard_Speelpenning_Convolutions.Link_to_Circuit;

  begin
    Hyperplane_Solution_Scaling.Scale(sol,idz,m);
    for i in 1..m loop
      crc := hom.crc(hom.crc'last-m+i);
      Adjust(crc.cff,crc.cst,sol,idz,m,i);
    end loop;
  end Scale_and_Adjust;

  procedure Scale_and_Adjust
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                m : in integer32 ) is

    crc : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;

  begin
    Hyperplane_Solution_Scaling.Scale(sol,idz,m);
    for i in 1..m loop
      crc := hom.crc(hom.crc'last-m+i);
      Adjust(crc.cff,crc.cst,sol,idz,m,i);
    end loop;
  end Scale_and_Adjust;

  procedure Scale_and_Adjust
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                m : in integer32 ) is

    crc : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;

  begin
    Hyperplane_Solution_Scaling.Scale(sol,idz,m);
    for i in 1..m loop
      crc := hom.crc(hom.crc'last-m+i);
      Adjust(crc.cff,crc.cst,sol,idz,m,i);
    end loop;
  end Scale_and_Adjust;

end Hyperplane_Convolution_Scaling;
