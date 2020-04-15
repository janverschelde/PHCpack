with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;

package body Hyperplane_Convolution_Scaling is

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

end Hyperplane_Convolution_Scaling;
