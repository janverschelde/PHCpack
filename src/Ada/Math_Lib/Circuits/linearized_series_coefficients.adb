with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;

package body Linearized_Series_Coefficients is

  procedure Delinearize ( vy,yv : in Standard_Complex_VecVecs.VecVec ) is
  begin
    for k in vy'range loop
      declare
        vyk : constant Standard_Complex_Vectors.Link_to_Vector := vy(k);
        left : Standard_Complex_Vectors.Link_to_Vector;
      begin
        for i in yv'range loop  -- vyk holds k-th coefficient of all series
          left := yv(i);        -- so we assign to coefficients of series i
          left(k) := vyk(i);    -- at position k the i-th value of vyk
        end loop;
      end;
    end loop;
  end Delinearize;

  procedure Delinearize ( vy,yv : in DoblDobl_Complex_VecVecs.VecVec ) is
  begin
    for k in vy'range loop
      declare
        vyk : constant DoblDobl_Complex_Vectors.Link_to_Vector := vy(k);
        left : DoblDobl_Complex_Vectors.Link_to_Vector;
      begin
        for i in yv'range loop  -- vyk holds k-th coefficient of all series
          left := yv(i);        -- so we assign to coefficients of series i
          left(k) := vyk(i);    -- at position k the i-th value of vyk
        end loop;
      end;
    end loop;
  end Delinearize;

  procedure Delinearize ( vy,yv : in QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    for k in vy'range loop
      declare
        vyk : constant QuadDobl_Complex_Vectors.Link_to_Vector := vy(k);
        left : QuadDobl_Complex_Vectors.Link_to_Vector;
      begin
        for i in yv'range loop  -- vyk holds k-th coefficient of all series
          left := yv(i);        -- so we assign to coefficients of series i
          left(k) := vyk(i);    -- at position k the i-th value of vyk
        end loop;
      end;
    end loop;
  end Delinearize;

end Linearized_Series_Coefficients;
