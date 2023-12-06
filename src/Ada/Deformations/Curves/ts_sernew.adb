with Communications_with_User;           use Communications_with_User;
with Test_SD_Newton_Matrix_Series;
with Test_DD_Newton_Matrix_Series;
with Test_QD_Newton_Matrix_Series;

procedure ts_sernew is

-- DESCRIPTION :
--   Tests the application of Newton's method to compute series solutions.

  prc : constant character := Prompt_for_Precision;

begin
  case prc is 
    when '0' => Test_SD_Newton_Matrix_Series.Main;
    when '1' => Test_DD_Newton_Matrix_Series.Main;
    when '2' => Test_QD_Newton_Matrix_Series.Main;
    when others => null;
  end case;
end ts_sernew;
