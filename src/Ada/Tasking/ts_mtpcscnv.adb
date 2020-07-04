with Communications_with_User;            use Communications_with_User;
with Multitasked_Path_Convolutions;      use Multitasked_Path_Convolutions;

procedure ts_mtpcscnv is

-- DESCRIPTION :
--   Tests multitasked tracking with predictor-corrector-shift
--   loops on homotopy systems of convolution circuits.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then launches
  --   the test corresponding to the selected precision.

    precision : constant character := Prompt_for_Precision;

  begin
    case precision is
      when '0' => Standard_Main(0);
      when '1' => DoblDobl_Main(0);
      when '2' => QuadDobl_Main(0);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtpcscnv;
