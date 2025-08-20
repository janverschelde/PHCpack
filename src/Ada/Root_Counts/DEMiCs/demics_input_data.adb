package body demics_input_data is

  package body class_dataSet is

    function new_dataSet return dataSet is

      res : dataSet;

    begin
      res.dim := 0;
      return res;
    end new_dataSet;

    procedure delete_dataSet ( this : access dataSet ) is
    begin
      null;
    end delete_dataSet;

    procedure support_in ( this : access dataSet;
                           rowIdx : integer32; colIdx : integer32;
                           elem : double_float ) is
    begin
      null;
    end support_in;

    function support_out ( this : access dataSet;
                           rowIdx : integer32; colIdx : integer32 )
                         return double_float is

      res : double_float := 0.0;

    begin
      return res;
    end support_out;

    procedure getInputFile ( this : access dataSet;
                             inputFile : Link_to_String ) is
    begin
      null;
    end getInputFile;

    procedure info_preamble ( this : access dataSet ) is
    begin
      null;
    end info_preamble;

    procedure info_supports ( this : access dataSet ) is
    begin
      null;
    end info_supports;

  end class_dataSet;

end demics_input_data;
