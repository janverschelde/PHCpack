package body demics_itest is

  package body class_uData is

    function new_uData return uData is

      res : uData;

    begin
      return res;
    end new_uData;

    procedure delete_uData ( this : access uData ) is
    begin
      null;
    end delete_uData;

    procedure create ( this : access uData;
                       depth : in integer32; dim : in integer32 ) is
    begin
      null;
    end create;

    procedure init ( this : access uData ) is
    begin
      null;
    end init;

    procedure getDir ( this : access uData; val : in double_float;
                       idx : in integer32 ) is
    begin
      null;
    end getDir;

    procedure getRed ( this : access uData; val : out double_float;
                       idx : in integer32 ) is
    begin
      null;
    end getRed;

    procedure info_dirRed ( this : access uData ) is
    begin
      null;
    end info_dirRed;

  end class_uData;

  package body class_inifData is

    function new_inifData return inifData is

      res : inifData;

    begin
      return res;
    end new_inifdata;

    procedure delete_inifData ( this : access inifData ) is
    begin
      null;
    end delete_inifData;

    procedure create ( this : access inifData;
                       length : in integer32; depth : in integer32;
                       dim : in integer32 ) is
    begin
      null;
    end create;

    procedure get_info
      ( this : access inifData;
        data : access demics_input_data.class_dataSet.dataSet;
        lifting : in Standard_Floating_Vectors.Link_to_Vector;
        termSet : in Standard_Integer_Vectors.Link_to_Vector;
        termStart : in Standard_Integer_Vectors.Link_to_Vector;
        depth : in integer32; dim : in integer32; supN : in integer32 ) is

    begin
      null;
    end get_info;

    procedure info_all_dirRed ( this : access inifData ) is
    begin
      null;
    end info_all_dirRed;

    procedure info_feasIdx ( this : access inifData ) is
    begin
      null;
    end info_feasIdx;

    procedure info_fNext ( this : access inifData ) is
    begin
      null;
    end info_fNext;

    procedure info_next ( this : access inifData ) is
    begin
      null;
    end info_next;

    procedure info_prev ( this : access inifData ) is
    begin
      null;
    end info_prev;

  end class_inifData;

  package body class_iLvData is

    function new_iLvData return iLvData is

      res : iLvData;

    begin
      return res;
    end new_iLvData;

    procedure delete_iLvData ( this : access iLvData ) is
    begin
      null;
    end delete_iLvData;

    procedure create ( this : access iLvData;
                       depth : in integer32; supN : in integer32;
                       dim : in integer32; termMax : in integer32 ) is
    begin
      null;
    end create;

    procedure getInit
      ( this : access iLvData;
        Data : access demics_input_data.class_dataSet.dataSet;
        lifting : in Standard_Floating_Vectors.Link_to_Vector;
        termSet : in Standard_Integer_Vectors.Link_to_Vector;
        termStart : in Standard_Integer_Vectors.Link_to_Vector;
        dim : in integer32; supN : in integer32 ) is
    begin
      null;
    end getInit;

    procedure init
      ( this : access iLvData; supN : in integer32; depth : in integer32;
        preRsp : in Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      null;
    end init;

    procedure info_rsp ( this : access iLvData ) is
    begin
      null;
    end info_rsp;

    procedure info_all_dirRed ( this : access iLvData ) is
    begin
      null;
    end info_all_dirRed;

    procedure info_feasIdx ( this : access iLvData; depth : in integer32 ) is
    begin
      null;
    end info_feasIdx;

    procedure info_all_feasIdx ( this : access iLvData ) is
    begin
      null;
    end info_all_feasIdx;

  end class_iLvData;

end demics_itest;
