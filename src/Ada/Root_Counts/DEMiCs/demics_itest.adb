package body demics_itest is

  package body class_uData is

    function new_uData return uData is

      res : uData;

    begin
      return res;
    end new_uData;

    procedure delete_uData ( this : in Link_to_uData ) is
    begin
      null;
    end delete_uData;

    procedure create ( this : in Link_to_uData;
                       depth : in integer32; dim : in integer32 ) is
    begin
      null;
    end create;

    procedure init ( this : in Link_to_uData ) is
    begin
      null;
    end init;

    procedure getDir ( this : in Link_to_uData; val : in double_float;
                       idx : in integer32 ) is
    begin
      null;
    end getDir;

    procedure getRed ( this : in Link_to_uData; val : out double_float;
                       idx : in integer32 ) is
    begin
      null;
    end getRed;

    procedure info_dirRed ( this : in Link_to_uData ) is
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

    procedure delete_inifData ( this : in Link_to_inifData ) is
    begin
      null;
    end delete_inifData;

    procedure create ( this : in Link_to_inifData;
                       length : in integer32; depth : in integer32;
                       dim : in integer32 ) is
    begin
      null;
    end create;

    procedure get_info
                ( this : in Link_to_inifData;
                  data : in demics_input_data.class_dataSet.dataSet;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector;
                  termSet : in Standard_Integer_Vectors.Link_to_Vector;
                  termStart : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32; dim : in integer32;
                  supN : in integer32 ) is

    begin
      null;
    end get_info;

    procedure info_all_dirRed ( this : in Link_to_inifData ) is
    begin
      null;
    end info_all_dirRed;

    procedure info_feasIdx ( this : in Link_to_inifData ) is
    begin
      null;
    end info_feasIdx;

    procedure info_fNext ( this : in Link_to_inifData ) is
    begin
      null;
    end info_fNext;

    procedure info_next ( this : in Link_to_inifData ) is
    begin
      null;
    end info_next;

    procedure info_prev ( this : in Link_to_inifData ) is
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

    procedure delete_iLvData ( this : in Link_to_iLvData ) is
    begin
      null;
    end delete_iLvData;

    procedure create ( this : in Link_to_iLvData;
                       depth : in integer32; supN : in integer32;
                       dim : in integer32; termMax : in integer32 ) is
    begin
      null;
    end create;

    procedure getInit
                ( this : in Link_to_iLvData;
                  data : in demics_input_data.class_dataSet.dataSet;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector;
                  termSet : in Standard_Integer_Vectors.Link_to_Vector;
                  termStart : in Standard_Integer_Vectors.Link_to_Vector;
                  dim : in integer32; supN : in integer32 ) is
    begin
      null;
    end getInit;

    procedure init ( this : in Link_to_iLvData;
                     supN : in integer32; depth : in integer32;
                     preRsp : in Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      null;
    end init;

    procedure info_rsp ( this : in Link_to_iLvData ) is
    begin
      null;
    end info_rsp;

    procedure info_all_dirRed ( this : in Link_to_iLvData ) is
    begin
      null;
    end info_all_dirRed;

    procedure info_feasIdx ( this : in Link_to_iLvData;
                             depth : in integer32 ) is
    begin
      null;
    end info_feasIdx;

    procedure info_all_feasIdx ( this : in Link_to_iLvData ) is
    begin
      null;
    end info_all_feasIdx;

  end class_iLvData;

end demics_itest;
