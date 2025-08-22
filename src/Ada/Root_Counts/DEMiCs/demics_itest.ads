with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with demics_input_data;

package demics_itest is

-- DESCRIPTION :
--   Defines classes to store integer data of linear programs.
--   Translation initiated by g++ -c -fdump-ada-spec iTest.h.

  package class_uData is

    type uData;

    type uData is record
      nfN : integer32;
      next : access uData;
      prev : access uData;
      fNext : access uData;
      supLab : integer32;
      red : double_float;
      dir : Standard_Floating_Vectors.Link_to_Vector;
    end record;

    function new_uData return uData;

    -- DESCRIPTION :
    --   Returns an empty data structure.

    procedure delete_uData ( this : access uData );

    -- DESCRIPTION :
    --   Deallocates the this.dir vector.

    procedure create ( this : access uData;
                       depth : in integer32; dim : in integer32 );

    procedure init ( this : access uData );

    procedure getDir ( this : access uData; val : in double_float;
                       idx : in integer32 );

    procedure getRed ( this : access uData; val : out double_float;
                       idx : in integer32 );

    procedure info_dirRed ( this : access uData );

  end class_uData;

  use class_uData;

  package class_inifData is

    type inifData is record
      head : access uData;
      fHead : access uData;
      last : access uData;
    end record;

    function new_inifData return inifData;

    procedure delete_inifData ( this : access inifData );

    procedure create ( this : access inifData;
                       length : in integer32; depth : in integer32;
                       dim : in integer32 );

    procedure get_info
      ( this : access inifData;
        data : access demics_input_data.class_dataSet.dataSet;
        lifting : in Standard_Floating_Vectors.Link_to_Vector;
        termSet : in Standard_Integer_Vectors.Link_to_Vector;
        termStart : in Standard_Integer_Vectors.Link_to_Vector;
        depth : in integer32; dim : in integer32; supN : in integer32 );

    procedure info_all_dirRed ( this : access inifData );

    procedure info_feasIdx ( this : access inifData );

    procedure info_fNext ( this : access inifData );

    procedure info_next ( this : access inifData );

    procedure info_prev ( this : access inifData );

  end class_inifData;

  use class_inifData;

  package class_iLvData is

    type iLvData is record
      rspLen : integer32;
      inifLen : integer32;
      inif : access inifData;
      rsp : Standard_Integer_Vectors.Link_to_Vector;
    end record;

    function new_iLvData return iLvData;

    procedure delete_iLvData ( this : access iLvData );

    procedure create ( this : access iLvData;
                       depth : in integer32; supN : in integer32;
                       dim : in integer32; termMax : in integer32 );

    procedure getInit
      ( this : access iLvData;
        Data : access demics_input_data.class_dataSet.dataSet;
        lifting : in Standard_Floating_Vectors.Link_to_Vector;
        termSet : in Standard_Integer_Vectors.Link_to_Vector;
        termStart : in Standard_Integer_Vectors.Link_to_Vector;
        dim : in integer32; supN : in integer32 );

    procedure init
      ( this : access iLvData; supN : in integer32; depth : in integer32;
        preRsp : in Standard_Integer_Vectors.Link_to_Vector );

    procedure info_rsp ( this : access iLvData );

    procedure info_all_dirRed ( this : access iLvData );

    procedure info_feasIdx ( this : access iLvData; depth : in integer32 );

    procedure info_all_feasIdx ( this : access iLvData );

  end class_iLvData;

end demics_itest;
