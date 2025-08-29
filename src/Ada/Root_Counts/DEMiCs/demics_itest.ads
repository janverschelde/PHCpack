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
    type Link_to_uData is access uData;
    type Array_of_uData is array ( integer32 range <> ) of Link_to_uData;
    type Link_to_Array_of_uData is access Array_of_uData;

    type uData is record
      nfN : integer32;
      next : Link_to_uData;
      prev : Link_to_uData;
      fNext : Link_to_uData;
      supLab : integer32;
      red : double_float;
      dir : Standard_Floating_Vectors.Link_to_Vector;
    end record;

    function new_uData return uData;

    -- DESCRIPTION :
    --   Returns an empty data structure.

    procedure delete_uData ( this : in Link_to_uData );

    -- DESCRIPTION :
    --   Deallocates the this.dir vector.

    procedure create ( this : in Link_to_uData;
                       depth : in integer32; dim : in integer32 );

    procedure init ( this : in Link_to_uData );

    procedure getDir ( this : in Link_to_uData; val : in double_float;
                       idx : in integer32 );

    procedure getRed ( this : in Link_to_uData; val : out double_float;
                       idx : in integer32 );

    procedure info_dirRed ( this : in Link_to_uData );

  end class_uData;

  use class_uData;

  package class_inifData is

    type inifData is record
      head : Link_to_uData;
      fHead : Link_to_uData;
      last : Link_to_uData;
    end record;

    type Link_to_inifData is access inifData;
    type Array_of_inifData is array ( integer32 range <> ) of Link_to_inifData;
    type Link_to_Array_of_inifData is access Array_of_inifData;

    function new_inifData return inifData;

    procedure delete_inifData ( this : in Link_to_inifData );

    procedure create ( this : in Link_to_inifData;
                       length : in integer32; depth : in integer32;
                       dim : in integer32 );

    procedure get_info
                ( this : in Link_to_inifData;
                  data : in demics_input_data.class_dataSet.dataSet;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector;
                  termSet : in Standard_Integer_Vectors.Link_to_Vector;
                  termStart : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32; dim : in integer32;
                  supN : in integer32 );

    procedure info_all_dirRed ( this : in Link_to_inifData );

    procedure info_feasIdx ( this : in Link_to_inifData );

    procedure info_fNext ( this : in Link_to_inifData );

    procedure info_next ( this : in Link_to_inifData );

    procedure info_prev ( this : in Link_to_inifData );

  end class_inifData;

  use class_inifData;

  package class_iLvData is

    type iLvData is record
      rspLen : integer32;
      inifLen : integer32;
      inif : Link_to_inifData;
      rsp : Standard_Integer_Vectors.Link_to_Vector;
    end record;

    type Link_to_iLvData is access iLvData;
    type Array_of_iLvData is array ( integer32 range <> ) of Link_to_iLvData;
    type Link_to_Array_of_iLvData is access Array_of_iLvData;

    function new_iLvData return iLvData;

    procedure delete_iLvData ( this : in Link_to_iLvData );

    procedure create ( this : in Link_to_iLvData;
                       depth : in integer32; supN : in integer32;
                       dim : in integer32; termMax : in integer32 );

    procedure getInit
                ( this : in Link_to_iLvData;
                  Data : in demics_input_data.class_dataSet.dataSet;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector;
                  termSet : in Standard_Integer_Vectors.Link_to_Vector;
                  termStart : in Standard_Integer_Vectors.Link_to_Vector;
                  dim : in integer32; supN : in integer32 );

    procedure init ( this : in Link_to_iLvData;
                     supN : in integer32; depth : in integer32;
                     preRsp : in Standard_Integer_Vectors.Link_to_Vector );

    procedure info_rsp ( this : in Link_to_iLvData );

    procedure info_all_dirRed ( this : in Link_to_iLvData );

    procedure info_feasIdx ( this : in Link_to_iLvData; depth : in integer32 );

    procedure info_all_feasIdx ( this : in Link_to_iLvData );

  end class_iLvData;

end demics_itest;
