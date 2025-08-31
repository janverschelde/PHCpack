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
    --   Returns an empty data structure,
    --   with values set to zero or to null.

    procedure delete_uData ( this : in Link_to_uData );

    -- DESCRIPTION :
    --   Deallocates the this.dir vector.

    procedure create ( this : in Link_to_uData;
                       depth : in integer32; dim : in integer32;
                       vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Sets the value of this.nfN to dim 
    --   and allocates memory for this.dir.

    procedure init ( this : in Link_to_uData );

    -- DESCRIPTION :
    --   Sets all elements of this.dir to zero,
    --   and sets this.supLab and this.red to zero as well.

    procedure getDir ( this : in Link_to_uData; val : in double_float;
                       idx : in integer32 );

    -- DESCRIPTION :
    --   Sets the value of this.dir at index idx to val.

    procedure getRed ( this : in Link_to_uData; val : in double_float );

    -- DESCRIPTION :
    --   Sets this.red to the value val.
    --   The original getRed contained a superfluous parameter idx,
    --   which was unused and therefore has been removed.

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

    -- DESCRIPTION :
    --   Returns an inifData object will all pointers set to null,
    --   which must be invoked each time a new Link_to_uData is made.

    procedure delete_inifData ( this : in Link_to_inifData );

    procedure create ( this : in Link_to_inifData;
                       length : in integer32; depth : in integer32;
                       dim : in integer32; vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Makes space for a linked list of the given length
    --   with uData nodes for the given depth and dim.

    procedure get_info
                ( this : in Link_to_inifData;
                  data : in demics_input_data.class_dataSet.dataSet;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector;
                  termSet : in Standard_Integer_Vectors.Link_to_Vector;
                  termStart : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32; dim : in integer32;
                  supN : in integer32; vrblvl : in integer32 := 0 );

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
      inif : Link_to_Array_of_inifData;
      rsp : Standard_Integer_Vectors.Link_to_Vector;
    end record;

    type Link_to_iLvData is access iLvData;
    type Array_of_iLvData is array ( integer32 range <> ) of Link_to_iLvData;
    type Link_to_Array_of_iLvData is access Array_of_iLvData;

    function new_iLvData return iLvData;

    -- DESCRIPTION :
    --   Returns a record with zero and null values.

    procedure delete_iLvData ( this : in Link_to_iLvData );

    procedure create ( this : in Link_to_iLvData;
                       depth : in integer32; supN : in integer32;
                       dim : in integer32; termMax : in integer32;
                       vrblvl : in integer32 := 0 );

    -- DESCRIPTION :
    --   Sets the values of this.rspLen and this.inifLen, allocates memory
    --   for this.rsp and makes the values of the this.inif array.

    procedure getInit
                ( this : in Link_to_iLvData;
                  Data : in demics_input_data.class_dataSet.dataSet;
                  lifting : in Standard_Floating_Vectors.Link_to_Vector;
                  termSet : in Standard_Integer_Vectors.Link_to_Vector;
                  termStart : in Standard_Integer_Vectors.Link_to_Vector;
                  dim : in integer32; supN : in integer32;
                  vrblvl : in integer32 := 0 );

    procedure init ( this : in Link_to_iLvData;
                     supN : in integer32; depth : in integer32;
                     preRsp : in Standard_Integer_Vectors.Link_to_Vector );

    procedure info_rsp ( this : in Link_to_iLvData );

    procedure info_all_dirRed ( this : in Link_to_iLvData );

    procedure info_feasIdx ( this : in Link_to_iLvData; depth : in integer32 );

    procedure info_all_feasIdx ( this : in Link_to_iLvData );

  end class_iLvData;

end demics_itest;
