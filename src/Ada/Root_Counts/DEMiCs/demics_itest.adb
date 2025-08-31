with Ada.text_io;                       use Ada.text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;

package body demics_itest is

  package body class_uData is

    function new_uData return uData is

      res : uData;

    begin
      res.nfN := 0;
      res.supLab := 0;
      res.red := 0.0;
      res.dir := null;
      res.next := null;
      res.prev := null;
      res.fNext := null;
      return res;
    end new_uData;

    procedure delete_uData ( this : in Link_to_uData ) is
    begin
      Standard_Floating_Vectors.clear(this.dir);
    end delete_uData;

    procedure create ( this : in Link_to_uData;
                       depth : in integer32; dim : in integer32;
                       vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_itest.class_uData.create, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      this.nfN := dim;
      this.dir := new Standard_Floating_Vectors.Vector'(0..this.nfN-1 => 0.0);
    end create;

    procedure init ( this : in Link_to_uData ) is

      use Standard_Floating_Vectors;

    begin
      if this.dir /= null then
        for i in this.dir'range loop
          this.dir(i) := 0.0;
        end loop;
      end if;
      this.red := 0.0;
      this.supLab := 0;
    end init;

    procedure getDir ( this : in Link_to_uData; val : in double_float;
                       idx : in integer32 ) is
    begin
      this.dir(idx) := val;
    end getDir;

    procedure getRed ( this : in Link_to_uData; val : in double_float ) is
    begin
      this.red := val;
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
      res.head := null;
      res.fhead := null;
      res.last := null;
      return res;
    end new_inifdata;

    procedure delete_inifData ( this : in Link_to_inifData ) is
    begin
      null;
    end delete_inifData;

    procedure create ( this : in Link_to_inifData;
                       length : in integer32; depth : in integer32;
                       dim : in integer32; vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_itest.class_inifData.create, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      for i in 0..length-1 loop
        declare
          newData : constant Link_to_uData := new uData'(new_uData);
        begin
          class_uData.create(newData,depth,dim,vrblvl-1);
          if this.last /= null then
            this.last.next := newData;
            this.last.fNext := newData;
            newData.prev := this.last;
          else
            this.head := newData;
            this.fHead := newData;
          end if;
          this.last := newData;
        end;
      end loop;
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
      res.rspLen := 0;
      res.inifLen := 0;
      res.inif := null;
      res.rsp := null;
      return res;
    end new_iLvData;

    procedure delete_iLvData ( this : in Link_to_iLvData ) is
    begin
      null;
    end delete_iLvData;

    procedure create ( this : in Link_to_iLvData;
                       depth : in integer32; supN : in integer32;
                       dim : in integer32; termMax : in integer32;
                       vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_itest.class_iLvData.create, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      this.rspLen := supN - depth - 1;
      this.inifLen := supN;
      this.inif := new Array_of_inifData(0..supN-1);
      this.rsp := new Standard_Integer_Vectors.Vector(0..supN-1);
      for i in 0..this.inifLen-1 loop
        this.inif(i) := new inifData'(new_inifData);
        class_inifData.create(this.inif(i),termMax,depth,dim,vrblvl-1);
      end loop;
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
