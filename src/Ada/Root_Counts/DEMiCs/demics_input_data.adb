with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;

package body demics_input_data is

-- All indexing in the original code started at zero,
-- while all vectors start at one.

  package body class_dataSet is

    function new_dataSet return dataSet is

      res : dataSet;

    begin
      res.dim := 0;
      res.supN := 0;
      res.termSumNum := 0;
      res.termMax := 0;
      res.typeMax := 0;
      res.termSet := null;
      res.termStart := null;
      res.supType := null;
      res.support := null;
      res.coef := null;
      res.outFile := null;
      return res;
    end new_dataSet;

    procedure delete_dataSet ( this : in out dataSet ) is
    begin
      this.dim := 0;
      this.supN := 0;
      this.termSumNum := 0;
      this.termMax := 0;
      this.typeMax := 0;
      Standard_Integer_Vectors.Clear(this.termSet);
      Standard_Integer_Vectors.Clear(this.termStart);
      Standard_Integer_Vectors.Clear(this.supType);
      Standard_Floating_Vectors.Clear(this.support);
      Standard_Floating_Vectors.Clear(this.coef);
      String_Splitters.Clear(this.outFile);
    end delete_dataSet;

    procedure support_in ( this : in dataSet;
                           rowIdx : in integer32; colIdx : in integer32;
                           elem : in double_float ) is

      idx : constant integer32 := this.dim*rowIdx + colIdx;

    begin
      this.support(idx) := elem; 
    end support_in;

    function support_out ( this : in dataSet;
                           rowIdx : in integer32; colIdx : in integer32 )
                         return double_float is

      idx : constant integer32 := this.dim*rowIdx + colIdx;

    begin
      return this.support(idx);
    end support_out;

    function makeOutputFile ( inputFile : Link_to_String )
                            return Link_to_String is

      res : Link_to_String;
      idx : integer := inputFile'last;

    begin
      for i in inputFile'range loop
        if inputFile(i) = '.'
         then idx := i; exit;
        end if;  
      end loop;
      if inputFile(idx) = '.' then
        declare
          outstring : constant string := inputFile(1..idx) & "out";
        begin
          res := new string'(outstring); 
        end;
      else
        declare
          outstring : constant string := inputFile(1..idx) & ".out";
        begin
          res := new string'(outstring); 
        end;
      end if;
      return res;
    end makeOutputFile;

    procedure parse_preamble ( file : in file_type;
                               this : in out dataSet;
                               fail : out boolean ) is

      ch : character;

    begin
      fail := false;
      while not end_of_file(file) loop -- look for 'Dim ='
        get(file,ch);
        exit when ch = '=';
      end loop;
      if ch /= '=' then
        fail := true;
      else
        get(file,this.dim);
        while not end_of_file(file) loop -- look for 'Support ='
          get(file,ch);
          exit when ch = '=';
        end loop;
        if ch /= '=' then
          fail := true;
        else
          get(file,this.supN);
          while not end_of_file(file) loop -- look for 'Elem ='
            get(file,ch);
            exit when ch = '=';
          end loop;
          if ch /= '=' then
            fail := true;
          else
            this.termSet
              := new Standard_Integer_Vectors.Vector(0..this.supN-1);
            for i in 0..this.supN-1 loop
              get(file,this.termSet(i));
            end loop;
            while not end_of_file(file) loop -- look for 'Type ='
              get(file,ch);
              exit when ch = '=';
            end loop;
            if ch /= '=' then
              fail := true;
            else
              this.supType
                := new Standard_Integer_Vectors.Vector(0..this.supN-1);
              for i in 0..this.supN-1 loop
                get(file,this.supType(i));
              end loop;
            end if;
          end if;
        end if;
      end if;
      this.termMax := this.termSet(0);
      this.typeMax := this.supType(0);
      for i in 1..this.supN-1 loop
        if this.termSet(i) > this.termMax
         then this.termMax := this.termSet(i);
        end if;
        if this.supType(i) > this.typeMax
         then this.typeMax := this.supType(i);
        end if;
      end loop;
    end parse_preamble;

    procedure parse_supports ( file : in file_type;
                               this : in out dataSet; fail : out boolean ) is
      idx : integer32 := 0;
      size : integer32;

    begin
      fail := false;
      this.termSumNum := 0;
      this.termStart := new Standard_Integer_Vectors.Vector(0..this.supN);
      this.termStart(0) := 1;
      for k in 0..this.supN-1 loop
        this.termSumNum := this.termSumNum + this.termSet(k);
        this.termStart(k+1) := this.termSumNum;
      end loop;
      this.termStart(this.supN) := this.termSumNum;
      size := this.termSumNum*this.dim;
      this.support := new Standard_Floating_Vectors.Vector(0..size-1);
      for i in 0..this.supN-1 loop
        for j in 1..this.termSet(i) loop
          for k in 1..this.dim loop
            get(file,this.support(idx));
            idx := idx + 1;
          end loop;
        end loop;
      end loop;
    exception
      when others => fail := true; return;
    end parse_supports;

    procedure getInputFile ( this : in out dataSet;
                             inputFile : in Link_to_String;
                             fail : out boolean ) is

      file : file_type;

    begin
      fail := false;
      declare
      begin
        Open(file,in_file,inputFile.all);
      exception
        when others => fail := true; return;
      end;
      parse_preamble(file,this,fail);
      if not fail then
        parse_supports(file,this,fail);
        if not fail
         then this.outFile := makeOutputFile(inputFile);
        end if;
      end if;
    end getInputFile;

    procedure info_preamble ( this : in dataSet ) is
    begin
      put("Dim = "); put(this.dim,1); new_line;
      put("Support = "); put(this.supN,1); new_line; new_line;
      put("Elem = ");
      for i in 0..this.supN-1 loop
        put(this.termSet(i),1); put(" ");
      end loop;
      new_line;
      put("Type = ");
      for i in 0..this.supN-1 loop
        put(this.supType(i),1); put(" ");
      end loop;
      new_line;
    end info_preamble;

    procedure info_supports ( this : in dataSet ) is

      counter : integer32 := 0;
      top : integer32 := 0;

    begin
      for k in 0..this.supN-1 loop
        for j in top..this.termSet(k) + top - 1 loop
          for i in 0..this.dim-1 loop
            put(integer32(support_out(this,j,i)),1); put(" ");
          end loop;
          new_line;
          counter := counter + 1;
        end loop;
        new_line;
        top := counter;
      end loop;
    end info_supports;

  end class_dataSet;

end demics_input_data;
