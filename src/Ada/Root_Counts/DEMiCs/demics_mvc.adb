with Ada.text_io;                       use Ada.text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with DEMiCs_Global_Constants;

package body demics_mvc is

  package body class_mvc is

    procedure getMemory
                ( this : in Link_to_mvc;
                  depth : in integer32; lvl : in integer32;
                  length : in integer32; vrblvl : in integer32 := 0 ) is

      elemLen,div : integer32;

    begin
      if vrblvl > 0
       then put("-> in demics_mvc.getMemory, ");
      end if;
      if depth /= 0
       then div := 2;
       else div := 1;
      end if;
      if lvl /= length-1 
       then elemLen := this.termSet(depth)/((lvl + 1)*div);
       else elemLen := 1;
      end if;
      if vrblvl > 0 then
        put("depth : "); put(depth,1);
        put(", lvl : "); put(lvl,1);
        put(", elemLen : "); put(elemLen,1); new_line;
      end if;
      this.lv(depth).fTest(lvl)
        := new demics_fTest.class_ftData.ftData'
              (demics_fTest.class_ftData.new_ftData);
      for i in 0..this.termSet(depth)-1 loop
        demics_fTest.class_ftData.create_elem
          (this.lv(depth).fTest(lvl),this.row,this.col,
           this.termSet(depth),this.supType(depth),vrblvl-1);
        demics_fTest.class_ftData.add_elem(this.lv(depth).fTest(lvl),vrblvl-1);
      end loop;
      demics_fTest.class_ftData.mark(this.lv(depth).fTest(lvl),vrblvl-1);
      this.lv(depth).fTest(lvl).cur := this.lv(depth).fTest(lvl).head;
    end getMemory;

    procedure initMemoryCheck
                ( this : in Link_to_mvc;
                  data : in demics_fTest.class_ftData.Link_to_ftData;
                  depth : in integer32 ) is

      sn : constant integer32 := this.sp(depth);

      use demics_fTest.class_theData;

    begin
      if data.cur = null then
        demics_fTest.class_ftData.create_elem
          (data,this.row,this.col,this.termSet(sn),this.supType(sn));
        demics_fTest.class_ftData.add_elem(data);
      end if;
    end initMemoryCheck;

    procedure memoryCheck
                ( this : in Link_to_mvc;
                  data : in demics_fTest.class_ftData.Link_to_ftData;
                  depth : in integer32 ) is
    begin
      initMemoryCheck(this,data,depth); -- same code as initMemoryCheck
    end memoryCheck;

    procedure get_candIdx
                ( this : in Link_to_mvc;
                  curInif : in demics_itest.class_inifData.Link_to_inifData
                ) is

      num : integer32 := 0;
      n_curr : demics_iTest.class_uData.Link_to_uData := curInif.fHead;

      use demics_iTest.class_uData;

    begin
      while n_curr /= null loop
        this.candIdx(num+1) := n_curr.supLab;
        n_curr := n_curr.fNext;
        num := num + 1;
      end loop;
      this.candIdx(0) := num;
    end get_candIdx;

    function chooseSup
                ( this : Link_to_mvc; depth : integer32;
                  curNode : demics_fTest.class_theData.Link_to_theData;
             curInif : demics_itest.class_inifData.Link_to_Array_of_inifData;
             nextInif : demics_itest.class_inifData.Link_to_Array_of_inifData;
                  vrblvl : integer32 := 0 )
                return integer32 is

      flag : integer32;

      use Standard_Floating_Vectors;
      use demics_fTest.class_theData;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.chooseSup, depth : ");
        put(depth,1); put_line(" ...");
        if curNode = null then
          put_line("curNode = null, BUG!");
        else
          put_line("curNode /= null, okay");
          if curNode.transMat = null
           then put_line("curNode.transMat = null");
           else put_line("curNode.transMat /= null");
          end if;
          if curNode.transMat_ptr = null
           then put_line("curNode.transMat_ptr = null");
           else put_line("curNode.transMat_ptr /= null");
          end if;
        end if;
      end if;
      case depth is
        when 0 => fUpdateDirRed(this,curInif,nextInif,curNode,
                                this.iLv(depth).rsp,depth,vrblvl-1);
        when others => updateDirRed(this,curInif,nextInif,curNode,
                                    this.iLv(depth).rsp,depth,vrblvl-1);
      end case;
      case curNode.artV is
        when 0 =>
          flag := findUnbDir(this,nextInif,curNode,this.iLv(depth+1).rsp,
                             this.iLv(depth).rsp,depth,vrblvl-1);    
        when 1 =>
          flag := findUnbDir_art(this,nextInif,curNode,this.iLv(depth+1).rsp,
                                 this.iLv(depth).rsp,depth,vrblvl-1);
        when others => null;
      end case;
      return flag;
    end chooseSup;

    procedure fUpdateDirRed
                ( this : in Link_to_mvc;
            curInif : in demics_itest.class_inifData.Link_to_Array_of_inifData;
           nextInif : in demics_itest.class_inifData.Link_to_Array_of_inifData;
                  curNode : in demics_fTest.Class_theData.Link_to_theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32; vrblvl : in integer32 := 0 ) is

      num,nfPos,idx,fIdx,nfN,flag : integer32;
      length,lvl,pivOutNum,colPos,rowPos : integer32;
      nf_pos,pivOutList : Standard_Integer_Vectors.Link_to_Vector;
      val,preRed : double_float;
      transRed : Standard_Floating_Vectors.Link_to_Vector;
      c_curr,n_curr : demics_iTest.class_uData.Link_to_uData;

      use Standard_Floating_Vectors;
      use demics_iTest.class_uData;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.fUpdateDirRed, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      if curNode.transRed_ptr /= null
       then transRed := curNode.transRed_ptr;
       else put_line("applying patch ...");
            transRed := curNode.transRed;
      end if;
      nf_pos := curNode.nf_pos_ptr;
      nfN := curNode.nfN;
      pivOutNum := curNode.pivOutNum;
      pivOutList := curNode.pivOutList;
      length := this.supN - depth - 1;
      fIdx := this.firIdx(depth);
      colPos := this.termStart(this.sp(depth));
      if curNode.transMat_ptr /= null then
        for i in 0..this.dim*this.dim-1 loop
         -- put("accessing at i = "); put(i,1); new_line;
          this.trMat(i) := curNode.transMat_ptr(i);
        end loop;
      else
        put_line("applying patch ...");
        for i in 0..this.dim*this.dim-1 loop
         -- put("accessing at i = "); put(i,1); new_line;
          this.trMat(i) := curNode.transMat(i);
        end loop;
      end if;
      for j in 0..this.dim-1 loop
        this.trMat(j + j*this.dim) := this.trMat(j + j*this.dim) - 1.0;
        for i in 0..this.dim-1 loop
          this.trMat(i + j*this.dim)
            := this.trMat(i + j*this.dim)*double_float(this.trNeg(fIdx)(i));
        end loop;
      end loop;
      for j in 0..length-1 loop
        lvl := curRsp(j);
        rowPos := this.termStart(lvl);
        c_curr := curInif(lvl).fHead;
        n_curr := nextInif(lvl).fHead;
        num := 0;
        while c_curr /= null loop
          flag := DEMiCs_Global_Constants.CONTINUE;
          for i in 0..curNode.polyDim loop
            if table_out(this,colPos + curNode.nodeLabel(i),
                              rowPos + c_curr.supLab)
                    = DEMiCs_Global_Constants.UNBOUNDED then
              flag := DEMiCs_Global_Constants.UNBOUNDED; exit;
            end if; 
          end loop;
          if flag = DEMiCs_Global_Constants.CONTINUE then
            n_curr.supLab := c_curr.supLab;
           -- for direction
            for k in 0..nfN-1 loop
              val := 0.0;
              nfPos := nf_pos(k);
              for i in 0..pivOutNum-1 loop
                idx := pivOutList(i);
                val := val + this.trMat(idx + nfPos*this.dim)*c_curr.dir(idx);
              end loop;
              n_curr.dir(nfPos) := val
                + double_float(this.trNeg(fIdx)(nfPos))*c_curr.dir(nfPos);
            end loop; 
           --  for reduced cost
            val := 0.0;
            preRed := 0.0;
            for i in 0..this.dim - 1 loop
              val := val
                - double_float(this.trNeg(fIdx)(i))*transRed(i)*c_curr.dir(i);
              preRed := preRed
                + double_float(this.trNeg(fIdx)(i))*c_curr.dir(i);
            end loop;
            n_curr.red := val - preRed + c_curr.red;
          else
           -- skipPtr(this,n_curr,nextInif(lvl).fHead,vrblvl-1);
            skipPtr(n_curr,nextInif(lvl).fHead,vrblvl-1);
          end if;
          c_curr := c_curr.fNext;
          n_curr := n_curr.fNext;
          num := num + 1;
        end loop;
        if n_curr /= null
         then n_curr.prev.fNext := null;
        end if;
      end loop;
    end fUpdateDirRed;

    procedure updateDirRed
                ( this : in Link_to_mvc;
            curInif : in demics_itest.class_inifData.Link_to_Array_of_inifData;
           nextInif : in demics_itest.class_inifData.Link_to_Array_of_inifData;
                  curNode : in demics_fTest.class_theData.Link_to_theData;
                  curRsp : in Standard_Integer_Vectors.Link_to_Vector;
                  depth : in integer32; vrblvl : in integer32 := 0 ) is

      num,nfPos,idx,nfN,flag : integer32;
      length,pivOutNum,lvl,colPos,rowPos : integer32;
      nf_pos,pivOutList : Standard_Integer_Vectors.Link_to_Vector;
      val : double_float;
      transRed : Standard_Floating_Vectors.Link_to_Vector;
      c_curr,n_curr : demics_iTest.class_uData.Link_to_uData;

      use demics_iTest.class_uData;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.updateDirRed, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      transRed := curNode.transRed_ptr;
      nf_pos := curNode.nf_pos_ptr;
      nfN := curNode.nfN;
      pivOutNum := curNode.pivOutNum;
      pivOutList := curNode.pivOutList;
      length := this.supN - depth - 1;
      colPos := this.termStart(this.sp(depth));
      for i in 0..this.dim*this.dim-1 loop
        this.trMat(i) := curNode.transMat_ptr(i);
      end loop;
      for i in 0..this.dim-1 loop
        this.trMat(i + i * this.dim) := this.trMat(i + i * this.dim) - 1.0;
      end loop;
      for j in 0..length-1 loop
        lvl := curRsp(j);
        rowPos := this.termStart(lvl);
        c_curr := curInif(lvl).fHead;
        n_curr := nextInif(lvl).fHead;
        num := 0;
        while c_curr /= null loop
          flag := DEMiCs_Global_Constants.CONTINUE;
          for i in 0..curNode.polyDim loop
            if table_out(this,colPos + curNode.nodeLabel(i), 
                              rowPos + c_curr.supLab)
                    = DEMiCs_Global_Constants.UNBOUNDED then
              flag := DEMiCs_Global_Constants.UNBOUNDED; exit;
            end if; 
          end loop;
          if flag = DEMiCs_Global_Constants.CONTINUE then
            n_curr.supLab := c_curr.supLab;
           -- for direction
            for k in 0..nfN-1 loop
              val := 0.0;
              nfPos := nf_pos(k);
              for i in 0..pivOutNum-1 loop
                idx := pivOutList(i);
                val := val + this.trMat(idx + nfPos*this.dim)*c_curr.dir(idx);
              end loop;
              n_curr.dir(nfPos) := val + c_curr.dir(nfPos);
            end loop;
           -- for reduced cost
            val := 0.0;
            for i in 0..pivOutNum-1 loop
              idx := pivOutList(i);
              val := val - transRed(idx)*c_curr.dir(idx);
            end loop;
            n_curr.red := val + c_curr.red;
          else
           -- skipPtr(this,n_curr,nextInif(lvl).fHead,vrblvl-1);
            skipPtr(n_curr,nextInif(lvl).fHead,vrblvl-1);
          end if;
          c_curr := c_curr.fNext;
          n_curr := n_curr.fNext;
          num := num + 1;
        end loop;
        if n_curr /= null
         then n_curr.prev.fNext := null;
        end if;
      end loop;
    end updateDirRed;

    function findUnbDir
               ( this : Link_to_mvc;
             nextInif : demics_itest.class_inifData.Link_to_Array_of_inifData;
                 curNode : demics_fTest.class_theData.Link_to_theData;
                 nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                 curRsp : Standard_Integer_Vectors.Link_to_Vector;
                 depth : integer32; vrblvl : integer32 := 0 )
               return integer32 is

     nfN,flag,lvl,length,cnt,feasNum : integer32;
     min_feasNum : integer32 := DEMiCs_Global_Constants.BIGINT;
     min_lvl : integer32 := 0;
    -- basisIdx : Standard_Integer_Vectors.Link_to_Vector;
     nf_pos : Standard_Integer_Vectors.Link_to_Vector;
     n_curr,fHead,cor_ptr : demics_iTest.class_uData.Link_to_uData;

     use demics_iTest.class_uData;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.findUnbDir, depth : ");
        put(depth,1); put_line(" ...");
      end if;
     -- basisIdx := curNode.basisIdx_ptr;
      nf_pos := curNode.nf_pos_ptr;
      nfN := curNode.nfN;
      length := this.supN - depth - 1;
      for i in 0..length-1 loop
        lvl := curRsp(i);
        if vrblvl > 0 then -- #if DBG_FINDUNB
          put("-------- Support : "); put(lvl+1,1); put_line(" --------");
        end if;
        fHead := nextInif(lvl).fHead;
        n_curr := fHead;
        feasNum := 0;
        while n_curr /= null loop
          if vrblvl > 0 then -- #if DBG_FINDUNB
            put("-- tarIdx "); put(n_curr.supLab + 1,1); put_line(" --");
          end if;
          cor_ptr := nextInif(lvl).fHead;
         -- checkDir(this,cor_ptr,n_curr,n_curr.dir,n_curr.red,nf_pos,
         --          basisIdx,nfN,flag,vrblvl-1);
          checkDir(cor_ptr,n_curr,n_curr.dir,n_curr.red,nf_pos,
                   nfN,flag,vrblvl-1);
          if flag = DEMiCs_Global_Constants.UNB_TAR then
	   -- skipPtr(this,n_curr,nextInif(lvl).fHead,vrblvl-1);
	    skipPtr(n_curr,nextInif(lvl).fHead,vrblvl-1);
            if vrblvl > 0 then -- #if DBG_FINDUNB
              put_line("UNB_TAR");
            end if;
          elsif flag = DEMiCs_Global_Constants.UNB_COR then
           -- skipPtr(this,cor_ptr,nextInif(lvl).fHead,vrblvl-1);
            skipPtr(cor_ptr,nextInif(lvl).fHead,vrblvl-1);
            feasNum := feasNum + 1;
            if vrblvl > 0 then -- #if DBG_FINDUNB
              put_line("UNB_COR");
            end if;
          else
            feasNum := feasNum + 1;
            if vrblvl > 0 then -- #if DBG_FINDUNB
              put_line("CONTINUE");
            end if;
          end if;
          n_curr := n_curr.fNext;
        end loop;
        if feasNum < min_feasNum then
          min_feasNum := feasNum;
          min_lvl := lvl;
        end if;
      end loop;
      if vrblvl > 0 then -- #if DBG_FINDUNB
        put("min_lvl : "); put(min_lvl+1,1); new_line;
        put("min_feasNum : "); put(min_feasNum,1); new_line;
      end if;
      this.sp(depth + 1) := min_lvl;
      cnt := 0;
      for i in 0..length-1 loop
        if curRsp(i) /= min_lvl then
          nextRsp(cnt) := curRsp(i);
          cnt := cnt + 1;
        end if;
      end loop;
      if min_feasNum <= 1
       then return DEMiCs_Global_Constants.STOP;
       else return DEMiCs_Global_Constants.CONTINUE;
      end if;
    end findUnbDir;

    function findUnbDir_art
               ( this : Link_to_mvc;
             nextInif : demics_itest.class_inifData.Link_to_Array_of_inifData;
                 curNode : demics_fTest.class_theData.Link_to_theData;
                 nextRsp : Standard_Integer_Vectors.Link_to_Vector;
                 curRsp : Standard_Integer_Vectors.Link_to_Vector;
                 depth : integer32; vrblvl : integer32 := 0 )
               return integer32 is

      nfN,flag,lvl,length,cnt,feasNum : integer32;
      min_feasNum : integer32 := DEMiCs_Global_Constants.BIGINT;
      min_lvl : integer32 := 0;
      basisIdx,nf_pos : Standard_Integer_Vectors.Link_to_Vector;
      n_curr,fHead,cor_ptr : demics_iTest.class_uData.Link_to_uData;

      use demics_iTest.class_uData;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.findUnbDir_art, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      basisIdx := curNode.basisIdx_ptr;
      nf_pos := curNode.nf_pos_ptr;
      nfN := curNode.nfN;
      length := this.supN - depth - 1;
      for i in 0..length-1 loop
        lvl := curRsp(i);
        if vrblvl > 0 then -- #if DBG_FINDUNB
          put("-------- Support : "); put(lvl+1,1); put_line(" --------");
        end if;
        fHead := nextInif(lvl).fHead;
        n_curr := fHead;
        feasNum := 0;
        while n_curr /= null loop
          if vrblvl > 0 then -- #if DBG_FINDUNB
            put("-- tarIdx "); put(n_curr.supLab+1,1); put_line(" --");
          end if;
          cor_ptr := nextInif(lvl).fHead;
          checkDir_art(this,cor_ptr,n_curr,n_curr.dir,n_curr.red, 
                       nf_pos,basisIdx,nfN,flag,vrblvl-1);
          if flag = DEMiCs_Global_Constants.UNB_TAR then
           -- skipPtr(this,n_curr,nextInif(lvl).fHead,vrblvl-1);
            skipPtr(n_curr,nextInif(lvl).fHead,vrblvl-1);
            if vrblvl > 0 then -- #if DBG_FINDUNB
              put_line("UNB");
            end if;
          elsif flag = DEMiCs_Global_Constants.UNB_COR then
           -- skipPtr(this,cor_ptr,nextInif(lvl).fHead,vrblvl-1);
            skipPtr(cor_ptr,nextInif(lvl).fHead,vrblvl-1);
            feasNum := feasNum + 1;
          else
            feasNum := feasNum + 1;
            if vrblvl > 0 then -- #if DBG_FINDUNB
              put_line("CONTINUE");
            end if;
          end if;
          n_curr := n_curr.fNext;
        end loop; 
        if feasNum < min_feasNum then
          min_feasNum := feasNum;
          min_lvl := lvl;
        end if;
      end loop;
      if vrblvl > 0 then -- #if DBG_FINDUNB
        put("min_lvl : "); put(min_lvl+1,1); new_line;
        put("min_feasNum : "); put(min_feasNum,1); new_line;
      end if;
      this.sp(depth+1) := min_lvl;
      cnt := 0;
      for i in 0..length-1 loop
        if curRsp(i) /= min_lvl then
          nextRsp(cnt) := curRsp(i);
          cnt := cnt + 1;
        end if;
      end loop;
      if min_feasNum <= 1
       then return DEMiCs_Global_Constants.STOP;
       else return DEMiCs_Global_Constants.CONTINUE;
      end if;
    end findUnbDir_art;

    procedure checkDir
               -- ( this : in Link_to_mvc;
                ( corPtr : in out demics_iTest.class_uData.Link_to_uData;
                  tarPtr : in demics_iTest.class_uData.Link_to_uData;
                  tar_dir : in Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : in double_float;
                  nf_pos : in Standard_Integer_Vectors.Link_to_Vector;
                 -- basisIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  nfN : in integer32; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      nfPos,ans,sign : integer32;

      use demics_iTest.class_uData;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_mvc.class_mvc.checkDir ...");
      end if;
      while corPtr /= null loop
        if corPtr /= tarPtr then
          if vrblvl > 0 then -- #if DBG_FINDUNB
            put(corPtr.supLab+1,1); put(" : ");
            put(corPtr.red - tar_red); put(" : ");
          end if;
         -- sign := checkSign_red(this,corPtr.red,tar_red);
          sign := checkSign_red(corPtr.red,tar_red);
          if sign = DEMiCs_Global_Constants.NEGATIVE then
            flag := DEMiCs_Global_Constants.UNBOUNDED;
            for i in 0..nfN-1 loop
              nfPos := nf_pos(i);
              if vrblvl > 0 then -- #if DBG_FINDUNB
                put(corPtr.dir(nfPos) - tar_dir(nfPos)); put(" ");
              end if;
             -- ans := checkNonNeg_dir(this,corPtr.dir(nfPos),tar_dir(nfPos));
              ans := checkNonNeg_dir(corPtr.dir(nfPos),tar_dir(nfPos));
              if ans = DEMiCs_Global_Constants.FALSE then
                if vrblvl > 0 then -- #if DBG_FINDUNB
                  new_line;
                end if;
                flag := DEMiCs_Global_Constants.CONTINUE;
                exit;
              end if;
            end loop;
            if flag = DEMiCs_Global_Constants.UNBOUNDED then
              if vrblvl > 0 then -- #if DBG_FINDUNB	  
                new_line;
              end if;
              flag := DEMiCs_Global_Constants.UNB_TAR;
              return;
            end if;
          else
            flag := DEMiCs_Global_Constants.UNBOUNDED;
            for i in 0..nfN-1 loop
              nfPos := nf_pos(i);
              if vrblvl > 0 then -- #if DBG_FINDUNB
                put(corPtr.dir(nfPos) - tar_dir(nfPos)); put(" ");
              end if;
             -- ans := checkNonPos_dir(this,corPtr.dir(nfPos),tar_dir(nfPos));
              ans := checkNonPos_dir(corPtr.dir(nfPos),tar_dir(nfPos));
              if ans = DEMiCs_Global_Constants.FALSE then
                if vrblvl > 0 then -- #if DBG_FINDUNB
                  new_line;
                end if;
                flag := DEMiCs_Global_Constants.CONTINUE;
                exit;
              end if;
            end loop;
            if flag = DEMiCs_Global_Constants.UNBOUNDED then
              if vrblvl > 0 then -- #if DBG_FINDUNB	  
                new_line;
              end if;
              flag := DEMiCs_Global_Constants.UNB_COR;
              return;
            end if;
          end if;
        end if;
        corPtr := corPtr.fNext;
      end loop;
      if vrblvl > 0 then -- #if DBG_FINDUNB	  
        new_line;
      end if;
      flag := DEMiCs_Global_Constants.CONTINUE;
    end checkDir;

    procedure checkDir_art
                ( this : in Link_to_mvc;
                  corPtr : in out demics_iTest.class_uData.Link_to_uData;
                  tarPtr : in demics_itest.class_uData.Link_to_uData;
                  tar_dir : in Standard_Floating_Vectors.Link_to_Vector;
                  tar_red : in double_float;
                  nf_pos : in Standard_Integer_Vectors.Link_to_Vector;
                  basisIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  nfN : in integer32; flag : out integer32;
                  vrblvl : in integer32 := 0 ) is

      nfPos,ans,sign,nonNegVarNum : integer32;
     -- cnt : integer32;

      use demics_iTest.class_uData;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_mvc.class_mvc.checkbDir_art ...");
      end if;
      while corPtr /= null loop
        if corPtr /= tarPtr then
          if vrblvl > 0 then -- #if DBG_FINDUNB
            put(corPtr.supLab+1,1); put(" : ");
          end if;
          flag := DEMiCs_Global_Constants.CONTINUE;
          for i in 0..nfN-1 loop
            nfPos := nf_pos(i);
            if basisIdx(nfPos) >= this.termSumNum - this.supN  then
             -- ans := checkZero_dir(this,corPtr.dir(nfPos),tar_dir(nfPos));
              ans := checkZero_dir(corPtr.dir(nfPos),tar_dir(nfPos));
              if ans = DEMiCs_Global_Constants.FALSE then
                flag := DEMiCs_Global_Constants.STOP;
                exit;
              end if;
            end if;
          end loop;
          exit when (flag = DEMiCs_Global_Constants.STOP);
         -- sign := checkSign_red(this,corPtr.red,tar_red);
          sign := checkSign_red(corPtr.red,tar_red);
          if(sign = DEMiCs_Global_Constants.NEGATIVE) then
           -- cnt := 0;
            nonNegVarNum := 0;
            flag := DEMiCs_Global_Constants.STOP;
            for i in 0..nfN-1 loop
              nfPos := nf_pos(i);
              if vrblvl > 0 then -- #if DBG_FINDUNB
                put(corPtr.dir(nfPos) - tar_dir(nfPos)); put(" ");
              end if;
              if basisIdx(nfPos) < this.termSumNum - this.supN then
                nonNegVarNum := nonNegVarNum + 1;
               -- ans := checkNonNeg_dir
               --          (this,corPtr.dir(nfPos),tar_dir(nfPos));
                ans := checkNonNeg_dir(corPtr.dir(nfPos),tar_dir(nfPos));
                if ans = DEMiCs_Global_Constants.FALSE then
                  if vrblvl > 0 then -- #if DBG_FINDUNB
                    new_line;
                  end if;
                  flag := DEMiCs_Global_Constants.STOP; exit;
                else
                  flag := DEMiCs_Global_Constants.UNBOUNDED;
                end if;
              end if;
            end loop;
            if flag = DEMiCs_Global_Constants.UNBOUNDED then
              if vrblvl > 0 then -- #if DBG_FINDUNB	  
                new_line;
              end if;
              flag := DEMiCs_Global_Constants.UNB_TAR;
              return;
            end if;
          else
           -- cnt := 0;
            nonNegVarNum := 0;
            flag := DEMiCs_Global_Constants.STOP;
            for i in 0..nfN-1 loop
              nfPos := nf_pos(i);
              if vrblvl > 0 then -- #if DBG_FINDUNB
                put(corPtr.dir(nfPos) - tar_dir(nfPos)); put(" ");
              end if;
              if basisIdx(nfPos) < this.termSumNum - this.supN then
                nonNegVarNum := nonNegVarNum + 1;
               -- ans := checkNonPos_dir
               --          (this,corPtr.dir(nfPos),tar_dir(nfPos));
                ans := checkNonPos_dir(corPtr.dir(nfPos),tar_dir(nfPos));
                if ans = DEMiCs_Global_Constants.FALSE then
                  if vrblvl > 0 then -- #if DBG_FINDUNB
                    new_line;
                  end if;
                  flag := DEMiCs_Global_Constants.STOP; exit;
                else
                  flag := DEMiCs_Global_Constants.UNBOUNDED;
                end if;
              end if;
            end loop;
            if flag = DEMiCs_Global_Constants.UNBOUNDED then
              if vrblvl > 0 then -- #if DBG_FINDUNB	  
                new_line;
              end if;
              flag := DEMiCs_Global_Constants.UNB_COR;
              return;
            end if;
          end if;
        end if;
        corPtr := corPtr.fNext;
      end loop;
      if vrblvl > 0 then -- #if DBG_FINDUNB	  
        new_line;
      end if;
      flag := DEMiCs_Global_Constants.CONTINUE;
    end checkDir_art;

    procedure skipPtr
               -- ( this : in Link_to_mvc;
                ( curr : in demics_itest.class_uData.Link_to_uData;
                  fHead : in out demics_itest.class_uData.Link_to_uData;
                  vrblvl : in integer32 := 0 ) is

      use demics_iTest.class_uData;

    begin
      if vrblvl > 0
       then put_line("-> in demics_mvc.class_mvc.skipPtr ...");
      end if;
      if curr = fHead then
        fHead := curr.fNext;
      elsif curr.fNext /= null then
        curr.prev.fNext := curr.fNext;
        curr.fNext.prev := curr.prev;
      else
        curr.prev.fNext := curr.fNext;
      end if;
    end skipPtr;

    procedure get_tuple_index
               -- ( this : in Link_to_mvc;
                ( node : in demics_fTest.class_ftData.Link_to_ftData;
                  data : in demics_fTest.class_ftData.Link_to_Array_of_ftData;
                  length : in integer32; vrblvl : in integer32 := 0 ) is

      on : constant := DEMiCs_Global_Constants.ON;
      tmpIdx : integer32;

      use Standard_Integer_Vectors;
      use demics_fTest.class_theData;
      use demics_fTest.class_ftData;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.get_tuple_index, length : ");
        put(length,1); new_line;
        for i in 0..length-2 loop
          if data = null
           then put("data = null "); put_line("BUG!");
           else put_line("data /= null");
          end if;
          put("data("); put(i); put(")");
          if data(i) = null
           then put_line(" = null BUG!");
           else put_line(" /= null");
          end if;
          put("data("); put(i); put(").parent");
          if data(i).parent = null
           then put_line(" = null BUG!");
           else put_line(" /= null");
          end if;
          if node = null
           then put_line("node = null BUG!");
           else put_line("node /= null");
          end if;
          if node.parent = null
           then put_line("node.parent = null BUG!");
           else put_line("node.parent /= null");
          end if;
          if node.parent.nodeLabel = null
           then put_line("node.parent.nodeLabel = null BUG!");
           else put_line("node.parent.nodeLabel /= null");
          end if;
        end loop;
      end if;
      for i in 0..length-2 loop
        node.parent.nodeLabel(i) := data(i).parent.fIdx;
      end loop;
      node.parent.nodeLabel(length-1) := data(length-1).cur.fIdx;
      if data(1).parent.sw = on then
        tmpIdx := node.parent.nodeLabel(1);
        node.parent.nodeLabel(1) := node.parent.nodeLabel(0);
        node.parent.nodeLabel(0) := tmpIdx;
      end if;
      if vrblvl > 0
       then put_line("exiting get_tuple_index");
      end if;
    end get_tuple_index;

--  procedure dbg_init_transMat
--              ( this : in Link_to_mvc;
--                curNode : in demics_fTest.class_theData.Link_to_theData ) is
--  begin
--    null;
--  end dbg_init_transMat;

--  procedure dbg_transMat
--              ( this : in Link_to_mvc;
--                preNode : in demics_fTest.class_theData.Link_to_theData;
--                curNode : in demics_fTest.class_theData.Link_to_theData ) is
--  begin
--    null;
--  end dbg_transMat;

--  procedure check_transMat
--              ( this : in Link_to_mvc;
--                preNode : in demics_fTest.class_theData.Link_to_theData;
--                curNode : in demics_fTest.class_theData.Link_to_theData ) is
--  begin
--    null;
--  end check_transMat;

--  procedure check_init_transRed
--              ( this : in Link_to_mvc;
--                curNode : in demics_fTest.class_theData.Link_to_theData ) is
--  begin
--    null;
--  end check_init_transRed;

    function checkSign_red
               -- ( this : Link_to_mvc;
                ( curRed : double_float;
                  tarRed : double_float ) return integer32 is

      flag : integer32 := 0;
      pluszero : constant := DEMiCs_Global_Constants.PLUSZERO;
      minuszero : constant := DEMiCs_Global_Constants.MINUSZERO;
      positive : constant := DEMiCs_Global_Constants.POSITIVE;
      negative : constant := DEMiCs_Global_Constants.NEGATIVE;

    begin
      if curRed > tarRed + pluszero then
        flag := positive;
      elsif curRed < tarRed + minuszero then
        flag := negative;
      end if;
      return flag;
    end checkSign_red;

    function checkNonNeg_dir
               -- ( this : Link_to_mvc;
                ( curDirElem : double_float;
                  tarDirElem : double_float ) return integer32 is
    begin
      if curDirElem < tarDirElem + DEMiCs_Global_Constants.PLUSZERO
       then return DEMiCs_Global_Constants.TRUE;
       else return DEMiCs_Global_Constants.FALSE;
      end if;
    end checkNonNeg_dir;

    function checkNonPos_dir
              -- ( this : Link_to_mvc;
               ( curDirElem : double_float;
                 tarDirElem : double_float ) return integer32 is
    begin
      if curDirElem > tarDirElem + DEMiCs_Global_Constants.MINUSZERO
       then return DEMiCs_Global_Constants.TRUE;
       else return DEMiCs_Global_Constants.FALSE;
      end if;
    end checkNonPos_dir;

    function checkZero_dir
               -- ( this : Link_to_mvc;
                ( curDirElem : double_float;
                  tarDirElem : double_float ) return integer32 is

      val : constant double_float := curDirElem - tarDirElem;

    begin
      if (DEMiCs_Global_Constants.MINUSZERO < val) and 
         (val < DEMiCs_Global_Constants.PLUSZERO)
       then return DEMiCs_Global_Constants.TRUE;
       else return DEMiCs_Global_Constants.FALSE;
      end if;
    end checkZero_dir;

    function table_out
                ( this : in Link_to_mvc;
                  row : integer32; col : integer32 ) return integer32 is
    begin
      return this.table(row + col*this.termSumNum);
    end table_out;

    procedure info_neg
                ( this : in Link_to_mvc; termSet : in integer32;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec ) is
    begin
      put_line("<< trNeg >>");
      for j in 0..termSet-1 loop
        for i in 0..this.row-1 loop
          put(this.trNeg(j)(i),1); put(" ");
        end loop;
        new_line;
      end loop;
      put_line("<< negIdx >>");
      for j in 0..termSet-1 loop
        for i in 0..negIdx(j)(0) loop
          put(negIdx(j)(i+1),1); put(" ");
        end loop;
        new_line;
      end loop;
    end info_neg;

    procedure info_sp ( this : in Link_to_mvc; depth : in integer32 ) is
    begin
      put("sp :");
      for i in 0..depth-1 loop
        put(" "); put(this.sp(i),1);
      end loop;
      new_line;
    end info_sp;

    procedure info_parent_node ( this : in Link_to_mvc;
                                 depth : in integer32 ) is
    begin
      put("Node : ");
      for i in 0..depth-1 loop
        put(this.sp(i),1); put(" : ");
        demics_fTest.class_ftData.info_parent_node(this.lv(this.sp(i)).node);
      end loop;
      new_line;
    end info_parent_node;

    procedure info_tuple
                ( this : in Link_to_mvc; lvl : in integer32 ) is
               -- ; depth : in integer32 ) is
    begin
      put("( ");
      for i in 0..lvl-1 loop
        put(this.mFeaIdx(i)(this.mRepN(i))+1,1); put(" ");
      end loop;
      put_line(")");
    end info_tuple;

    procedure info_all_dirRed
                ( this : in Link_to_mvc;
                  depth : in integer32;
                  node : in demics_fTest.class_ftData.Link_to_ftData;
                  nextInif : in demics_itest.class_inifData.Array_of_inifData
                ) is

      nfPos,nfN : integer32;
      nf_pos : Standard_Integer_Vectors.Link_to_Vector;
      val : double_float;
      n_curr : demics_iTest.class_uData.Link_to_uData;

      use demics_iTest.class_uData;

    begin
      put_line("<< info_all_dirRed >>");
      nf_pos := node.parent.nf_pos_ptr;
      nfN := node.parent.nfN;
      for j in depth+1..this.supN-1 loop
         n_curr := nextInif(j).fHead;
         put("--- Support : "); put(j+1,1); put_line(" ---");
         while n_curr /= null loop
            put(n_curr.supLab+1,1); put(" : ");
            for k in 0..nfN-1 loop
              nfPos := nf_pos(k);
              val := n_curr.dir(nfPos);
              if (val < DEMiCs_Global_Constants.PLUSZERO) and
                 (val > DEMiCs_Global_Constants.MINUSZERO) then
                put("0 ");
              else
                put(val); put(" ");
              end if;
            end loop;
            put(" : "); put(n_curr.red);
            n_curr := n_curr.fNext;
            new_line;
         end loop;
         new_line;
      end loop;
    end info_all_dirRed;

    procedure info_mFea ( this : in Link_to_mvc; length : in integer32 ) is
    begin
      put("mFea :");
      for i in 0..length-1 loop
        put(" "); put(this.mFea(i),1);
      end loop;
      new_line;
      put("mRepN :");
      for i in 0..length-1 loop
        put(" "); put(this.mRepN(i),1);
      end loop;
      new_line;
    end info_mFea;

    procedure info_firIdx ( this : in Link_to_mvc; length : in integer32 ) is
    begin
      put_line("<< firIdx >>");
      for i in 0..length loop
        put(this.firIdx(i),1); put(" ");
      end loop;
      new_line;
    end info_firIdx;

    procedure info_fIdx
               -- ( this : in Link_to_mvc;
                ( data : in demics_fTest.class_ftData.Link_to_ftData ) is
    begin
      put("First Index : "); put(data.parent.fIdx+1,1); new_line;
    end info_fIdx;

    procedure info_candIdx ( this : in Link_to_mvc ) is
    begin
      put("candIdx :");
      for i in 0..this.candIdx(0)-1 loop
        put(this.candIdx(i+1),1); put(" ");
      end loop;
      new_line;
    end info_candIdx;

    procedure info_elemNum
               -- ( this : in Link_to_mvc;
                ( length : in integer32;
                  data : in demics_fTest.class_ftData.Link_to_Array_of_ftData;
                  node : in demics_fTest.class_ftData.Link_to_ftData ) is
    begin
      put("numElem : ");
      for i in 0..length-2 loop
        put(data(i).elemNum); put(" ");
      end loop;
      put(node.elemNum); new_line;
    end info_elemNum;

    procedure info_prop_elemNum
               -- ( this : in Link_to_mvc;
                ( length : in integer32;
                  data : in demics_fTest.class_ftData.Link_to_Array_of_ftData;
                  node : in demics_fTest.class_ftData.Link_to_ftData ) is
    begin
      put("prop_numElem : ");
      for i in 0..length - 2 loop
        demics_fTest.class_ftData.info_numElem(data(i));
      end loop;
      demics_fTest.class_ftData.info_numElem(node);
      new_line;
    end info_prop_elemNum;

    procedure info_table ( this : in Link_to_mvc ) is

    -- NOTE :
    --   This is essentially the same code as the procedure
    --   demics_reltab.class_reltab.info_table.

      u_cnt : integer32 := 0;
      t_cnt : integer32 := 0;

    begin
      put_line("<< Relation table >>");
      for i in 0..this.termSumNum-1 loop
        for k in 0..i-1 loop
          put("  ");
        end loop;
        for j in i+1..this.termSumNum-1 loop
          if table_out(this,j,i) = DEMiCs_Global_Constants.UNBOUNDED
           then u_cnt := u_cnt + 1;
          end if;
          put(table_out(this,j,i),1); put(" ");
          t_cnt := t_cnt + 1;
        end loop;
        new_line;
      end loop;
      new_line;
      put("# Unb. LPs : "); put(u_cnt,1); new_line;
      put("# Elem. : "); put(t_cnt,1); new_line;
      put("Ratio : "); put(double_float(u_cnt)/double_float(t_cnt));
      new_line;
    end info_table;

    function new_mvc return mvc is
 
      res : mvc;

    begin
      res.dim := 0;
      res.supN := 0;
      res.row := 0;
      res.col := 0;
      res.termSumNum := 0;
      res.termMax := 0;
      res.maxLength := 0;
      res.total_iter := 0.0;
      res.total_feasLP := 0.0;
      res.total_LPs := 0.0;
      res.total_1PT := 0.0;
      res.total_2PT := 0.0;
      res.total_triLPs_mLP := 0.0;
      res.total_unbLP_tab := 0.0;
      res.lvl_1PT := null;
      res.lvl_2PT := null;
      res.actNode := null;
      res.mfNum := null;
      res.termSet := null;
      res.termStart := null;
      res.re_termStart := null;
      res.supType := null;
      res.mRepN := null;
      res.mFeaIdx := null;
      res.mFea := null;
      res.trNeg := null;
      res.firIdx := null;
      res.repN := null;
      res.sp := null;
      res.candIdx := null;
      res.trMat := null;
      res.table := null;
      res.lv := null;
      res.iLv := null;
      return res;
    end new_mvc;

    procedure delete_mvc ( this : in Link_to_mvc ) is

      use Standard_Integer_VecVecs;

    begin
      if this.trNeg /= null then
        for i in 0..this.termSet(this.sp(0))-1 loop
          Standard_Integer_Vectors.clear(this.trNeg(i));
        end loop;
        Standard_Integer_VecVecs.shallow_clear(this.trNeg);
        this.trNeg := null;
      end if;
      Standard_Integer_Vectors.clear(this.re_termStart);
      Standard_Integer_Vectors.clear(this.mfNum);
      Standard_Floating_Vectors.clear(this.lvl_1PT);
      Standard_Floating_Vectors.clear(this.lvl_2PT);
      Standard_Floating_Vectors.clear(this.actNode);
      Standard_Integer_Vectors.clear(this.firIdx);
      Standard_Integer_Vectors.clear(this.repN);
      Standard_Integer_Vectors.clear(this.sp);
      Standard_Integer_Vectors.clear(this.candIdx);
      Standard_Floating_Vectors.clear(this.trMat);
     -- delete [] this.lv;
     -- delete [] this.iLv;
    end delete_mvc;

    procedure allocateAndIni
                ( this : in Link_to_mvc;
                  data : in demics_input_data.class_dataSet.dataSet;
                  seedNum : in integer32; output : in integer32;
                  vrblvl : in integer32 := 0 ) is

      length : integer32;

    begin
      if vrblvl > 0
       then put_line("-> in demics_mvc.allocateAndIni ...");
      end if;
      this.dim := data.dim;
      this.supN := data.supN;
      this.row := data.dim;
      this.termSumNum := data.termSumNum;
      this.termMax := data.termMax;
      this.maxLength := data.typeMax + 1;
      this.col := this.termSumNum - this.supN + this.dim;
      this.termSet := data.termSet;
      this.termStart := data.termStart;
      this.supType := data.supType;
      this.mfNum := new Standard_Integer_Vectors.Vector(0..this.supN-1);
      this.lvl_1PT
        := new Standard_Floating_Vectors.Vector'(0..this.supN-1 => 0.0);
      this.lvl_2PT
        := new Standard_Floating_Vectors.Vector'(0..this.supN-1 => 0.0);
      this.actnode
        := new Standard_Floating_Vectors.Vector'(0..this.supN-1 => 0.0);
      this.firIdx
        := new Standard_Integer_Vectors.Vector'(0..this.supN => 0);
      this.re_termStart := new Standard_Integer_Vectors.Vector(0..this.supN);
      this.re_termStart(0) := 0;
      this.repN
        := new Standard_Integer_Vectors.Vector'(0..this.supN-1 => 0);
      this.sp := new Standard_Integer_Vectors.Vector(0..this.supN-1);
      this.candIdx := new Standard_Integer_Vectors.Vector(0..this.termMax);
      this.trMat
        := new Standard_Floating_Vectors.Vector(0..this.dim*this.dim-1);
      this.lv := new demics_fTest.class_lvData.Array_of_lvData(0..this.supN-1);
      this.iLv
        := new demics_itest.class_iLvData.Array_of_iLvData(0..this.supN-1);
      for i in 0..this.supN-1 loop
        this.re_termStart(i+1) := this.termStart(i+1) - i - 1;
        this.sp(i) := i;
        this.lv(i) := new demics_fTest.class_lvData.lvData'
                         (demics_fTest.class_lvData.new_lvData);
        demics_fTest.class_lvData.create
          (this.lv(i),i,--this.supN,this.dim,
           this.supType(i)+1,this.termMax,vrblvl-1);
        this.iLv(i) := new demics_itest.class_iLvData.ilvData'
                          (demics_itest.class_iLvData.new_iLvData);
        demics_itest.class_iLvData.create
          (this.iLv(i),i,this.supN,this.dim,this.termMax,vrblvl-1);
        length := this.supType(i) + 1;
        for j in 0..length-1 loop
          getMemory(this,i,j,length,vrblvl-1);
        end loop;
      end loop;
      this.the_Simplex := new demics_simplex.class_simplex.simplex'
                             (demics_simplex.class_simplex.new_simplex);
      demics_simplex.class_simplex.allocateAndIni
        (this.the_Simplex,data,this.firIdx,seedNum,output,vrblvl-1);
      this.the_Reltab := new demics_reltab.class_reltab.reltab'
                            (demics_reltab.class_reltab.new_reltab);
      demics_reltab.class_reltab.allocateAndIni
        (this.the_Reltab,this.the_Simplex,this.firIdx,this.dim,this.supN,
         this.termSumNum,this.termSet,this.termStart,this.re_termStart,
         vrblvl-1);
      demics_itest.class_iLvData.getInit
        (this.iLv(0),data,this.the_Simplex.lifting,this.termSet,
         this.termStart,this.dim,this.supN,vrblvl-1);
    end allocateAndIni;

    procedure initFeasTest
                ( this : in Link_to_mvc; depth : in integer32;
                  vrblvl : in integer32 := 0 ) is

      lvl,sn,feaNum,flag : integer32;
     -- length : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.initFeasTest, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      sn := this.sp(depth);
     -- length := this.supType(sn) + 1;
      demics_fTest.class_lvData.get_info
        (this.lv(sn),this.mRepN,this.mFeaIdx,this.mFea);
      lvl := 0;
      initCheck(this,depth,this.lv(sn).fTest(lvl),vrblvl-1);
      lvl := lvl + 1;
      feaNum := 0;
     -- flag := DEMiCs_Global_Constants.CONTINUE;
     -- if flag = DEMiCs_Global_Constants.CONTINUE then
        findNode(this,depth,lvl,feaNum,this.lv(sn).ftest,flag,vrblvl-1);
     -- else
     --   findNextNode(this,depth,lvl,feaNum,this.lv(sn).fTest,flag);
     -- end if;
    end initFeasTest;

    procedure initCheck
                ( this : in Link_to_mvc; depth : in integer32;
                  data : in demics_fTest.class_ftData.Link_to_ftData;
                  vrblvl : in integer32 := 0 ) is

      sn : constant integer32 := this.sp(depth);
      negIdx : Standard_Integer_VecVecs.Link_to_VecVec;
      negNum : integer32;
      feaNum : integer32 := 0;
      val : Standard_Floating_Vectors.Link_to_Vector;
      elem : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.initCheck, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      val := new Standard_Floating_Vectors.Vector(0..this.termSet(sn)-2);
      negIdx := new Standard_Integer_VecVecs.VecVec(0..this.termSet(sn)-1);
      this.trNeg
        := new Standard_Integer_VecVecs.VecVec(0..this.termSet(sn)-1);
      for i in 0..this.termSet(sn)-1 loop
        negIdx(i) := new Standard_Integer_Vectors.Vector(0..this.dim);
        this.trNeg(i) := new Standard_Integer_Vectors.Vector(0..this.dim-1);
      end loop;
      Standard_Random_Numbers.Set_Seed(12);
      for i in 0..this.termSet(sn)-2 loop
        val(i) := abs(Standard_Random_Numbers.Random);
      end loop;
      this.firIdx(this.supN) := 0;
      for idx_one in 0..this.termSet(sn)-1 loop
        if vrblvl > 0 then -- #if DBG_NODE
          put("------------------- Idx : "); put(idx_one+1,1);
          put_line(" -------------------");
        end if;
        initMemoryCheck(this,data,depth);
        demics_fTest.class_theData.clear(data.cur);
        this.firIdx(sn) := idx_one;
        demics_simplex.class_simplex.get_iNbN_nfN
          (this.the_Simplex,data.cur,this.termSet(sn)-1+this.dim,this.dim);
        negNum := 0;
        for i in 0..this.dim-1 loop
          elem := 0.0;
          for j in 0..this.termSet(sn) - 2 loop
            elem := elem + val(j)*demics_simplex.class_simplex.put_elem_supp
                                    (this.the_Simplex,sn,idx_one,i,j);
          end loop;
          if elem < DEMiCs_Global_Constants.MINUSZERO then
            data.cur.p_sol(this.termSumNum - this.supN + i) := -elem;
            negIdx(idx_one)(negNum + 1) := i;
            this.trNeg(idx_one)(i) := -1;
            negNum := negNum + 1;
            for k in 0..this.termSet(sn) - 2 loop
              demics_simplex.class_simplex.mult_elem_supp
                (this.the_Simplex,sn,idx_one,i,k);
            end loop;
          elsif elem > DEMiCs_Global_Constants.PLUSZERO then
            data.cur.p_sol(this.termSumNum - this.supN + i) := elem;
            this.trNeg(idx_one)(i) := 1;
          else
            data.cur.p_sol(this.termSumNum - this.supN + i) := 0.0;
            this.trNeg(idx_one)(i) := 1;
          end if; 
        end loop;
        negIdx(idx_one)(0) := negNum;
        demics_fTest.class_ftData.make_init_data
          (data,this.termSumNum,this.supN,this.termSet(sn),
           this.re_termStart(sn),vrblvl-1);
        initLP(this,data,negIdx,depth,idx_one,feaNum,vrblvl-1);
      end loop;
      Standard_Floating_Vectors.Clear(val);
      for i in 0..this.termSet(sn)-1 loop
        Standard_Integer_Vectors.clear(negIdx(i));
      end loop;
      Standard_Integer_VecVecs.Shallow_Clear(negIdx);
    end initCheck;

    procedure initLP
                ( this : in Link_to_mvc;
                  data : in demics_fTest.class_ftData.Link_to_ftData;
                  negIdx : in Standard_Integer_VecVecs.Link_to_VecVec;
                  depth : in integer32; idx : in integer32;
                  feaNum : in out integer32; vrblvl : in integer32 := 0 ) is

      sn : constant integer32 := this.sp(depth);
      flag,iter : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.initLP, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      demics_simplex.class_simplex.get_cur(this.the_Simplex,data.cur);
      demics_simplex.class_simplex.copy_p1_d_sol(this.the_Simplex,data.cur);
      iter := 0;
      demics_simplex.class_simplex.fSolLP
        (this.the_Simplex,this.termSet(sn),this.re_termStart(sn),
         iter,flag,vrblvl-1);
      this.total_LPs := this.total_LPs + 1.0;
      this.total_1PT := this.total_1PT + 1.0;
      this.lvl_1PT(depth) := this.lvl_1PT(depth) + 1.0;
      if flag = DEMiCs_Global_Constants.OPT then
        if vrblvl > 0 then -- #if DBG_FEA
          put_line("OPT");
        end if;
        this.total_iter := this.total_iter + double_float(iter);
        this.total_feasLP := this.total_feasLP + 1.0;
        demics_fTest.class_theData.joint(data.cur);
        data.cur.fIdx := idx;
        demics_simplex.class_simplex.get_res(this.the_Simplex,data);
        demics_simplex.class_simplex.get_pivOutNum
          (this.the_Simplex,data.cur);
        this.mFeaIdx(depth)(feaNum) := idx;
        this.mFea(depth) := this.mFea(depth) + 1;
        feaNum := feaNum + 1;
        for j in 0..negIdx(idx)(0) - 1 loop
          for k in 0..this.termSet(sn) - 2 loop
            demics_simplex.class_simplex.mult_elem_supp
              (this.the_Simplex,sn,idx,negIdx(idx)(j + 1),k);
          end loop;
          for k in 0..this.dim-1 loop
            data.cur.invB(negIdx(idx)(j + 1) + k*this.dim)
              := -data.cur.invB(negIdx(idx)(j + 1) + k*this.dim);
          end loop; 
          data.cur.d_sol(negIdx(idx)(j + 1)) 
            := -data.cur.d_sol(negIdx(idx)(j + 1));
        end loop;
        if vrblvl > 0 then -- #if DBG_INI_CUR_INFO 
          demics_fTest.class_ftData.info_cur(data);
        end if;
        data.cur := data.cur.next;
      elsif flag = DEMiCs_Global_Constants.UNBOUNDED then
        demics_fTest.class_ftData.init_info(data);
        if vrblvl > 0 then -- #if DBG_FEA
          put_line("UNB");
        end if;
      else
        put_line("Error: too many iterations at initLP");
        put("( "); put(idx,1); put_line(" )");
      end if;
    end initLP;

    function feasTest
                ( this : Link_to_mvc; depth : integer32;
                  parent : demics_fTest.class_theData.Link_to_theData;
                  vrblvl : in integer32 := 0 ) return integer32 is

      lvl,feaNum,sn,flag : integer32;
     -- length : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.feasTest, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      sn := this.sp(depth);
     -- length := this.supType(sn) + 1;
      demics_fTest.class_lvData.get_info
        (this.lv(sn),this.mRepN,this.mFeaIdx,this.mFea);
      lvl := 0;
      flag := iCheck(this,depth,parent,this.lv(sn).fTest(lvl),
                     this.iLv(depth).inif(sn),vrblvl-1);
      lvl := lvl + 1;
      feaNum := 0;
      flag := DEMiCs_Global_Constants.CONTINUE;
      if flag = DEMiCs_Global_Constants.CONTINUE then
        findNode(this,depth,lvl,feaNum,this.lv(sn).fTest,flag,vrblvl-1);
      else
        findNextNode(this,depth,lvl,feaNum,this.lv(sn).fTest,flag,vrblvl-1);
      end if;
      return flag;
    end feasTest;

    procedure upFeasTest
                ( this : in Link_to_mvc; depth : in out integer32;
                  flag : out integer32; vrblvl : in integer32 := 0 ) is

      lvl,feaNum : integer32;
     -- length : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.upFeasTest, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      loop
        demics_iTest.class_iLvData.init
          (this.iLv(depth),this.supN,depth-1,this.iLv(depth - 1).rsp,
           vrblvl-1);
        demics_fTest.class_ftData.delete_addedElem 
          (this.lv(this.sp(depth-1)).node);
        demics_fTest.class_ftData.init_ptr(this.lv(this.sp(depth-1)).node);
        demics_fTest.class_lvData.init_ptr(this.lv(this.sp(depth)));
        demics_fTest.class_lvData.get_info
          (this.lv(this.sp(depth-1)),this.mRepN,this.mFeaIdx,this.mFea);
       -- length := this.supType(this.sp(depth-1)) + 1;
        feaNum := 0;
        lvl := this.supType(this.sp(depth-1));
        findNextNode(this,depth-1,lvl,feaNum,this.lv(this.sp(depth-1)).fTest,
                     flag,vrblvl-1);
        if flag = DEMiCs_Global_Constants.CONTINUE then
          findNode(this,depth-1,lvl,feaNum,this.lv(this.sp(depth-1)).fTest,
                   flag,vrblvl-1);
        end if;
        -- STOP or FNN
        depth := depth - 1;
        exit when ((flag = DEMiCs_Global_Constants.FNN) or (depth = 0));
      end loop; 
    end upFeasTest;

    procedure findMixedCell
                ( this : in Link_to_mvc; depth : in integer32;
                  parent : in demics_fTest.class_theData.Link_to_theData;
                  vrblvl : in integer32 := 0 ) is

      lvl,feaNum,sn,flag : integer32;
     -- length : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.findMixedCell, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      sn := this.sp(depth);
     -- length := this.supType(sn) + 1;
      demics_fTest.class_lvData.get_info
        (this.lv(sn),this.mRepN,this.mFeaIdx,this.mFea);
      lvl := 0;
      flag := iCheck(this,depth,parent,this.lv(sn).fTest(lvl),
                     this.iLv(depth).inif(sn));
      lvl := lvl + 1;
      feaNum := 0;
      flag := DEMiCs_Global_Constants.CONTINUE;
      loop
        if flag = DEMiCs_Global_Constants.CONTINUE then
          findNode(this,depth,lvl,feaNum,this.lv(sn).fTest,
                   flag,vrblvl-1);
          if flag = DEMiCs_Global_Constants.STOP
           then exit;
          end if;
        else
          findNextNode(this,depth,lvl,feaNum,this.lv(sn).fTest,
                       flag,vrblvl-1);
          if flag = DEMiCs_Global_Constants.STOP
           then exit;
          end if;
        end if;
      end loop;
    end findMixedCell;

    procedure findAllMixedCells
                ( this : in Link_to_mvc; depth : in integer32;
                  vrblvl : in integer32 := 0 ) is

      lvl,feaNum,flag : integer32;
     -- length : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.findAllMixedCells, depth : ");
        put(depth,1); put_line(" ...");
      end if;
     -- length := this.supType(depth) + 1;
      demics_fTest.class_lvData.get_info
        (this.lv(depth),this.mRepN,this.mFeaIdx,this.mFea);
      lvl := 0;
      initCheck(this,depth,this.lv(depth).fTest(lvl),vrblvl-1);
      lvl := lvl + 1;
      feaNum := 0;
      flag := DEMiCs_Global_Constants.CONTINUE;
      loop
        if flag = DEMiCs_Global_Constants.CONTINUE then
          findNode(this,depth,lvl,feaNum,this.lv(depth).fTest,flag,vrblvl-1);
          demics_fTest.class_ftData.delete_addedElem(this.lv(depth).node);
          demics_fTest.class_ftData.init_ptr(this.lv(depth).node);
          if flag = DEMiCs_Global_Constants.STOP
           then exit;
          end if;
        else
          findNextNode(this,depth,lvl,feaNum,this.lv(depth).fTest,
                       flag,vrblvl-1);
          demics_fTest.class_ftData.delete_addedElem(this.lv(depth).node);
          demics_fTest.class_ftData.init_ptr(this.lv(depth).node);
          if flag = DEMiCs_Global_Constants.STOP
           then exit;
          end if;
        end if;
      end loop;
    end findAllMixedCells;

    function iCheck
                ( this : Link_to_mvc; depth : integer32;
                  parent : demics_fTest.class_theData.Link_to_theData;
                  data : demics_fTest.class_ftData.Link_to_ftData;
                  curInif : demics_itest.class_inifData.Link_to_inifData;
                  vrblvl : integer32 := 0 ) return integer32 is

      feaNum : integer32 := 0;
      preNbN,flag : integer32;
     -- sn : integer32;
      fst_pivInIdx, sub_fst_pivInIdx : integer32;
      curr : demics_iTest.class_uData.Link_to_uData;

      use demics_iTest.class_uData;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.iCheck, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      if vrblvl > 0 then -- #if DBG_NODE
        put_line("+++++++++++++++++<< iCheck >>+++++++++++++++++");
        info_parent_node(this,depth);
      end if;
      if vrblvl > 0 then -- #if DBG_PRE_INFO
        put_line("<< PRE >>");
        demics_fTest.class_theData.info_p_sol_ptr(parent);
        demics_fTest.class_theData.info_d_sol_ptr(parent);
        demics_fTest.class_theData.info_invB_ptr(parent);
        demics_fTest.class_theData.info_basisIdx_ptr(parent);
        demics_fTest.class_theData.info_nf_pos_ptr(parent);
        demics_fTest.class_theData.info_nbIdx_ptr(parent);
        demics_fTest.class_theData.info_redVec_ptr(parent);
      end if;
      flag := parent.artV;
     -- sn := this.sp(depth);
      case flag is
        when DEMiCs_Global_Constants.TRUE =>
          if vrblvl > 0 then -- #if DBG_NODE
            put_line("<< Art >>");
          end if;
          get_candIdx(this,curInif);
          preNbN := parent.nbN;
          curr := curInif.fHead;
          while curr /= NULL loop
            if vrblvl > 0 then -- #if DBG_NODE
              put("------------------- Idx : ");
              put(curr.supLab+1,1);
              put_line(" -------------------");
            end if;
            iLP_art(this,parent,data,depth,curr.supLab,
                   -- fst_pivInIdx,sub_fst_pivInIdx,
                    preNbN,feaNum,vrblvl-1);
            curr := curr.fNext;
          end loop;
          if vrblvl > 0 then -- #if DBG_NODE
            new_line;
          end if;
          if feaNum <= 1 then
            this.mFea(0) := 0;
            demics_fTest.class_ftData.delete_addedElem(data);
            demics_fTest.class_ftData.init_ptr(data);
            return DEMiCs_Global_Constants.SLIDE;
          else
            return DEMiCs_Global_Constants.CONTINUE;
          end if;
        when DEMiCs_Global_Constants.FALSE =>
          if vrblvl > 0 then -- #if DBG_NODE
            put_line("<< NoArt >>");
          end if;
          demics_simplex.class_simplex.fstRed_candIdx
            (this.the_Simplex,curInif,this.candIdx,fst_pivInIdx,
             sub_fst_pivInIdx,vrblvl-1);
          preNbN := parent.nbN;
          curr := curInif.fHead;
          while curr /= null loop
            if vrblvl > 0 then -- #if DBG_NODE
              put("------------------- Idx : ");
              put(curr.supLab+1,1);
              put_line(" -------------------");
            end if;
            iLP(this,parent,data,depth,curr.supLab,fst_pivInIdx,
                sub_fst_pivInIdx,preNbN,feaNum,vrblvl-1);
            curr := curr.fNext;
          end loop;
          if vrblvl > 0 then -- #if DBG_NODE
            new_line;
          end if;
        when others => null;
      end case;
      if feaNum <= 1 then
        this.mFea(0) := 0;
        demics_fTest.class_ftData.delete_addedElem(data);
        demics_fTest.class_ftData.init_ptr(data);
        return DEMiCs_Global_Constants.SLIDE;
      else
        return DEMiCs_Global_Constants.CONTINUE;
      end if;
    end iCheck;

    procedure iLP ( this : in Link_to_mvc;
                    parent : in demics_fTest.class_theData.Link_to_theData;
                    data : in demics_fTest.class_ftData.Link_to_ftData;
                    depth : in integer32; idx_one : in integer32;
                    fst_pivInIdx : in integer32; 
                    sub_fst_pivInIdx : in integer32; preNbN : in integer32;
                    feaNum : in out integer32; vrblvl : in integer32 := 0 ) is

      flag,fst_pivInIdx2,sub_fst_pivInIdx2,sn,iter : integer32;
     -- termS : integer32;
      reTermS,repIdx,lNfN : integer32;
      fst_redCost : double_float;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.iLP, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      sn := this.sp(depth);
      initMemoryCheck(this,data,depth);
      repIdx := idx_one;
      this.firIdx(sn) := repIdx;
     -- termS := this.termStart(sn);
      reTermS := this.re_termStart(sn);
      lNfN := parent.nfN;
      demics_simplex.class_simplex.get_iNbN_nfN
        (this.the_Simplex,Data.cur,preNbN+this.candIdx(0)-1,lNfN);
      if repIdx < fst_pivInIdx then
        fst_pivInIdx2 := reTermS + fst_pivInIdx - 1;
        sub_fst_pivInIdx2 := preNbN - this.dim + sub_fst_pivInIdx - 1;
      elsif repIdx > fst_pivInIdx then
        fst_pivInIdx2 := reTermS + fst_pivInIdx;
        sub_fst_pivInIdx2 := preNbN - this.dim + sub_fst_pivInIdx;
      end if;
      if repIdx /= fst_pivInIdx then
        demics_fTest.class_ftData.init_info(data);
        demics_ftest.class_ftData.create_rIdx
          (data,preNbN,repIdx,this.candIdx);
        demics_simplex.class_simplex.get_repIdx_candIdx
          (this.the_Simplex,this.candIdx,repIdx);
        demics_simplex.class_simplex.get_parent(this.the_Simplex,parent);
        demics_simplex.class_simplex.get_cur(this.the_Simplex,data.cur);
        fst_redCost := demics_simplex.class_simplex.put_redCost
                         (this.the_Simplex,fst_pivInIdx);
        iter := 0;
        demics_simplex.class_simplex.solLP
          (this.the_Simplex,depth,fst_pivInIdx2,sub_fst_pivInIdx2,
           fst_redCost,DEMiCs_Global_Constants.ICHECK,
           this.termSet(sn),reTermS,preNbN,iter,flag,vrblvl-1);
        this.total_LPs := this.total_LPs + 1.0;
        this.total_1PT := this.total_1PT + 1.0;
        this.lvl_1PT(depth) := this.lvl_1PT(depth) + 1.0;
        if flag = DEMiCs_Global_Constants.OPT then
          if vrblvl > 0 then -- #if DBG_FEA
            put_line("OPT-1");
          end if;
          this.total_iter := this.total_iter + double_float(iter);
          this.total_feasLP := this.total_feasLP + 1.0;
          demics_simplex.class_simplex.get_pivOutNum
            (this.the_Simplex,data.cur);
          demics_fTest.class_theData.joint(data.cur);
          data.cur.fIdx := idx_one;
          if vrblvl > 0 then -- #if DBG_CUR_INFO
            put_line("<< Cur >>");
            demics_fTest.class_ftData.info_cur(data);
          end if;
          this.mFeaIdx(0)(feaNum) := repIdx;
          this.mFea(0) := this.mFea(0) + 1;
          feaNum := feaNum + 1;
          data.cur := data.cur.next;
        elsif flag = DEMiCs_Global_Constants.UNBOUNDED then
          if vrblvl > 0 then -- #if DBG_FEA
            put_line("UNB-1");
          end if;
        else
          put_line("Error: too many iterations at iLP");
          info_parent_node(this,depth);
          put("( "); put(idx_one+1,1); put_line(" )");
          return; -- exit(EXIT_FAILURE);
        end if;
      else
        if vrblvl > 0 then -- #if DBG_FEA
          put_line("OPT-2");
        end if;
        this.mFeaIdx(0)(feaNum) := repIdx;
        this.mFea(0) := this.mFea(0) + 1;
        feaNum := feaNum + 1;
        data.cur.fIdx := idx_one;
        demics_simplex.class_simplex.copy_eye(this.the_Simplex,data.cur); 
        demics_simplex.class_simplex.cal_redVec
          (this.the_Simplex,this.termSet(sn),reTermS,fst_pivInIdx,data.cur,
           vrblvl-1); 
        demics_fTest.class_ftData.iGetPtr(data,parent);
        demics_fTest.class_ftData.get_nbIdx_rIdx
          (data,preNbN,repIdx,this.candIdx,reTermS,parent);
        demics_fTest.class_ftData.init_info(data);
        demics_fTest.class_theData.iJoint(data.cur);
        if vrblvl > 0 then -- #if DBG_CUR_INFO
          put_line("<< Cur_ptr >>");
          demics_ftest.class_ftData.info_cur_ptr(data);
          put_line("<< Cur >>");
          demics_ftest.class_ftData.info_cur_rIdx(data);
        end if;
        data.cur := data.cur.next;
      end if;
    end iLP;

    procedure iLP_art
                ( this : in Link_to_mvc;
                  parent : in demics_fTest.class_theData.Link_to_theData;
                  data : in demics_fTest.class_ftData.Link_to_ftData;
                  depth : in integer32; idx_one : in integer32;
                 -- fst_pivInIdx : in integer32;
                 -- sub_fst_pivInIdx : in integer32;
                  preNbN : in integer32; feaNum : in out integer32;
                  vrblvl : in integer32 := 0 ) is

     -- termS : integer32;
      flag,reTermS,repIdx,sn,iter,lNfN : integer32;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.iLP_art, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      sn := this.sp(depth);
      initMemoryCheck(this,data,depth);
      repIdx := idx_one;
      this.firIdx(sn) := repIdx;
     -- termS := this.termStart(sn);
      reTermS := this.re_termStart(sn);
      lNfN := parent.nfN;
      demics_simplex.class_simplex.get_iNbN_nfN
        (this.the_Simplex,data.cur,preNbN+this.candIdx(0)-1, lNfN);
      demics_simplex.class_simplex.copy_eye(this.the_Simplex,data.cur);
      demics_fTest.class_ftData.copy(data,this.col,parent);
      demics_fTest.class_ftData.iCopy
        (data,preNbN,lNfN,repIdx, --this.termSet(sn),
         reTermS,this.candIdx,parent);
      demics_fTest.class_ftData.init_info(data);
      demics_simplex.class_simplex.get_cur(this.the_Simplex,data.cur);
      iter := 0;
      demics_simplex.class_simplex.solLP_art
        (this.the_Simplex,depth, -- repIdx,fst_pivInIdx,
         preNbN,this.termSet(sn),reTermS,iter,flag,vrblvl-1);
      this.total_LPs := this.total_LPs + 1.0;
      this.total_1PT := this.total_1PT + 1.0;
      this.lvl_1PT(depth) := this.lvl_1PT(depth) + 1.0;
      if flag = DEMiCs_Global_Constants.OPT then
        if vrblvl > 0 then -- #if DBG_FEA
          put_line("OPT-1");
        end if;
        this.total_iter := this.total_iter + double_float(iter);
        this.total_feasLP := this.total_feasLP + 1.0;
        demics_fTest.class_theData.joint(data.cur);
        data.cur.fIdx := idx_one;
        demics_simplex.class_simplex.get_res(this.the_Simplex,data);
        demics_simplex.class_simplex.get_pivOutNum(this.the_Simplex,data.cur);
        this.mFeaIdx(0)(feaNum) := repIdx;
        this.mFea(0) := this.mFea(0) + 1;
        feaNum := feaNum + 1;
        if vrblvl > 0 then -- #if DBG_CUR_INFO
          put_line("Cur");
          demics_fTest.class_ftData.info_cur_ptr(data);
        end if;
        data.cur := data.cur.next;
      elsif flag = DEMiCs_Global_Constants.UNBOUNDED then
        if vrblvl > 0 then -- #if DBG_FEA
          put_line("UNB-1");
        end if;
      else
        put_line("Error: too much iterations at iLP_art");
        info_parent_node(this,depth);
        put("( "); put(idx_one+1,1); put_line(" )");
        return; --  exit(EXIT_FAILURE);
      end if;
    end iLP_art;

    procedure findNode
                ( this : in Link_to_mvc; depth : in integer32;
                  lvl : in out integer32; feaNum : in out integer32;
                  data : in demics_fTest.class_ftData.Link_to_Array_of_ftData;
                  flag : out integer32; vrblvl : in integer32 := 0 ) is

      sn : constant integer32 := this.sp(depth);
      polyDim : constant integer32 := this.supType(sn);
      pre : demics_fTest.class_ftdata.Link_to_ftData := data(lvl-1);
      cur : demics_fTest.class_ftData.Link_to_ftData := data(lvl);

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.findNode, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      if vrblvl > 0 then -- #if DBG_SUC
        info_parent_node(this,depth);
      end if;
      loop
        if vrblvl > 0 then -- #if DBG_SUC
          put("++++++++++ lvl : "); put(lvl,1); put_line(" ++++++++++");
          info_tuple(this,lvl);
        end if;
        mLP(this,pre,cur,data,this.mFeaIdx(lvl-1),this.mFeaIdx(lvl),
            this.mRepN(lvl-1),this.mRepN,this.mFea(lvl-1),depth,
            this.mFea(lvl),lvl,polyDim+1,flag,vrblvl-1);
        if flag = DEMiCs_Global_Constants.NODE then
          feaNum := feaNum + 1;
          flag := DEMiCs_Global_Constants.FNN;
          exit;
        end if;
        if (this.mFea(lvl) <  polyDim - lvl + 1) or
           (lvl = polyDim) then -- SLIDE
          if vrblvl > 0 then -- #if DBG_SUC	
            put_line("-- SLIDE --");
          end if;
          demics_fTest.class_ftData.next_data(pre);
          if lvl /= polyDim
           then demics_fTest.class_ftData.delete_addedElem(cur);
          end if;
          demics_fTest.class_ftData.init_ptr(cur);
          this.mRepN(lvl-1) := this.mRepN(lvl-1) + 1;
          this.mFea(lvl) := 0;
          findUpNode(this,data,pre,cur,lvl,polyDim,depth,vrblvl-1); -- UP
        else -- DOWN
          if vrblvl > 0 then --#if DBG_SUC
            put_line("-- DOWN --");
          end if;
          lvl := lvl + 1;
          pre := Data(lvl - 1);
          cur := Data(lvl);
        end if; 
        if lvl = 0 then
          flag := DEMiCs_Global_Constants.STOP;
          exit;
        end if;
      end loop;
    end findNode;

    procedure findNextNode
                ( this : in Link_to_mvc; depth : in integer32;
                  lvl : in out integer32; feaNum : in out integer32;
                  Data : in demics_fTest.class_ftData.Link_to_Array_of_ftData;
                  flag : out integer32; vrblvl : in integer32 := 0 ) is

      sn : constant integer32 := this.sp(depth);
      polyDim : constant integer32 := this.supType(sn);
      pre,cur : demics_fTest.class_ftData.Link_to_ftData;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.findNextNode, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      if vrblvl > 0 then -- #if DBG_SUC
        info_parent_node(this,depth);
        put("++++++++++ lvl : "); put(lvl,1); put_line(" ++++++++++");
        info_tuple(this,lvl); -- , depth);
      end if;
      mLP(this,data(lvl-1),data(lvl),data,
          this.mFeaIdx(lvl-1),this.mFeaIdx(lvl),
          this.mRepN(polyDim) + this.mRepN(lvl-1),this.mRepN,
          this.mFea(lvl-1),depth,this.mFea(lvl),lvl,polyDim+1,
          flag,vrblvl-1);
      if flag = DEMiCs_Global_Constants.CONTINUE then
        if vrblvl > 0 then -- #if DBG_SUC
          put_line("-- SLIDE or UP --");
        end if;
        demics_fTest.class_ftData.next_data(data(lvl-1));
        this.mRepN(lvl-1) := this.mRepN(lvl-1) + 1;
        this.mRepN(lvl) := 0;
        this.mFea(lvl) := 0;
        findUpNode(this,data,pre,cur,lvl,polyDim,depth,vrblvl-1); -- UP
        if lvl = 0
         then flag := DEMiCs_Global_Constants.STOP;
        end if;
      else
        feaNum := feaNum + 1;
        flag := DEMiCs_Global_Constants.FNN;
      end if;
    end findNextNode;

    procedure findUpNode
                ( this : in Link_to_mvc;
                  data : in demics_fTest.class_ftData.Link_to_Array_of_ftData;
                  pre : out demics_fTest.class_ftData.Link_to_ftData;
                  cur : out demics_fTest.class_ftData.Link_to_ftData;
                  lvl : in out integer32; polyDim : in integer32;
                  depth : in integer32; vrblvl : in integer32 := 0 ) is
    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.findUpNode, depth : ");
        put(depth,1); put_line(" ...");
      end if;
      loop
        if this.mFea(lvl-1)-this.mRepN(lvl-1)-1 < polyDim-lvl+1 then -- UP
          if vrblvl > 0 then -- #if DBG_SUC
            put_line("-- UP --");
          end if;
          this.mFea(lvl) := 0;
          this.mRepN(lvl-1) := 0;
          lvl := lvl - 1;
          demics_fTest.class_ftData.delete_addedElem(data(lvl));
          demics_fTest.class_ftData.init_ptr(data(lvl));
          if lvl = 0 then
            this.mFea(lvl) := 0;
            this.mRepN(lvl) := 0;
            exit;
          else
            this.mFea(lvl) := 0;
            this.mRepN(lvl-1) := this.mRepN(lvl-1) + 1;
            demics_fTest.class_ftData.next_data(data(lvl-1));
            pre := data(lvl-1);
            cur := data(lvl);
          end if;
        else
          exit;
        end if;
      end loop;
    end findUpNode;

    procedure mLP ( this : in Link_to_mvc;
                    pre : in demics_fTest.class_ftData.Link_to_ftData;
                    cur : in demics_fTest.class_ftData.Link_to_ftData;
                  data : in demics_fTest.class_ftData.Link_to_Array_of_ftData;
                    repIdx : in Standard_Integer_Vectors.Link_to_Vector;
                    feaIdx : in Standard_Integer_Vectors.Link_to_Vector;
                    tarIdx : in integer32;
                    mRepN : in Standard_Integer_Vectors.Link_to_Vector;
                    totalN : in integer32; depth : in integer32;
                    feaNum : in out integer32; lvl : in integer32;
                    length : in integer32; flag : out integer32;
                    vrblvl : in integer32 := 0 ) is

      idx2,iter : integer32;
      sn : constant integer32 := this.sp(depth);
      lNbN,lNfN,fst_pivInIdx,fst_sub_pivInIdx : integer32;
      sub_firIdx,sub_tarIdx,colPos,rowStartPos : integer32;
      fst_redCost : double_float;
      target : demics_fTest.class_theData.Link_to_theData;

      use demics_fTest.class_theData;
      use demics_fTest.class_ftData;

    begin
      if vrblvl > 0 then
        put("-> in demics_mvc.class_mvc.mLP, depth : ");
        put(depth,1); put(", lvl : "); put(lvl,1); put_line(" ...");
        if this.lv(sn).node = null
         then put("this.lv("); put(sn,1); put_line(").node = null");
         else put("this.lv("); put(sn,1); put_line(").node /= null");
        end if;
      end if;
      sub_firIdx := mRepN(0);
      sub_tarIdx := mRepN(lvl-1);
      rowStartPos := this.termStart(sn);
      colPos := repIdx(sub_tarIdx) + this.termStart(sn);
      for i in tarIdx+1..totalN-1 loop
        if vrblvl > 0 then -- #if DBG_SUC
          put("-- ( "); put(repIdx(i)+1,1); put_line(" ) --");
        end if;
        if vrblvl > 0 then -- #if DBG_S_PRE_INFO
          put_line("<< Pre_ptr >>");
          demics_fTest.class_ftData.info_parent_ptr(pre);
          put_line("<< Pre >>");
          demics_fTest.class_ftData.info_parent_rIdx(pre);
        end if;
        memoryCheck(this,cur,depth);
        if (table_out(this,colPos,rowStartPos + repIdx(i))
            = DEMiCs_Global_Constants.OPT) then
          if lvl > 0
           then get_firIdx(this,data(0).all,data(1).all,sn,lvl);
          end if;
          target := pre.parent;
         -- flag := checkBasis(this,target,repIdx(i));
          flag := checkBasis(target,repIdx(i));
          demics_simplex.class_simplex.get_mNbN_nfN
            (this.the_Simplex,target,cur.cur);
          demics_fTest.class_theData.put_info
            (target,repIdx(i)-1,idx2,lNbN,lNfN);
          cur.cur.sw := DEMiCs_Global_Constants.OFF;
          if(flag = DEMiCs_Global_Constants.CONTINUE) and (lvl = 1) then
           -- checkAnotherBasis
           --   (this,repIdx(sub_firIdx),i-sub_firIdx,target,flag);
            checkAnotherBasis(repIdx(sub_firIdx),i-sub_firIdx,target,flag);
            if flag = DEMiCs_Global_Constants.OPT then
              demics_simplex.class_simplex.get_mNbN_nfN
                (this.the_Simplex,target,cur.cur);
              demics_fTest.class_theData.put_info
                (target,repIdx(sub_firIdx),idx2,lNbN,lNfN);
              cur.cur.sw := DEMiCs_Global_Constants.ON;
              this.firIdx(sn) := repIdx(i);
            end if;
          end if;
          if flag = DEMiCs_Global_Constants.OPT then
            if vrblvl > 0 then -- #if DBG_FEA
              put_line("OPT-1");
            end if;
            this.total_triLPs_mLP := this.total_triLPs_mLP + 1.0;
            this.actNode(depth) := this.actNode(depth) + 1.0;
            feaIdx(feaNum) := repIdx(i);
            feaNum := feaNum + 1;
            cur.cur.fIdx := repIdx(i);
            demics_fTest.class_ftdata.mGetPtr(cur,target);
            demics_fTest.class_ftData.get_nf_pos(cur,target,lNfN,idx2);
            demics_fTest.class_theData.mJoint(cur.cur);
            demics_fTest.class_ftData.copy_rIdx(cur,target,this.termSet(sn));
            demics_fTest.class_ftData.copy_pivOutIdx(cur,target);      
            if vrblvl > 0 then -- #if DBG_S_CUR_INFO
             -- put_line("<< Cur_ptr >>");
             -- demics_fTest.class_ftData.info_cur_ptr(cur); => crash!
              put_line("<< Cur >>");
              demics_fTest.class_ftData.info_cur_rIdx(cur);
            end if;
            if lvl = length - 1 then
             -- get_tuple_index(this,this.lv(sn).node,data,length);
              get_tuple_index(this.lv(sn).node,data,length,vrblvl-1);
              if depth = this.supN - 1 then
                demics_simplex.class_simplex.calMixedVol
                  (this.the_Simplex,this.lv.all,this.sp,this.supN,vrblvl-1);
              end if;
              mRepN(lvl) := mRepN(lvl) + (i - tarIdx);
              cur.cur := cur.cur.next;
              flag := DEMiCs_Global_Constants.NODE;
              return;
            else
              cur.cur := cur.cur.next;
            end if;
          else -- if flag /= opt
            demics_fTest.class_ftData.copy_rIdx
              (cur,pre.parent,this.termSet(sn));
            demics_fTest.class_ftData.copy_pivOutIdx(cur,pre.parent);
            demics_simplex.class_simplex.get_parent
              (this.the_Simplex,pre.parent);
            demics_simplex.class_simplex.get_cur
              (this.the_Simplex,cur.cur);
            fst_sub_pivInIdx := -1 * idx2 - 1;
            fst_pivInIdx := pre.parent.nbIdx_ptr(fst_sub_pivInIdx);
            fst_redCost := pre.parent.redVec_ptr(fst_pivInIdx);
            iter := 0;
            demics_simplex.class_simplex.solLP
              (this.the_Simplex,depth,fst_pivInIdx,fst_sub_pivInIdx,
               fst_redCost,DEMiCs_Global_Constants.MCHECK,
               this.termSet(sn),this.re_termStart(sn),lNbN,
               iter,flag,vrblvl-1);
            this.total_LPs := this.total_LPs + 1.0;
            this.total_2PT := this.total_2PT + 1.0;
            this.lvl_2PT(depth) := this.lvl_2PT(depth) + 1.0;
            if(flag = DEMiCs_Global_Constants.OPT) then
              if vrblvl > 0 then -- #if DBG_FEA	  
                put_line("OPT-2");
              end if;
              this.total_iter := this.total_iter + double_float(iter);
              this.total_feasLP := this.total_feasLP + 1.0;
              this.actNode(depth) := this.actNode(depth) + 1.0;
              demics_simplex.class_simplex.get_pivOutNum
                (this.the_Simplex,cur.cur);
              demics_fTest.class_theData.joint(cur.cur);	
              demics_fTest.class_ftData.decrease_nfN(cur);
              cur.cur.fIdx := repIdx(i);
              feaIdx(feaNum) := repIdx(i);
              feaNum := feaNum + 1;
              if vrblvl > 0 then -- #if DBG_S_CUR_INFO 
                put_line("<< Cur >>");
                demics_fTest.class_ftData.info_cur(cur);
              end if;
              if lvl = length - 1 then
               -- get_tuple_index
               --   (this,this.lv(sn).Node,this.lv(sn).fTest,length);
                get_tuple_index
                  (this.lv(sn).Node,this.lv(sn).fTest,length,vrblvl-1);
                if depth = this.supN - 1 then
                  demics_simplex.class_simplex.calMixedVol
                    (this.the_Simplex,this.lv.all,this.sp,this.supN,vrblvl-1);
                end if;
                this.mRepN(lvl) := this.mRepN(lvl) + (i - tarIdx);
                cur.cur := cur.cur.next;
                flag := DEMiCs_Global_Constants.NODE;
                return;
              else
                cur.cur := cur.cur.next;
              end if;
            elsif flag = DEMiCs_Global_Constants.UNBOUNDED then
              if vrblvl > 0 then -- #if DBG_FEA
                put_line("UNB-1");
              end if;
              demics_fTest.class_ftData.init_info(cur);
            else
              put_line("Error: too much iterations at solLP");
              info_parent_node(this,depth);
              info_tuple(this,lvl); -- , depth);
              put("( "); put(repIdx(i)+1,1); put_line(" )");
              flag := -1; -- EXIT_FAILURE;
              return;
            end if;
          end if; -- ends else case if /= opt
        else -- if table_out() /= opt
          if vrblvl > 0 then -- #if DBG_FEA
            put_line("UNB-table");
          end if;
          demics_fTest.class_ftData.init_info(cur);
        end if; 
      end loop;
      flag := DEMiCs_Global_Constants.CONTINUE;
    end mLP;

    function checkBasis
               -- ( this : Link_to_mvc;
                ( target : demics_fTest.class_theData.Link_to_theData;
                  sub_sIdx : integer32 ) return integer32 is
    begin
      if target.rIdx(sub_sIdx - 1) >= 0
       then return DEMiCs_Global_Constants.OPT;
       else return DEMiCs_Global_Constants.CONTINUE;
      end if;
    end checkBasis;

    procedure checkAnotherBasis
               -- ( this : in Link_to_mvc;
                ( repIdx : in integer32; dist : in integer32;
                  target : in out demics_fTest.class_theData.Link_to_theData;
                  result : out integer32 ) is
    begin
      for i in 0..dist-1 loop
        target := target.next;
      end loop;
      if target.rIdx(repIdx) >= 0
       then result := DEMiCs_Global_Constants.OPT;
       else result := DEMiCs_Global_Constants.CONTINUE;
      end if;
    end checkAnotherBasis;

    procedure get_firIdx
                ( this : in Link_to_mvc;
                  data_a : in demics_fTest.class_ftData.ftData;
                  data_b : in demics_fTest.class_ftData.ftData;
                  sn : in integer32; lvl : in integer32 ) is

      off : constant := DEMiCs_Global_Constants.OFF;

    begin
      if lvl = 1 then
        this.firIdx(sn) := data_a.parent.fIdx;
      else
        if data_b.parent.sw = off
         then this.firIdx(sn) := data_a.parent.fIdx;
         else this.firIdx(sn) := data_b.parent.fIdx;
        end if;
      end if;
    end get_firIdx;

    procedure info_cpuTime
                ( this : in Link_to_mvc;
                  cpuTime_start : in double_float;
                  cpuTime_end : in double_float ) is
    begin
      null;
    end info_cpuTime;

    procedure info_final ( this : in Link_to_mvc ) is
    begin
      null;
    end info_final;

    procedure enum ( this : in Link_to_mvc; vrblvl : in integer32 := 0 ) is

      depth : integer32 := 0;
      flag : integer32;
      ftlast : integer32;

      use demics_fTest.class_theData;
      use demics_fTest.class_ftData;

    begin
      if vrblvl > 0 then
        put_line("-> in demics_mvc.class_mvc.enum ...");
      end if;
      demics_reltab.class_reltab.makeTable
        (this.the_Reltab,this.total_unbLP_tab,vrblvl-1);
      this.table := this.the_Reltab.table;
      for i in this.lv'range loop
        put("checking this.lv("); put(i,1); put(").node");
        if this.lv(i).node /= null then
          put_line(" /= null");
        else
          put_line(" = null, assigning to last component ...");
          ftlast := this.lv(i).fTest'last;
          this.lv(i).node := this.lv(i).fTest(ftlast);
        end if;
        put_line("checking fTest ...");
        if this.lv(i).fTest = null then
          put("this.lv("); put(i,1); put_line(").fTest = null");
        else
          put("this.lv("); put(i,1); put_line(").fTest /= null");
          for j in this.lv(i).fTest'range loop
            put("this.lv("); put(i,1); 
            put(").fTest("); put(j,1); put(")");
            if this.lv(i).fTest(j) = null
             then put_line(" = null");
             else put_line(" /= null");
            end if;
          end loop;
        end if;
      end loop;
      if this.supN = 1 then
        findAllMixedCells(this,depth,vrblvl-1);
      else
        initFeasTest(this,depth,vrblvl-1);
        depth := depth + 1;
        loop
          if depth = this.supN - 1 then
            flag := chooseSup
                      (this,depth-1,this.lv(this.sp(depth-1)).Node.parent,
                       this.iLv(depth-1).inif,this.iLv(depth).inif,vrblvl-1);
            if flag = DEMiCs_Global_Constants.CONTINUE then
              findMixedCell
                (this,depth,this.lv(this.sp(depth-1)).Node.parent,vrblvl-1);
              flag := DEMiCs_Global_Constants.STOP;
            end if;
          else
            flag := chooseSup
                      (this,depth-1,this.lv(this.sp(depth-1)).Node.parent,
                       this.iLv(depth-1).inif,this.iLv(depth).inif,vrblvl-1);
            if flag = DEMiCs_Global_Constants.CONTINUE then
              flag := feasTest(this,depth,
                               this.lv(this.sp(depth-1)).Node.parent,
                               vrblvl-1);
            end if;
          end if;
          if flag = DEMiCs_Global_Constants.STOP then -- SLIDE or UP
            upFeasTest(this,depth,flag,vrblvl-1);
            if flag = DEMiCs_Global_Constants.STOP
             then exit;
            end if;
          end if;
          depth := depth + 1; -- DOWN
        end loop;
      end if;
      demics_simplex.class_simplex.info_mv(this.the_Simplex);
    end enum;

  end class_mvc;

end demics_mvc;
