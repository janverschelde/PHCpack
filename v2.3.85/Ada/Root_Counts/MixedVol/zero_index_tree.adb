with unchecked_deallocation;
with text_io;

package body Zero_Index_Tree is

  procedure L0IdxNode_Init 
              ( p : in out Link_to_L0IdxNode; i : in integer32 ) is
  begin
    p.idx := i;
    p.R := null;
    p.D := null;
  end L0IdxNode_Init;

  procedure L0_IML_Init ( p : in out Link_to_L0_IML ) is
  begin
    p.L0head := new L0IdxNode;
    L0IdxNode_Init(p.L0head,-1);
    p.L0curr := p.L0head;
    p.L0prev := p.L0head;
    p.LP1 := new LPPL;
    LPPL_Init(p.LP1,null,null);
  end L0_IML_Init;

  procedure L0_Migrate
              ( p : in out Link_to_L0_IML;
                inp : in out Link_to_IndexNode;
                status : out integer32 ) is
  begin
    if inp = null
     then text_io.put_line("Migrate: inp = null!");
          status := 1; return;
    end if;
    if p.L0head.D /= null then
      p.L0prev := p.L0head.D;
      inp.idx := p.L0prev.idx;
      inp.S := p.L0prev.R;
      p.L0head.D := p.L0prev.D;
      Clear(p.L0prev);
      status := 1;
    else
      Clear(p.L0head);
      status := 0;
    end if;
  end L0_Migrate;

  procedure L0_FindInR
               ( p : in out Link_to_L0_IML;
                 IDX : in integer32; found : out boolean ) is
  begin
    p.curr := p.prev.S;
    while p.curr /= null loop
      if IDX <= p.curr.idx then
        if IDX = p.curr.idx
         then found := true; return;
         else found := false; return; -- prev.idx < IDX < curr.idx
        end if;
      end if;
      p.prev := p.prev.S;
      p.curr := p.curr.S;
    end loop;
    found := false; -- all indices are smaller than IDX,
                    -- i.e.: < prev.idx < IDX (prev.S = curr = 0)
  end L0_FindInR;

  procedure L0_FindInD
               ( p : in out Link_to_L0_IML;
                 IDX : in integer32; found : out boolean ) is
  begin
    p.L0curr := p.L0prev.D;
    while p.L0curr /= null loop
      if IDX <= p.L0curr.idx then
        if IDX = p.L0curr.idx
         then found := true; return;
         else found := false; return;  -- L0prev.idx < IDX < L0curr.idx
        end if;
      end if;
      p.L0prev := p.L0prev.D;
      p.L0curr := p.L0curr.D;
    end loop;
    found := false; -- all indices are smaller than IDX,
                    -- i.e.: < L0prev.idx < IDX (L0prev.D = L0curr = null)
  end L0_FindInD;

  procedure L0_Add1
               ( p : in out Link_to_L0_IML; n : in integer32; 
                 J : in Standard_Integer_Vectors.Link_to_Vector; 
                 d : in integer32; 
                 I : in Standard_Integer_Vectors.Link_to_Vector; 
                 X : in Standard_Floating_Vectors.Link_to_Vector;
                 A : in Standard_Floating_Matrices.Link_to_Matrix ) is

    LPnotused : boolean := true;
    lpp : Link_to_LPdata := new LPdata;
    k : integer32;
    found : boolean;

  begin
   -- text_io.put_line("executing L0_Add1");
    LPdata_Init(lpp,d,I,X,A);
    p.L0prev := p.L0head;
    for i in 0..(n-1) loop                -- run through all points in J
      k := i+2;
      L0_FindInD(p,J(i),found);
      if found then                               -- then L0curr /= null
        if i+1 < n then
          if p.L0curr.R /= null then
            p.L0prev := p.L0curr;
            p.curr := p.L0curr.R;
            if p.curr.idx > J(i+1) then
              p.curr := new IndexNode;
              IndexNode_Init(p.curr,J(i+1),lpp);
              p.curr.S := p.L0prev.R;
              p.L0prev.R := p.curr;
              LPnotused := false;
            elsif p.curr.idx < J(i+1) then
              k := i+1;
            end if;
            while k < n loop
              p.prev := p.curr;
              L0_FindInR(p,J(k),found);
              if not found then
                p.curr := new IndexNode;
                IndexNode_Init(p.curr,J(k),lpp);
                p.curr.S := p.prev.S;
                p.prev.S := p.curr;
                LPnotused := false;
              end if;
              k := k + 1;
            end loop;
          else
            p.curr := new IndexNode;
            IndexNode_Init(p.curr,J(k+1),lpp);
            p.L0curr.R := p.curr;
            LPnotused := false;
            while k < n loop
              p.prev := p.curr;                       -- Add2S(J[k],lpp)
              p.curr := new IndexNode;
              IndexNode_Init(p.curr,J(k),lpp);
              p.curr.S := p.prev.S;
              p.prev.S := p.curr;
              k := k + 1;
            end loop;
          end if;
        end if;
      else                                      -- add J[i] after L0prev
        p.L0curr := new L0IdxNode;
        L0IdxNode_Init(p.L0curr,J(i));
        LPnotused := false;
        p.L0curr.D := p.L0prev.D;
        p.L0prev.D := p.L0curr;
        if i+1 < n then
          p.curr := new IndexNode;
          IndexNode_Init(p.curr,J(i+1),lpp);
          p.L0curr.R := p.curr;
          while k < n loop
            p.prev := p.curr;                       -- Add2S(J[k],lpp)
            p.curr := new IndexNode;
            IndexNode_Init(p.curr,J(k),lpp);
            p.curr.S := p.prev.S;
            p.prev.S := p.curr;
            k := k + 1;
          end loop;
        end if;
      end if;
      p.L0prev := p.L0curr;
    end loop;
    if LPnotused then
      Standard_Integer_Vectors.Clear(lpp.JJJ);
      Standard_Floating_Vectors.Clear(lpp.xxx);
      Standard_Floating_Matrices.Clear(lpp.INV);
      Clear(lpp);
    else
      p.LP1.next := new LPPL;
      LPPL_Init(p.LP1.next,lpp,p.LP1.next);
    end if;
  end L0_Add1;

  procedure L0_Add2 
               ( p : in out Link_to_L0_IML;
                 J : in Standard_Integer_Vectors.Link_to_Vector;
                 d : in integer32;
                 I : in Standard_Integer_Vectors.Link_to_Vector;
                 X : in Standard_Floating_Vectors.Link_to_Vector;
                 A : in Standard_Floating_Matrices.Link_to_Matrix ) is

    lpp : Link_to_LPdata;
    found : boolean;

  begin
   -- text_io.put_line("executing L0_Add2");
    p.L0prev := p.L0head;
    for k in 0..integer32(1) loop            -- run through all points in J
      L0_FindInD(p,J(k),found);
      if found then                                  -- then L0curr /= null
        if k = 0 then
          if p.L0curr.R /= null then
            p.L0prev := p.L0curr;
            p.curr := p.L0curr.R;
            if p.curr.idx > J(1) then
              lpp := new LPdata;
              LPdata_Init(lpp,d,I,X,A);
              p.LP1.next := new LPPL;
              LPPL_Init(p.LP1.next,lpp,p.LP1.next);
              p.curr := new IndexNOde;
              IndexNode_Init(p.curr,J(1),lpp);
              p.curr.S := p.L0prev.R;
              p.L0prev.R := p.curr;
            elsif p.curr.idx < J(1) then
              p.prev := p.curr;
              L0_FindInR(p,J(1),found);
              if not found then                        -- Add2S(J[1],lpp)
                lpp := new LPdata;
                LPdata_Init(lpp,d,I,X,A);
                p.curr := new IndexNode;
                IndexNode_Init(p.curr,J(1),lpp);
                p.curr.S := p.prev.S;
                p.prev.S := p.curr;
              else
                return;             -- J[1] is in the right branch of J[0]
              end if;
            end if;
          else
            lpp := new LPdata;
            LPdata_Init(lpp,d,I,X,A);
            p.LP1.next := new LPPL;
            LPPL_Init(p.LP1.next,lpp,p.LP1.next);
            p.L0curr.R := new IndexNode;
            IndexNode_Init(p.L0curr.R,J(1),lpp);       -- Add2E(J[1],lpp)
          end if;
        end if;
      else
        if k = 0 then
          lpp := new LPdata;
          LPdata_Init(lpp,d,I,X,A);
          p.LP1.next := new LPPL;
          LPPL_Init(p.LP1.next,lpp,p.LP1.next);
          p.L0curr := new L0IdxNode;
          L0IdxNode_Init(p.L0curr,J(k));
          p.L0curr.D := p.L0prev.D;
          p.L0prev.D := p.L0curr;
          p.L0curr.R := new IndexNode;
          IndexNode_Init(p.L0curr.R,J(1),lpp);         -- Add2E(J[1],lpp)
        else
          p.L0curr := new L0IdxNode;
          L0IdxNode_Init(p.L0curr,J(k));
          p.L0curr.D := p.L0prev.D;
          p.L0prev.D := p.L0curr;
        end if;
      end if;
      p.L0prev := p.L0curr;
    end loop;
  end L0_Add2;

  procedure L0_IML_Del ( p : in out Link_to_L0_IML ) is
  begin
    L0_Free(p);
    Clear(p.LP1);   -- L0head is removed by Migrate when empty
  end L0_IML_Del;

  procedure L0_Free ( li : in out Link_to_L0_IML ) is

    P : Link_to_LPPL := li.LP1.next;

  begin
    while P /= null loop
      li.LP1.next := P.next;
      Standard_Integer_Vectors.Clear(P.addr.JJJ);
      Standard_Floating_Vectors.Clear(P.addr.xxx);
      Standard_Floating_Matrices.clear(P.addr.INV);
      Clear(P);
      P := li.LP1.next;
    end loop;  -- Multi index links and L0head would be removed by Migrate() !
  end L0_Free;

  procedure Clear ( p : in out Link_to_L0IdxNode ) is

    procedure free is new unchecked_deallocation(L0IdxNode,Link_to_L0IdxNode);

  begin
    free(p);
  end Clear;

  procedure Clear ( p : in out Link_to_L0_IML ) is

    procedure free is new unchecked_deallocation(L0_IML,Link_to_L0_IML);

  begin
    free(p);
  end Clear;

end Zero_Index_Tree;
