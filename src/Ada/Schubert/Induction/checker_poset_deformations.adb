with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_Matrices;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Checker_Moves;
with Checker_Posets_io;
with Checker_Localization_Patterns;
with Checker_Homotopies;
with Setup_Flag_Homotopies;             use Setup_Flag_Homotopies;
with Moving_Flag_Continuation;          use Moving_Flag_Continuation;

package body Checker_Poset_Deformations is

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in out Standard_Complex_Matrices.Matrix;
                ls : in out Standard_Complex_Solutions.Link_to_Solution; 
                tol : in double_float; unhappy : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    start_mf : Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    dim,ptr,homtp,ctr,ind,fc : integer32;
    stay_child : boolean;
    fail : boolean := false;

    use Standard_Complex_Matrices;

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_Path_in_Poset 1 ...");
    end if;
    new_line(file);
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      put(file,"No tracking for path "); put(file,count,1);
      put(file," because "); Checker_Posets_io.Write(file,leaf.cols,cnd);
      put_line(file," is not happy.");
      unhappy := true;
    else
      unhappy := false;
      put(file,"Tracking path "); put(file,count,1);
      put(file," in poset starting at a happy ");
      Checker_Posets_io.Write(file,leaf.cols,cnd); put_line(file," ...");
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        Checker_Posets_io.Write_Node_in_Path(file,n,k,ps,path,i);
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        fc := Checker_Moves.Falling_Checker(p);
        put(file,"Calling Transformation with index = ");
        put(file,integer32(q(fc)),1);
        put(file," where fc = "); put(file,fc,1);
        put(file," and q = "); put(file,q); new_line(file);
        t := Checker_Localization_Patterns.Transformation(n,integer32(q(fc)));
       -- put_line(file,"The pattern t for the numerical transformation ");
       -- put(file,t);
       -- put_line(file,"The numerical transformation :");
       -- Write_Standard_Moving_Flag(file,Numeric_Transformation(t));
       -- put_line(file,"The moving flag before the update :");
       -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,mf);
        start_mf := mf;
        mf := mf*Numeric_Transformation(t);
       -- put_line(file,"The moving flag after the update :");
       -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,mf);
        Checker_Homotopies.Define_Generalizing_Homotopy
           (file,n,q,qr,qc,stay_child,homtp,ctr);
        Initialize_Symbol_Table(n,k,q,qr,qc,dim);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay
            (file,n,k,ctr,ind,q,p,qr,qc,pr,pc,verify,minrep,cond,mf,vf,ls,
             fail,vrblvl-1);
        elsif homtp = 1 then
          Stay_Homotopy(file,n,k,ctr,ind,q,p,qr,qc,pr,pc,verify,minrep,tosqr,
                        cond,vf,mf,start_mf,ls,tol,fail,rpt,vrblvl-1);
        else -- homtp = 2
          Setup_Flag_Homotopies.Add_t_Symbol;
          Swap_Homotopy(file,n,k,ctr,ind,q,p,qr,qc,pr,pc,verify,minrep,tosqr,
                        cond,mf,start_mf,vf,ls,tol,fail,rpt,vrblvl-1);
        end if;
        if fail then
          put_line(file,"no longer a valid solution, abort tracking");
          new_line(file);
          exit;
        end if;
      end loop;
    end if;
  end Track_Path_in_Poset;

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf : in out DoblDobl_Complex_Matrices.Matrix;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution; 
                tol : in double_float; unhappy : out boolean;
                vrblvl : in integer32 := 0 ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    start_mf : DoblDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    dim,ptr,homtp,ctr,ind,fc : integer32;
    stay_child : boolean;
    fail : boolean := false;

    use DoblDobl_Complex_Matrices;

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_Path_in_Poset 2 ...");
    end if;
    new_line(file);
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      put(file,"No tracking for path "); put(file,count,1);
      put(file," because "); Checker_Posets_io.Write(file,leaf.cols,cnd);
      put_line(file," is not happy.");
      unhappy := true;
    else
      unhappy := false;
      put(file,"Tracking path "); put(file,count,1);
      put(file," in poset starting at a happy ");
      Checker_Posets_io.Write(file,leaf.cols,cnd); put_line(file," ...");
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        Checker_Posets_io.Write_Node_in_Path(file,n,k,ps,path,i);
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        fc := Checker_Moves.Falling_Checker(p);
        put(file,"Calling Transformation with index = ");
        put(file,integer32(q(fc)),1);
        put(file," where fc = "); put(file,fc,1);
        put(file," and q = "); put(file,q); new_line(file);
        t := Checker_Localization_Patterns.Transformation(n,integer32(q(fc)));
       -- put_line(file,"The pattern t for the numerical transformation ");
       -- put(file,t);
       -- put_line(file,"The numerical transformation :");
       -- Write_Standard_Moving_Flag(file,Numeric_Transformation(t));
       -- put_line(file,"The moving flag before the update :");
       -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,mf);
        start_mf := mf;
        mf := mf*Numeric_Transformation(t);
       -- put_line(file,"The moving flag after the update :");
       -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,mf);
        Checker_Homotopies.Define_Generalizing_Homotopy
           (file,n,q,qr,qc,stay_child,homtp,ctr);
        Initialize_Symbol_Table(n,k,q,qr,qc,dim);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay
            (file,n,k,ctr,ind,q,p,qr,qc,pr,pc,verify,minrep,cond,mf,vf,ls,
             fail,vrblvl-1);
        elsif homtp = 1 then
          Stay_Homotopy(file,n,k,ctr,ind,q,p,qr,qc,pr,pc,verify,minrep,tosqr,
                        cond,vf,mf,start_mf,ls,tol,fail,vrblvl-1);
        else -- homtp = 2
          Setup_Flag_Homotopies.Add_t_Symbol;
          Swap_Homotopy(file,n,k,ctr,ind,q,p,qr,qc,pr,pc,verify,minrep,tosqr,
                        cond,mf,start_mf,vf,ls,tol,fail,vrblvl-1);
        end if;
        if fail then
          put_line(file,"no longer a valid solution, abort tracking");
          new_line(file);
          exit;
        end if;
      end loop;
    end if;
  end Track_Path_in_Poset;

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf : in out QuadDobl_Complex_Matrices.Matrix;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution; 
                tol : in double_float; unhappy : out boolean;
                vrblvl : in integer32 := 0 ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    start_mf : QuadDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    dim,ptr,homtp,ctr,ind,fc : integer32;
    stay_child : boolean;
    fail : boolean := false;

    use QuadDobl_Complex_Matrices;

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_Path_in_Poset 3 ...");
    end if;
    new_line(file);
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      put(file,"No tracking for path "); put(file,count,1);
      put(file," because "); Checker_Posets_io.Write(file,leaf.cols,cnd);
      put_line(file," is not happy.");
      unhappy := true;
    else
      unhappy := false;
      put(file,"Tracking path "); put(file,count,1);
      put(file," in poset starting at a happy ");
      Checker_Posets_io.Write(file,leaf.cols,cnd); put_line(file," ...");
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        Checker_Posets_io.Write_Node_in_Path(file,n,k,ps,path,i);
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        fc := Checker_Moves.Falling_Checker(p);
        put(file,"Calling Transformation with index = ");
        put(file,integer32(q(fc)),1);
        put(file," where fc = "); put(file,fc,1);
        put(file," and q = "); put(file,q); new_line(file);
        t := Checker_Localization_Patterns.Transformation(n,integer32(q(fc)));
       -- put_line(file,"The pattern t for the numerical transformation ");
       -- put(file,t);
       -- put_line(file,"The numerical transformation :");
       -- Write_Standard_Moving_Flag(file,Numeric_Transformation(t));
       -- put_line(file,"The moving flag before the update :");
       -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,mf);
        start_mf := mf;
        mf := mf*Numeric_Transformation(t);
       -- put_line(file,"The moving flag after the update :");
       -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,mf);
        Checker_Homotopies.Define_Generalizing_Homotopy
          (file,n,q,qr,qc,stay_child,homtp,ctr);
        Initialize_Symbol_Table(n,k,q,qr,qc,dim);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay
            (file,n,k,ctr,ind,q,p,qr,qc,pr,pc,verify,minrep,cond,mf,vf,ls,
             fail,vrblvl-1);
        elsif homtp = 1 then
          Stay_Homotopy(file,n,k,ctr,ind,q,p,qr,qc,pr,pc,verify,minrep,tosqr,
                        cond,vf,mf,start_mf,ls,tol,fail,vrblvl-1);
        else -- homtp = 2
          Setup_Flag_Homotopies.Add_t_Symbol;
          Swap_Homotopy(file,n,k,ctr,ind,q,p,qr,qc,pr,pc,verify,minrep,tosqr,
                        cond,mf,start_mf,vf,ls,tol,fail,vrblvl-1);
        end if;
        if fail then
          put_line(file,"no longer a valid solution, abort tracking");
          new_line(file);
          exit;
        end if;
      end loop;
    end if;
  end Track_Path_in_Poset;

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k,nt : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in out Standard_Complex_Matrices.Matrix;
                start : in Standard_Complex_Solutions.Solution_List;
                sols : out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; unhappy : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    start_mf : Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    dim,ptr,homtp,ctr,ind,fc : integer32;
    stay_child : boolean;
    fail : boolean := false;

    use Standard_Complex_Matrices;
    use Standard_Complex_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_Path_in_Poset 4 ...");
    end if;
    new_line(file);
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      put(file,"No tracking for path "); put(file,count,1);
      put(file," because "); Checker_Posets_io.Write(file,leaf.cols,cnd);
      put_line(file," is not happy.");
      unhappy := true;
    else
      if Is_Null(start) then
        put_line(file,"No start solutions ?  Abort track path in poset.");
        return;
      else
        put(file,"In track path in poset, number of start solutions : ");
        put(file,Length_Of(start),1); put_line(file,".");
        put_line(file,"THE START SOLUTIONS :");
        put(file,Length_Of(start),natural32(Head_Of(start).n),start);
      end if;
      unhappy := false;
      put(file,"Tracking path "); put(file,count,1);
      put(file," in poset starting at a happy ");
      Checker_Posets_io.Write(file,leaf.cols,cnd); put_line(file," ...");
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      Copy(start,sols);
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        Checker_Posets_io.Write_Node_in_Path(file,n,k,ps,path,i);
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        fc := Checker_Moves.Falling_Checker(p);
        put(file,"Calling Transformation with index = ");
        put(file,integer32(q(fc)),1);
        put(file," where fc = "); put(file,fc,1);
        put(file," and q = "); put(file,q); new_line(file);
        t := Checker_Localization_Patterns.Transformation(n,integer32(q(fc)));
       -- put_line(file,"The pattern t for the numerical transformation ");
       -- put(file,t);
       -- put_line(file,"The numerical transformation :");
       -- Write_Standard_Moving_Flag(file,Numeric_Transformation(t));
       -- put_line(file,"The moving flag before the update :");
       -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,mf);
        start_mf := mf;
        mf := mf*Numeric_Transformation(t);
       -- put_line(file,"The moving flag after the update :");
       -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,mf);
        Checker_Homotopies.Define_Generalizing_Homotopy
          (file,n,q,qr,qc,stay_child,homtp,ctr);
        Initialize_Symbol_Table(n,k,q,qr,qc,dim);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay
            (file,n,k,ctr,q,p,qr,qc,pr,pc,verify,minrep,cond,
             mf,vf,sols,tol,fail,vrblvl-1);
        elsif homtp = 1 then
          Stay_Homotopy(file,n,k,ctr,nt,q,p,qr,qc,pr,pc,verify,minrep,
                        tosqr,cond,vf,mf,start_mf,sols,tol,fail,rpt,vrblvl-1);
        else -- homtp = 2
          Setup_Flag_Homotopies.Add_t_Symbol;
          Swap_Homotopy(file,n,k,ctr,nt,q,p,qr,qc,pr,pc,verify,minrep,
                        tosqr,cond,mf,start_mf,vf,sols,tol,fail,rpt,vrblvl-1);
        end if;
        if fail then
          put_line(file,"no longer a valid solution, abort tracking");
          new_line(file);
          unhappy := true; -- prevent from being concatenated
          exit;
        end if;
      end loop;
    end if;
  end Track_Path_in_Poset;

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k,nt : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf : in out DoblDobl_Complex_Matrices.Matrix;
                start : in DoblDobl_Complex_Solutions.Solution_List;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; unhappy : out boolean;
                vrblvl : in integer32 := 0 ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    start_mf : DoblDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    dim,ptr,homtp,ctr,ind,fc : integer32;
    stay_child : boolean;
    fail : boolean := false;

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_Path_in_Poset 5 ...");
    end if;
    new_line(file);
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      put(file,"No tracking for path "); put(file,count,1);
      put(file," because "); Checker_Posets_io.Write(file,leaf.cols,cnd);
      put_line(file," is not happy.");
      unhappy := true;
    else
      if Is_Null(start) then
        put_line(file,"No start solutions ?  Abort track path in poset.");
        return;
      else
        put(file,"In track path in poset, number of start solutions : ");
        put(file,Length_Of(start),1); put_line(file,".");
        put_line(file,"THE START SOLUTIONS :");
        put(file,Length_Of(start),natural32(Head_Of(start).n),start);
      end if;
      unhappy := false;
      put(file,"Tracking path "); put(file,count,1);
      put(file," in poset starting at a happy ");
      Checker_Posets_io.Write(file,leaf.cols,cnd); put_line(file," ...");
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      Copy(start,sols);
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        Checker_Posets_io.Write_Node_in_Path(file,n,k,ps,path,i);
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        fc := Checker_Moves.Falling_Checker(p);
        put(file,"Calling Transformation with index = ");
        put(file,integer32(q(fc)),1);
        put(file," where fc = "); put(file,fc,1);
        put(file," and q = "); put(file,q); new_line(file);
        t := Checker_Localization_Patterns.Transformation(n,integer32(q(fc)));
       -- put_line(file,"The pattern t for the numerical transformation ");
       -- put(file,t);
       -- put_line(file,"The numerical transformation :");
       -- Write_Standard_Moving_Flag(file,Numeric_Transformation(t));
       -- put_line(file,"The moving flag before the update :");
       -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,mf);
        start_mf := mf;
        mf := mf*Numeric_Transformation(t);
       -- put_line(file,"The moving flag after the update :");
       -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,mf);
        Checker_Homotopies.Define_Generalizing_Homotopy
          (file,n,q,qr,qc,stay_child,homtp,ctr);
        Initialize_Symbol_Table(n,k,q,qr,qc,dim);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay
            (file,n,k,ctr,q,p,qr,qc,pr,pc,verify,minrep,cond,
             mf,vf,sols,tol,fail,vrblvl-1);
        elsif homtp = 1 then
          Stay_Homotopy(file,n,k,ctr,nt,q,p,qr,qc,pr,pc,verify,minrep,
                        tosqr,cond,vf,mf,start_mf,sols,tol,fail,vrblvl-1);
        else -- homtp = 2
          Setup_Flag_Homotopies.Add_t_Symbol;
          Swap_Homotopy(file,n,k,ctr,nt,q,p,qr,qc,pr,pc,verify,minrep,
                        tosqr,cond,mf,start_mf,vf,sols,tol,fail,vrblvl-1);
        end if;
        if fail then
          put_line(file,"no longer a valid solution, abort tracking");
          new_line(file);
          unhappy := true; -- prevent from being concatenated
          exit;
        end if;
      end loop;
    end if;
  end Track_Path_in_Poset;

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k,nt : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf : in out QuadDobl_Complex_Matrices.Matrix;
                start : in QuadDobl_Complex_Solutions.Solution_List;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; unhappy : out boolean;
                vrblvl : in integer32 := 0 ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    start_mf : QuadDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    dim,ptr,homtp,ctr,ind,fc : integer32;
    stay_child : boolean;
    fail : boolean := false;

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_Path_in_Poset 6 ...");
    end if;
    new_line(file);
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      put(file,"No tracking for path "); put(file,count,1);
      put(file," because "); Checker_Posets_io.Write(file,leaf.cols,cnd);
      put_line(file," is not happy.");
      unhappy := true;
    else
      if Is_Null(start) then
        put_line(file,"No start solutions ?  Abort track path in poset.");
        return;
      else
        put(file,"In track path in poset, number of start solutions : ");
        put(file,Length_Of(start),1); put_line(file,".");
        put_line(file,"THE START SOLUTIONS :");
        put(file,Length_Of(start),natural32(Head_Of(start).n),start);
      end if;
      unhappy := false;
      put(file,"Tracking path "); put(file,count,1);
      put(file," in poset starting at a happy ");
      Checker_Posets_io.Write(file,leaf.cols,cnd); put_line(file," ...");
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      Copy(start,sols);
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        Checker_Posets_io.Write_Node_in_Path(file,n,k,ps,path,i);
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        fc := Checker_Moves.Falling_Checker(p);
        put(file,"Calling Transformation with index = ");
        put(file,integer32(q(fc)),1);
        put(file," where fc = "); put(file,fc,1);
        put(file," and q = "); put(file,q); new_line(file);
        t := Checker_Localization_Patterns.Transformation(n,integer32(q(fc)));
       -- put_line(file,"The pattern t for the numerical transformation ");
       -- put(file,t);
       -- put_line(file,"The numerical transformation :");
       -- Write_Standard_Moving_Flag(file,Numeric_Transformation(t));
       -- put_line(file,"The moving flag before the update :");
       -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,mf);
        start_mf := mf;
        mf := mf*Numeric_Transformation(t);
       -- put_line(file,"The moving flag after the update :");
       -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,mf);
        Checker_Homotopies.Define_Generalizing_Homotopy
          (file,n,q,qr,qc,stay_child,homtp,ctr);
        Initialize_Symbol_Table(n,k,q,qr,qc,dim);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay
            (file,n,k,ctr,q,p,qr,qc,pr,pc,verify,minrep,cond,
             mf,vf,sols,tol,fail,vrblvl-1);
        elsif homtp = 1 then
          Stay_Homotopy(file,n,k,ctr,nt,q,p,qr,qc,pr,pc,verify,minrep,
                        tosqr,cond,vf,mf,start_mf,sols,tol,fail,vrblvl-1);
        else -- homtp = 2
          Setup_Flag_Homotopies.Add_t_Symbol;
          Swap_Homotopy(file,n,k,ctr,nt,q,p,qr,qc,pr,pc,verify,minrep,
                        tosqr,cond,mf,start_mf,vf,sols,tol,fail,vrblvl-1);
        end if;
        if fail then
          put_line(file,"no longer a valid solution, abort tracking");
          new_line(file);
          unhappy := true; -- prevent from being concatenated
          exit;
        end if;
      end loop;
    end if;
  end Track_Path_in_Poset;

  procedure Track_Path_in_Poset
              ( n,k,nt : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in out Standard_Complex_Matrices.Matrix;
                start : in Standard_Complex_Solutions.Solution_List;
                sols : out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; unhappy : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    start_mf : Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    ptr,homtp,ctr,ind,fc : integer32;
    stay_child : boolean;
    fail : boolean := false;

    use Standard_Complex_Matrices;
    use Standard_Complex_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_Path_in_Poset 7 ...");
    end if;
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      unhappy := true;
    else
      if Is_Null(start) then
        return;
      end if;
      unhappy := false;
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      Copy(start,sols);
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        fc := Checker_Moves.Falling_Checker(p);
        t := Checker_Localization_Patterns.Transformation(n,integer32(q(fc)));
        start_mf := mf;
        mf := mf*Numeric_Transformation(t);
        Checker_Homotopies.Define_Generalizing_Homotopy
          (n,q,qr,qc,stay_child,homtp,ctr);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay(n,k,ctr,q,p,qr,qc,pr,pc,sols,fail,vrblvl-1);
        elsif homtp = 1 then
          Stay_Homotopy(n,k,ctr,nt,q,p,qr,qc,pr,pc,minrep,tosqr,cond,
                        vf,mf,start_mf,sols,tol,fail,rpt,vrblvl-1);
        else -- homtp = 2
          Swap_Homotopy(n,k,ctr,nt,q,p,qr,qc,pr,pc,minrep,tosqr,cond,
                        mf,start_mf,vf,sols,tol,fail,rpt,vrblvl-1);
        end if;
        if fail then
          unhappy := true; -- prevent from being concatenated
          exit;
        end if;
      end loop;
    end if;
  end Track_Path_in_Poset;

  procedure Track_Path_in_Poset
              ( n,k,nt : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf : in out DoblDobl_Complex_Matrices.Matrix;
                start : in DoblDobl_Complex_Solutions.Solution_List;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; unhappy : out boolean;
                vrblvl : in integer32 := 0 ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    start_mf : DoblDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    ptr,homtp,ctr,ind,fc : integer32;
    stay_child : boolean;
    fail : boolean := false;

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_Path_in_Poset 8 ...");
    end if;
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      unhappy := true;
    else
      if Is_Null(start) then
        return;
      end if;
      unhappy := false;
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      Copy(start,sols);
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        fc := Checker_Moves.Falling_Checker(p);
        t := Checker_Localization_Patterns.Transformation(n,integer32(q(fc)));
        start_mf := mf;
        mf := mf*Numeric_Transformation(t);
        Checker_Homotopies.Define_Generalizing_Homotopy
          (n,q,qr,qc,stay_child,homtp,ctr);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay(n,k,ctr,q,p,qr,qc,pr,pc,sols,fail,vrblvl-1);
        elsif homtp = 1 then
          Stay_Homotopy(n,k,ctr,nt,q,p,qr,qc,pr,pc,minrep,tosqr,cond,
                        vf,mf,start_mf,sols,tol,fail,vrblvl-1);
        else -- homtp = 2
          Swap_Homotopy(n,k,ctr,nt,q,p,qr,qc,pr,pc,minrep,tosqr,cond,
                        mf,start_mf,vf,sols,tol,fail,vrblvl-1);
        end if;
        if fail then
          unhappy := true; -- prevent from being concatenated
          exit;
        end if;
      end loop;
    end if;
  end Track_Path_in_Poset;

  procedure Track_Path_in_Poset
              ( n,k,nt : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf : in out QuadDobl_Complex_Matrices.Matrix;
                start : in QuadDobl_Complex_Solutions.Solution_List;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; unhappy : out boolean;
                vrblvl : in integer32 := 0 ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    start_mf : QuadDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    ptr,homtp,ctr,ind,fc : integer32;
    stay_child : boolean;
    fail : boolean := false;

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_Path_in_Poset 9 ...");
    end if;
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      unhappy := true;
    else
      if Is_Null(start) then
        return;
      end if;
      unhappy := false;
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      Copy(start,sols);
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        fc := Checker_Moves.Falling_Checker(p);
        t := Checker_Localization_Patterns.Transformation(n,integer32(q(fc)));
        start_mf := mf;
        mf := mf*Numeric_Transformation(t);
        Checker_Homotopies.Define_Generalizing_Homotopy
          (n,q,qr,qc,stay_child,homtp,ctr);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay(n,k,ctr,q,p,qr,qc,pr,pc,sols,fail,vrblvl-1);
        elsif homtp = 1 then
          Stay_Homotopy(n,k,ctr,nt,q,p,qr,qc,pr,pc,minrep,tosqr,cond,
                        vf,mf,start_mf,sols,tol,fail,vrblvl-1);
        else -- homtp = 2
          Swap_Homotopy(n,k,ctr,nt,q,p,qr,qc,pr,pc,minrep,tosqr,cond,
                        mf,start_mf,vf,sols,tol,fail,vrblvl-1);
        end if;
        if fail then
          unhappy := true; -- prevent from being concatenated
          exit;
        end if;
      end loop;
    end if;
  end Track_Path_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k,nt : in integer32; ps : in Poset;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                tol : in double_float;
                sols : out Standard_Complex_Solutions.Solution_List;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      mf : Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      ls : Link_to_Solution;
      fail : boolean;

    begin
      cnt := cnt + 1;
      Track_Path_in_Poset
        (file,n,k,ps,nds,cnt,verify,minrep,tosqr,cond,vf,mf,ls,tol,fail,
         rpt,vrblvl-1);
      if not fail
       then Append(sols,sols_last,ls.all);
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_All_Paths_in_Poset 1 ...");
    end if;
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k,nt : in integer32; ps : in Poset;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                tol : in double_float;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      mf : DoblDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      ls : Link_to_Solution;
      fail : boolean;

    begin
      cnt := cnt + 1;
      Track_Path_in_Poset
        (file,n,k,ps,nds,cnt,verify,minrep,tosqr,cond,vf,mf,ls,tol,fail,
         vrblvl-1);
      if not fail
       then Append(sols,sols_last,ls.all);
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_All_Paths_in_Poset 2 ...");
    end if;
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k,nt : in integer32; ps : in Poset;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                tol : in double_float;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      mf : QuadDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      ls : Link_to_Solution;
      fail : boolean;

    begin
      cnt := cnt + 1;
      Track_Path_in_Poset
        (file,n,k,ps,nds,cnt,verify,minrep,tosqr,cond,vf,mf,ls,tol,fail,
         vrblvl-1);
      if not fail
       then Append(sols,sols_last,ls.all);
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_All_Paths_in_Poset 3 ...");
    end if;
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k,nt : in integer32; ps : in Poset;
                child : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                tol : in double_float;
                start : in Standard_Complex_Solutions.Solution_List;
                sols : out Standard_Complex_Solutions.Solution_List;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      mf : Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      pp_sols : Solution_List; -- solutions on path in poset
      fail : boolean;
      leaf : constant Standard_Natural_Vectors.Vector := nds(nds'first).cols;

    begin
      cnt := cnt + 1;
      new_line(file);
      put(file,"Examining match for path "); put(file,cnt,1);
      put_line(file," ...");
      put(file,"-> leaf : "); put(file,leaf);
      put(file,"  child : "); put(file,child);
      if not Standard_Natural_Vectors.Equal(leaf,child) then
        put(file," no match, skip path "); put(file,cnt,1); new_line(file);
      else
        put(file," match at path "); put(file,cnt,1); new_line(file);
        Track_Path_in_Poset
          (file,n,k,nt,ps,nds,cnt,verify,minrep,tosqr,
           cond,vf,mf,start,pp_sols,tol,fail,rpt,vrblvl-1);
        if not fail
         then Concat(sols,sols_last,pp_sols);
        end if;
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_All_Paths_in_Poset 4 ...");
    end if;
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k,nt : in integer32; ps : in Poset;
                child : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                tol : in double_float;
                start : in DoblDobl_Complex_Solutions.Solution_List;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      mf : DoblDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      pp_sols : Solution_List; -- solutions on path in poset
      fail : boolean;
      leaf : constant Standard_Natural_Vectors.Vector := nds(nds'first).cols;

    begin
      cnt := cnt + 1;
      new_line(file);
      put(file,"Examining match for path "); put(file,cnt,1);
      put_line(file," ...");
      put(file,"-> leaf : "); put(file,leaf);
      put(file,"  child : "); put(file,child);
      if not Standard_Natural_Vectors.Equal(leaf,child) then
        put(file," no match, skip path "); put(file,cnt,1); new_line(file);
      else
        put(file," match at path "); put(file,cnt,1); new_line(file);
        Track_Path_in_Poset
          (file,n,k,nt,ps,nds,cnt,verify,minrep,tosqr,
           cond,vf,mf,start,pp_sols,tol,fail,vrblvl-1);
        if not fail
         then Concat(sols,sols_last,pp_sols);
        end if;
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_All_Paths_in_Poset 5 ...");
    end if;
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k,nt : in integer32; ps : in Poset;
                child : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                tol : in double_float;
                start : in QuadDobl_Complex_Solutions.Solution_List;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      mf : QuadDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      pp_sols : Solution_List; -- solutions on path in poset
      fail : boolean;
      leaf : constant Standard_Natural_Vectors.Vector := nds(nds'first).cols;

    begin
      cnt := cnt + 1;
      new_line(file);
      put(file,"Examining match for path "); put(file,cnt,1);
      put_line(file," ...");
      put(file,"-> leaf : "); put(file,leaf);
      put(file,"  child : "); put(file,child);
      if not Standard_Natural_Vectors.Equal(leaf,child) then
        put(file," no match, skip path "); put(file,cnt,1); new_line(file);
      else
        put(file," match at path "); put(file,cnt,1); new_line(file);
        Track_Path_in_Poset
          (file,n,k,nt,ps,nds,cnt,verify,minrep,tosqr,
           cond,vf,mf,start,pp_sols,tol,fail,vrblvl-1);
        if not fail
         then Concat(sols,sols_last,pp_sols);
        end if;
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_All_Paths_in_Poset 6 ...");
    end if;
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( n,k,nt : in integer32; ps : in Poset;
                child : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                tol : in double_float;
                start : in Standard_Complex_Solutions.Solution_List;
                sols : out Standard_Complex_Solutions.Solution_List;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      mf : Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      pp_sols : Solution_List; -- solutions on path in poset
      fail : boolean;
      leaf : constant Standard_Natural_Vectors.Vector := nds(nds'first).cols;

    begin
      cnt := cnt + 1;
      if Standard_Natural_Vectors.Equal(leaf,child) then
        Track_Path_in_Poset
          (n,k,nt,ps,nds,cnt,minrep,tosqr,cond,vf,mf,start,pp_sols,tol,fail,
           rpt,vrblvl-1);
        if not fail
         then Concat(sols,sols_last,pp_sols);
        end if;
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_All_Paths_in_Poset 7 ...");
    end if;
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( n,k,nt : in integer32; ps : in Poset;
                child : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                tol : in double_float;
                start : in DoblDobl_Complex_Solutions.Solution_List;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      mf : DoblDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      pp_sols : Solution_List; -- solutions on path in poset
      fail : boolean;
      leaf : constant Standard_Natural_Vectors.Vector := nds(nds'first).cols;

    begin
      cnt := cnt + 1;
      if Standard_Natural_Vectors.Equal(leaf,child) then
        Track_Path_in_Poset
          (n,k,nt,ps,nds,cnt,minrep,tosqr,cond,vf,mf,start,pp_sols,tol,fail,
           vrblvl-1);
        if not fail
         then Concat(sols,sols_last,pp_sols);
        end if;
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_All_Paths_in_Poset 8 ...");
    end if;
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( n,k,nt : in integer32; ps : in Poset;
                child : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                tol : in double_float;
                start : in QuadDobl_Complex_Solutions.Solution_List;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      mf : QuadDobl_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      pp_sols : Solution_List; -- solutions on path in poset
      fail : boolean;
      leaf : constant Standard_Natural_Vectors.Vector := nds(nds'first).cols;

    begin
      cnt := cnt + 1;
      if Standard_Natural_Vectors.Equal(leaf,child) then
        Track_Path_in_Poset
          (n,k,nt,ps,nds,cnt,minrep,tosqr,cond,vf,mf,start,pp_sols,tol,fail,
           vrblvl-1);
        if not fail
         then Concat(sols,sols_last,pp_sols);
        end if;
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    if vrblvl > 0 then
      put("-> in checker_poset_deformations.");
      put_line("Track_All_Paths_in_Poset 9 ...");
    end if;
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

end Checker_Poset_Deformations;
