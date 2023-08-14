with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Double_Double_Vectors;
with Double_Double_Vectors_io;           use Double_Double_Vectors_io;

package body Double_Double_Theta_Algorithm is

  procedure Allocate
              ( tab : out Double_Double_VecVecs.VecVec;
                idx : out Standard_Integer_Vectors.Vector;
                dim : in integer32 ) is

    zero : constant double_double := create(0.0);

  begin
    for i in idx'range loop
      idx(i) := -1;
    end loop;
    for i in tab'range loop
       declare
         col : constant Double_Double_Vectors.Vector(0..dim)
             := (0..dim => zero);
         use Double_Double_Vectors;
       begin
         if tab(i) = null then
           tab(i) := new Double_Double_Vectors.Vector'(col);
         else -- the dimension may not match, must clear first
           Double_Double_Vectors.Clear(tab(i));
           tab(i) := new Double_Double_Vectors.Vector'(col);
         end if;
       end;
    end loop;
  end Allocate;

  procedure Initialize
              ( tab : in Double_Double_VecVecs.VecVec;
                dim : in out integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                nbr : in double_double;
                verbose : in boolean := true ) is

    one : constant double_double := Double_Double_Numbers.create(1.0);
    d1 : double_double;

  begin
    if dim = 0 then
      dim := 1;
      idx(0) := 0;
      tab(0)(idx(0)) := nbr;
      if verbose then
        put_line("theta 0 : "); put_line(tab(0)(0..0));
      end if;
    else
      idx(0) := idx(0) + 1;
      tab(0)(idx(0)) := nbr;
      if verbose then
        put_line("theta 0 : "); put_line(tab(0)(0..idx(0)));
      end if;
      if dim = 1 then
        dim := 2;
        d1 := one/(nbr - tab(0)(idx(0)-1));
        idx(1) := 0;
        tab(1)(idx(1)) := d1;
        if verbose then
           put_line("theta 1 : "); put_line(tab(1)(0..0));
        end if;
      else
        d1 := one/(nbr - tab(0)(idx(0)-1));
        idx(1) := idx(1) + 1;
        tab(1)(idx(1)) := d1;
        if verbose then
           put_line("theta 1 : "); put_line(tab(1)(0..idx(1)));
        end if;
      end if;
    end if;
  end Initialize;

  procedure Columns 
              ( tab : in Double_Double_VecVecs.VecVec;
                dim : in integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                nbr : in double_double;
                verbose : in boolean := true ) is

    one : constant double_double := Double_Double_Numbers.create(1.0);
    two : constant double_double := Double_Double_Numbers.create(2.0);
    d0,d1,d2,newtheta : double_double;
    prvcol1,prvcol2 : Double_Double_Vectors.Link_to_Vector;

  begin
    if verbose then
      put("in Columns, dim = "); put(dim,1); put_line(" ...");
    end if;
    idx(0) := idx(0) + 1;
    tab(0)(idx(0)) := nbr;
    if verbose then
      put_line("theta 0 : "); put_line(tab(0)(0..idx(0)));
    end if;
    d1 := one/(nbr - tab(0)(idx(0)-1));
    idx(1) := idx(1) + 1;
    tab(1)(idx(1)) := d1;
    if verbose then
      put_line("theta 1 : "); put_line(tab(1)(0..idx(1)));
    end if;
    for col in 2..dim-1 loop
      prvcol1 := tab(col-1);
      prvcol2 := tab(col-2);
      idx(col) := idx(col) + 1;
      if col mod 2 = 1 then
        d0 := prvcol1(idx(col)+1) - prvcol1(idx(col));
        newtheta := prvcol2(idx(col)+1) + one/d0;
      else
        d0 := prvcol2(idx(col)+2) - prvcol2(idx(col)+1);
        d1 := prvcol1(idx(col)+2) - prvcol1(idx(col)+1);
        d2 := prvcol1(idx(col)+2) - two*prvcol1(idx(col)+1)
            + prvcol1(idx(col));
        newtheta := prvcol2(idx(col)+1) + (d0*d1)/d2;
      end if;
      tab(col)(idx(col)) := newtheta;
      if verbose then
        put("theta "); put(col,1); put_line(" : ");
        put_line(tab(col)(0..idx(col)));
      end if;
    end loop;
  end Columns;

  procedure New_Column
              ( tab : in Double_Double_VecVecs.VecVec;
                dim : in out integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    one : constant double_double := Double_Double_Numbers.create(1.0);
    two : constant double_double := Double_Double_Numbers.create(2.0);
    d0,d1,d2,newtheta : double_double;
    lastcol1,lastcol2 : Double_Double_Vectors.Link_to_Vector;

  begin
    if verbose then
      put("in New_Column, dim = "); put(dim,1);
      put("  idx = "); put(idx(0..dim-1)); new_line;
    end if;
    if dim mod 2 = 0 then
      if idx(dim-1) < 2 then
        if verbose
         then put_line("too few elements to add even indexed column");
        end if;
      else
        if verbose
         then put_line("adding even indexed column");
        end if;
        lastcol1 := tab(dim-1);
        lastcol2 := tab(dim-2);
        d0 := lastcol2(2) - lastcol2(1);
        d1 := lastcol1(2) - lastcol1(1);
        d2 := lastcol1(2) - two*lastcol1(1) + lastcol1(0);
        newtheta := lastcol2(1) + (d0*d1)/d2;
        idx(dim) := 0;
        tab(dim)(idx(dim)) := newtheta;
        if verbose then
          put("theta "); put(dim,1); put_line(" : ");
          put_line(tab(dim)(0..idx(dim)));
        end if;
        dim := dim + 1;
      end if;
    else
      if idx(dim-1) < 1 then
        if verbose
         then put_line("too few elements to add odd indexed column");
        end if;
      else
        if verbose
         then put_line("adding odd indexed column");
        end if;
        lastcol1 := tab(dim-1);
        lastcol2 := tab(dim-2);
        d0 := lastcol1(1) - lastcol1(0);
        newtheta := lastcol2(1) + one/d0;
        idx(dim) := 0;
        tab(dim)(idx(dim)) := newtheta;
        if verbose then
          put("theta "); put(dim,1); put_line(" : ");
          put_line(tab(dim)(0..idx(dim)));
        end if;
        dim := dim + 1;
      end if;
    end if;
  end New_Column;

  procedure Extrapolate
              ( tab : in Double_Double_VecVecs.VecVec;
                dim : in out integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                nbr : in double_double;
                verbose : in boolean := true ) is
  begin
    if dim = 0 then
      Initialize(tab,dim,idx,nbr,verbose);
    elsif dim = 1 and idx(0) < 1 then
      Initialize(tab,dim,idx,nbr,verbose);
    else
      Columns(tab,dim,idx,nbr,verbose);
      New_Column(tab,dim,idx,verbose);
    end if;
  end Extrapolate;

end Double_Double_Theta_Algorithm;
