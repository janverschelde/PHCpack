with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;

package body Double_Rho_Algorithm is

  procedure Allocate
              ( tab : out Standard_Floating_VecVecs.VecVec;
                idx : out Standard_Integer_Vectors.Vector;
                dim : in integer32 ) is
  begin
    for i in idx'range loop
      idx(i) := -1;
    end loop;
    for i in tab'range loop
       declare
         col : constant Standard_Floating_Vectors.Vector(0..dim)
             := (0..dim => 0.0);
         use Standard_Floating_Vectors;
       begin
         if tab(i) = null then
           tab(i) := new Standard_Floating_Vectors.Vector'(col);
         else -- the dimension may not match, must clear first
           Standard_Floating_Vectors.Clear(tab(i));
           tab(i) := new Standard_Floating_Vectors.Vector'(col);
         end if;
       end;
    end loop;
  end Allocate;

  procedure Initialize
              ( tab : in Standard_Floating_VecVecs.VecVec;
                dim : in out integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                nbr : in double_float;
                verbose : in boolean := true ) is

    d1 : double_float;

  begin
    if dim = 0 then
      dim := 1;
      idx(0) := 0;
      tab(0)(idx(0)) := nbr;
      if verbose then
        put_line("rho 0 : "); put_line(tab(0)(0..0));
      end if;
    else
      idx(0) := idx(0) + 1;
      tab(0)(idx(0)) := nbr;
      if verbose then
        put_line("rho 0 : "); put_line(tab(0)(0..idx(0)));
      end if;
      if dim = 1 then
        dim := 2;
        d1 := 1.0/(nbr - tab(0)(idx(0)-1));
        idx(1) := 0;
        tab(1)(idx(1)) := d1;
        if verbose then
           put_line("rho 1 : "); put_line(tab(1)(0..0));
        end if;
      else
        d1 := 1.0/(nbr - tab(0)(idx(0)-1));
        idx(1) := idx(1) + 1;
        tab(1)(idx(1)) := d1;
        if verbose then
           put_line("rho 1 : "); put_line(tab(1)(0..idx(1)));
        end if;
      end if;
    end if;
  end Initialize;

  procedure Columns 
              ( tab : in Standard_Floating_VecVecs.VecVec;
                dim : in integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                nbr : in double_float;
                verbose : in boolean := true ) is

    d1,num,den,invrho1 : double_float;
    prvcol1,prvcol2 : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if verbose then
      put("in Columns, dim = "); put(dim,1); put_line(" ...");
    end if;
    idx(0) := idx(0) + 1;
    tab(0)(idx(0)) := nbr;
    if verbose then
      put_line("rho 0 : "); put_line(tab(0)(0..idx(0)));
    end if;
    d1 := 1.0/(nbr - tab(0)(idx(0)-1));
    idx(1) := idx(1) + 1;
    tab(1)(idx(1)) := d1;
    if verbose then
      put_line("rho 1 : "); put_line(tab(1)(0..idx(1)));
    end if;
    for col in 2..dim-1 loop
      prvcol1 := tab(col-1);
      prvcol2 := tab(col-2);
      num := prvcol2(idx(col-2)-1);
      den := prvcol1(idx(col-1)) - prvcol1(idx(col-1)-1);
      invrho1 := double_float(col)/den;
      idx(col) := idx(col) + 1;
      tab(col)(idx(col)) := num + invrho1;
      if verbose then
         put("rho "); put(col,1); put_line(" : ");
         put_line(tab(col)(0..idx(col)));
       end if;
    end loop;
  end Columns;

  procedure New_Column
              ( tab : in Standard_Floating_VecVecs.VecVec;
                dim : in out integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    d0,newrho : double_float;
    lastcol1,lastcol2 : Standard_Floating_Vectors.Link_to_Vector;

  begin
    lastcol1 := tab(dim-1);
    lastcol2 := tab(dim-2);
    d0 := lastcol1(1) - lastcol1(0);
    newrho := lastcol2(1) + double_float(dim)/d0;
    idx(dim) := 0;
    tab(dim)(idx(dim)) := newrho;
    if verbose then
       put("rho "); put(dim,1); put_line(" : ");
       put_line(tab(dim)(0..idx(dim)));
     end if;
    dim := dim + 1;
  end New_Column;

  procedure Extrapolate
              ( tab : in Standard_Floating_VecVecs.VecVec;
                dim : in out integer32;
                idx : in out Standard_Integer_Vectors.Vector;
                nbr : in double_float;
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

end Double_Rho_Algorithm;
