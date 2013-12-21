with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Multprec_Integer_Numbers;            use Multprec_Integer_Numbers;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Multprec_Integer_Vectors_io;         use Multprec_Integer_Vectors_io;
with Multprec_Lattice_Supports;

package body Multprec_Lattice_4d_Facets_io is

  procedure Check_Inner_Products
              ( A : in Multprec_Integer_Matrices.Matrix;
                v : in Multprec_Integer_Vectors.Vector;
                s : in Standard_Integer_Vectors.Vector ) is

    p : Multprec_Integer_Vectors.Vector(A'range(2))
      := Multprec_Lattice_Supports.Inner_Products(A,v);

  begin
    put("  IP :"); put(p);
    for i in s'first+1..s'last loop
      if not Equal(p(s(s'first)),p(s(i)))
       then put_line(" NOT okay!"); return;
      end if;
    end loop;
    put_line(" OK");
    Multprec_Integer_Vectors.Clear(p);
  end Check_Inner_Products;

  procedure Write_4D_Facet
              ( A : in Multprec_Integer_Matrices.Matrix;
                f : in Multprec_Lattice_4D_Facets.Link_to_4d_Facet ) is
 
    use Multprec_Lattice_4D_Facets;
  
  begin
    put("The facet normal :"); put(f.normal);
    put(", #ridges : "); put(f.m,1); new_line;
    put("  supported by "); put(f.n,1);
    put(" points :"); put(f.points); new_line;
    for i in f.ridges'range loop
      put("ridge "); put(i,1); put(" with normal");
      put(f.ridges(i).normal);
      put(" supports"); put(f.ridges(i).points); -- new_line;
      Check_Inner_Products(A,f.ridges(i).normal,f.ridges(i).points);
    end loop;
    for j in f.neighbors'range loop
      if f.neighbors(j) /= null then
        put("-> ridge "); put(j,1);
        put(" is shared with facet ");
        put(f.neighbors(j).label,1); new_line;
      end if;
    end loop;
  end Write_4D_Facet;

  procedure Write_4D_Facets
              ( A : in Multprec_Integer_Matrices.Matrix;
                f : in Multprec_Lattice_4D_Facets.Facet_4D_List ) is

    use Multprec_Lattice_4D_Facets;
    tmp : Facet_4D_List := f;
    lft : Link_to_4d_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      put("*** facet with label "); put(lft.label,1); put_line(" ***");
      Write_4D_Facet(A,lft);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_4D_Facets;

end Multprec_Lattice_4d_Facets_io;
