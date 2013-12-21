with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Standard_Integer64_Vectors_io;       use Standard_Integer64_Vectors_io;
with Standard_Lattice_4d_Facets_io;

package body Standard_Lattice_Facets_io is

  procedure Write_Facet
              ( A : in Standard_Integer64_Matrices.Matrix;
                f : in Standard_Lattice_Facets.Link_to_Facet ) is
 
    use Standard_Lattice_Facets;
    use Standard_Lattice_4d_Facets_io;
  
  begin
    put("The facet normal :"); put(f.normal);
    put(", #ridges : "); put(f.m,1); new_line;
    put("  supported by "); put(f.n,1);
    put(" points :"); put(f.points); new_line;
    if f.d = 4 then
      for i in f.ridges'range loop
        put("ridge "); put(i,1); put(" with normal");
        put(f.ridges(i).normal);
        put(" supports"); put(f.ridges(i).points); -- new_line;
        Check_Inner_Products(A,f.ridges(i).normal,f.ridges(i).points);
      end loop;
    else
      for i in f.faces'range loop
        put("face "); put(i,1); put(" with normal");
        put(f.faces(i).normal);
        put(" supports"); put(f.faces(i).points); -- new_line;
        Check_Inner_Products(A,f.faces(i).normal,f.faces(i).points);
      end loop;
    end if;
    for j in f.neighbors'range loop
      if f.neighbors(j) /= null then
        put("-> ridge "); put(j,1);
        put(" is shared with facet ");
        put(f.neighbors(j).label,1); new_line;
      end if;
    end loop;
  end Write_Facet;

  procedure Write_Facets
              ( A : in Standard_Integer64_Matrices.Matrix;
                f : in Standard_Lattice_Facets.Facet_List ) is

    use Standard_Lattice_Facets;
    tmp : Facet_List := f;
    lft : Link_to_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      put("*** facet with label "); put(lft.label,1); put_line(" ***");
      Write_Facet(A,lft);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Facets;

end Standard_Lattice_Facets_io;
