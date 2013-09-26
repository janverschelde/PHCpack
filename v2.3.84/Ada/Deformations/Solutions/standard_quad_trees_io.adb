with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;

package body Standard_Quad_Trees_io is

  procedure put ( root : Link_to_Quad_Node ) is
  begin
    put(standard_output,root);
  end put;

  procedure put ( file : file_type; root : Link_to_Quad_Node ) is
  begin
    for i in 1..root.depth loop
      put(file,"  ");
    end loop;
    if root.depth > 0
     then put(file,"|-- ");
    end if;
    put(file,root.size,1);
    new_line(file);
    if not root.leaf then
      put(file,root.ne);
      put(file,root.nw);
      put(file,root.sw);
      put(file,root.se);
    end if;
  end put;

end Standard_Quad_Trees_io;
