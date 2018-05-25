with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;
with Standard_Complex_VecVecs;
with Lists_of_Integer_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_SysFun;
with Exponent_Vectors;
with Standard_Complex_Laur_JacoMats;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Random_Coefficient_Systems;
with Lists_of_Strings;
with DEMiCs_Command_Line;
with DEMiCs_Output_Convertors;
with DEMiCs_Algorithm;                   use DEMiCs_Algorithm;
with DEMiCs_Output_Data;
with Pipelined_Cell_Indices;
with Semaphore;
with Multitasking;
with Polyhedral_Start_Systems;          use Polyhedral_Start_Systems;
with Pipelined_Cell_Trackers;           use Pipelined_Cell_Trackers;

package body Pipelined_Polyhedral_Homotopies is

  function Random_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return Standard_Floating_VecVecs.Link_to_VecVec is

    res : Standard_Floating_VecVecs.Link_to_VecVec;
    resrep : Standard_Floating_VecVecs.VecVec(mix'range);
    idx : integer32 := 1;
    len : integer32;

  begin
    for i in resrep'range loop
      len := integer32(Lists_of_Integer_Vectors.Length_Of(sup(idx)));
      declare
        vals : Standard_Floating_Vectors.Vector(1..len);
      begin
        for j in 1..len loop
          vals(j) := Standard_Random_Numbers.Random;
        end loop;
        resrep(i) := new Standard_Floating_Vectors.Vector'(vals);
      end;
      idx := idx + mix(i);
    end loop;
    res := new Standard_Floating_VecVecs.VecVec'(resrep);
    put_line("The random lifting : ");
    Standard_Floating_VecVecs_io.put(res.all);
    return res;
  end Random_Lifting;

  procedure Pipeline_Cells_to_Paths
              ( dim,nt : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                verbose : in boolean := true ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;

    lif : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Random_Lifting(mix,sup);
    lsp : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
        := DEMiCs_Output_Convertors.Apply_Lifting(mix.all,sup,lif.all);
    nbequ : constant integer32 := q'last;
    r : constant integer32 := mix'last;
    hom : Eval_Coeff_Laur_Sys(1..nbequ);
    cff : Standard_Complex_VecVecs.VecVec(hom'range);
    epv : Exponent_Vectors.Exponent_Vectors_Array(hom'range);
    ejf : Eval_Coeff_Jaco_Mat(hom'range,hom'first..hom'last+1);
    jmf : Mult_Factors(ejf'range(1),ejf'range(2));
    tmv : Standard_Natural_Vectors.Vector(2..nt) := (2..nt => 0);
    dpw : Standard_Floating_VecVecs.Array_of_VecVecs(2..nt);
    cft : Standard_Complex_VecVecs.Array_of_VecVecs(2..nt);
    tasksols,lastsols : Array_of_Solution_Lists(2..nt);
    sem : Semaphore.Lock;
    
    procedure Track ( idtask : in integer32; 
                      mtype : in Standard_Integer_Vectors.Link_to_Vector;
                      mic : in Mixed_Cell ) is
    begin
      Standard_Track_Cell(sem,idtask,nbequ,r,mix.all,mic,lsp,cff,
        dpw(idtask),cft(idtask),epv,hom,ejf,jmf,q,tmv(idtask),
        tasksols(idtask),lastsols(idtask));
    end Track;

  begin
    q := Random_Coefficient_Systems.Create(natural32(dim),sup);
    DEMiCs_Output_Data.allocate := true;
    DEMiCs_Output_Data.Store_Dimension_and_Mixture(dim,mix);
    DEMiCs_Output_Data.Initialize_Allocated_Cell_Pointer;
    cff := Coeff(q);
    epv := Exponent_Vectors.Create(q);
    hom := Create(q);
    Create(q,ejf,jmf);
    Allocate_Workspace_for_Exponents(epv,dpw);
    Allocate_Workspace_for_Coefficients(cff,cft);
    Pipelined_Cell_Indices.Pipelined_Mixed_Cells
      (nt,dim,mix,sup,lif,lsp,Track'access,verbose);
    for k in tasksols'range loop
      Standard_Complex_Solutions.Push(tasksols(k),qsols);
    end loop;
  end Pipeline_Cells_to_Paths;

end Pipelined_Polyhedral_Homotopies;
