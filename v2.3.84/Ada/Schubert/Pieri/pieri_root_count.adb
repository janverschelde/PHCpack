with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Localization_Posets;               use Localization_Posets;

function Pieri_Root_Count ( m,p,q : natural32 ) return natural32 is

  res : natural32;
  root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
  lnkroot : Link_to_Node := new Node'(root);
  nq : constant natural32 := m*p + q*(m+p);
  level_poset : Array_of_Nodes(0..integer32(nq));

begin
  Q_Top_Bottom_Create(lnkroot,root.bottom(integer32(p)),m+p);
  level_poset := Create_Leveled_Poset(lnkroot);
  Count_Roots(level_poset);
  res := natural32(level_poset(level_poset'last).roco);
  Clear(level_poset);
  return res;
end Pieri_Root_Count;
