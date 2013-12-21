with Localization_Posets;               use Localization_Posets;

function Pieri_Count ( m,p,q : integer ) return integer is

  root : Node(p) := Trivial_Root(m,p,q);
  lnkroot : Link_to_Node := new Node'(root);
  nq : constant natural := m*p + q*(m+p);
  level_poset : Array_of_Nodes(0..nq);

begin
  Q_Top_Bottom_Create(lnkroot,root.bottom(p),m+p);
  level_poset := Create_Leveled_Poset(lnkroot);
  Count_Roots(level_poset);
  return level_poset(level_poset'last).roco;
end Pieri_Count;
