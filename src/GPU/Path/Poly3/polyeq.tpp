template <class ComplexType, class RealType>
ComplexType PolyEq<ComplexType,RealType>::eval ( const ComplexType* x_val )
{
   ComplexType val = constant;

   for(int i=0; i<n_mon; i++) val += mon[i]->eval(x_val);

   return val;
}

template <class ComplexType, class RealType>
ComplexType PolyEq<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* deri )
{
   for(int i=0; i<dim; i++)
   {
      deri[i].init(0.0,0.0);
   }
   ComplexType val = constant;
   // std::cout << constant << std::endl;
    
   ComplexType* mon_deri = new ComplexType[dim];

   for(int i=0; i<n_mon; i++)
   {
      PolyMon<ComplexType,RealType>* m = mon[i];
      val += m->eval(x_val, mon_deri);

      for(int j=0; j<m->n_var; j++)
      {
         deri[m->pos[j]] += mon_deri[j];
      }
   }
   delete [] mon_deri;

   return val;
}

template <class ComplexType, class RealType>
ComplexType PolyEq<ComplexType,RealType>::eval
 ( const ComplexType* x_val, ComplexType* deri, ComplexType** deg_table )
{
   for(int i=0; i<dim; i++)
   {
      deri[i].init(0.0,0.0);
   }
   ComplexType val = constant;
   // std::cout << constant << std::endl;

   ComplexType* mon_deri = new ComplexType[dim];
   for(int mon_idx=0; mon_idx<n_mon; mon_idx++)
   {
      // std::cout << "mon " << mon_idx << std::endl;
      PolyMon<ComplexType,RealType>* m = mon[mon_idx];
      val += m->eval(x_val, mon_deri, deg_table);
      for(int j=0; j<m->n_var; j++)
      {
         deri[m->pos[j]] += mon_deri[j];
      }
   }
   delete [] mon_deri;

   return val;
}
