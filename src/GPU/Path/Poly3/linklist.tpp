// The file linklist.tpp defines the methods with prototypes in linklist.h.

template <class T>
void LinkList<T>::append(T data)
{
   link_last = link_last->append(data);
   n_node++;
}

template <class T>
void LinkList<T>::print()
{
   LinkNode<T>* link_tmp = link_header;
   link_tmp = link_tmp->next;
   while(link_tmp)
   {
      cout << link_tmp->data << endl;
      link_tmp = link_tmp->next;
   }
   cout << endl;
}

template <class T>
void LinkList<T>::print_class()
{
   LinkNode<T>* link_tmp = link_header;
   link_tmp = link_tmp->next;
   while(link_tmp)
   {
      link_tmp->data->print();
      link_tmp = link_tmp->next;
   }
   cout << endl;
}

template <class T>
void LinkList<T>::destroy()
{
   LinkNode<T>* item = link_header;
   while(item)
   {
      LinkNode<T>* old = item;
      item = item->next;
      delete old;
   }
}

template <class T>
void LinkList<T>::destroy_class()
{
   LinkNode<T>* item = link_header;
   while(item)
   {
      LinkNode<T>* old = item;
      item = item->next;
      delete old->data;
      delete old;
   }
}
