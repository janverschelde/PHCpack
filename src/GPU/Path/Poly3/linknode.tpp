// The file linknode.tpp defines the methods with prototypes in linknode.h.

template <class T>
LinkNode<T>* LinkNode<T>::append(T data)
{
   LinkNode<T>* newnode = new LinkNode<T>(data);
   next = newnode;
   return newnode;
}

template <class T>
void LinkNode<T>::print()
{
   cout << data << endl;
}
