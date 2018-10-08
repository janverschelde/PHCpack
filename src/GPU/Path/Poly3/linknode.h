// The file linknode.h defines the classes MonData and LinkNode.

#ifndef __LINKNODE_H__
#define __LINKNODE_H__

#include <iostream>

using namespace std;

// Locate the monomial from the string
/*
   Recording the monomial start & end positions, and coefficent.
 */

template <class T>
class MonData
{
   public:

      int start; // start position in the string
      int end;   // end position in the string
      T coef;    // coefficient of the monomial

      MonData ( int s, int e, T c ) // constructor
      /*
         @param s an integer argument.
         @param e a constant character pointer.
       */
      {
         start = s;
         end = e;
         coef = c;
      }

      void print() // prints the data in a MonData object
      {
         cout << start << " " << end << " " << coef << endl;
      }
};

// Linked node
/*
   LinkNode to be used in the LinkList
 */

template <class T>
class LinkNode
{
   public:

      T data; // data in the node
      LinkNode* next; // address of the next linked node

      LinkNode ( T new_data ) // constructor
      /*
         Creation of a new node of new_data.
         @param new_data a T type of data
       */
      {
         data = new_data;
         next = NULL;
      }

      LinkNode* append ( T new_data );
      // A constructor function
      /*
         Creates a new node of new_data and link to the current node.
         @param new_data a T type of data
         @return Address of the new node
         @sa LinkNode()
       */

      void print();
      void print_class();
};

#include "linknode.tpp"

#endif
