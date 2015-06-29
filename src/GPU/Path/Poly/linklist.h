#ifndef LINKLIST_H
#define LINKLIST_H

#include <iostream>

using namespace std;

/// Locate the monomial from the string
/**
Recording the monomial start & end positions, and coefficent.
*/
template <class T>
class MonData{
public:
	/// start position in the string
	int start;
	/// end position in the string
	int end;
	/// coefficient of the monomial
	T coef;


	/// A constructor.
	/**
	@param s an integer argument.
	@param e a constant character pointer.
	*/
	MonData(int s, int e, T c){
		start = s;
		end = e;
		coef = c;
	}

	/// A print function
	/**    
	Print the data of MonData instant. 
	*/
	void print(){
		cout << start << " " << end << " " << coef << endl;
	}
};

/// Linked node
/**
LinkNode to be used in the LinkList
*/
template <class T>
class LinkNode{
public:
	/// Data in the node
	T data;
	/// Address of the next linked node
	LinkNode* next;

	/// A constructor
	/**
	Create of a new node of new_data.
	@param new_data a T type of data
	*/
	LinkNode(T new_data){
		data = new_data;
		next = NULL;
	}

	/// A constructor function
	/**
	Create a new node of new_data and link to the current node.
	@param new_data a T type of data
	@return Address of the new node
	@sa LinkNode()
	*/
	LinkNode* append(T new_data);
	void print();
	void print_class();
};

/// linked list based on LinkNode
/**
A linked list based on LinkNode
*/

template <class T>
class LinkList{
public:
	/// Header of the LinkList
	LinkNode<T>* link_header;

	/// Last node of the LinkList
	LinkNode<T>* link_last;

	/// Number of Nodes in the LinkList
	int n_node;

	/// A constructor
	/**
	Create a zero node as the header;
	*/
	LinkList(){
		link_header = new LinkNode<T>(0);
		link_last = link_header;
		n_node = 0;
	}
	
	LinkNode<T>* header(){
		return link_header->next;
	}

	/// Append a new node with new_data
	/**
	Use LinkNode::append to creat a new node of new_data
	and append it to the NodeList. 
	Increase n_node by 1;
	/param new_data a data of T type
	/sa LinkNode::append
	*/
	void append(T new_data);

	/// Print the list
	/**
	Print the data of the list. Header is skipped.
	*/
	void print();

	/// Print the list of a certain class
	/**
	Use the print method in the class to print the list.
	Header is skipped.
	*/
	void print_class();

	/// Destructor for linked list of simple data
	/**
	Delete every node in the list
	*/
	void destroy();
	/// Destructor for linked list of certain class
	/**
	Delete every node and its data in the list
	*/
	void destroy_class();
};

#include "linklist.tpp"

#endif
