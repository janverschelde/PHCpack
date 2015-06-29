
#ifndef VARSET_H_
#define VARSET_H_

#include <iostream>
#include<algorithm>

#include "DefineType_Host.h"

#define IntSet GeneralSet<int>
#define ShortSet GeneralSet<short>
#define EqSet GeneralSet<EqIdxCoef>

class VarSet{
public:
	int n;
	VarSet* elements;

	void init(){
		this->n = 0;
		elements = NULL;

	}
	void init(const int n){
		this->n = n;
		elements = NULL;
	}

	void init(const int n, const VarSet* elements){
		this->n = n;
		if(elements != NULL){
			this->elements = new VarSet[n];
			for(int i=0; i<n; i++){
				this->elements[i] = elements[i];
			}
		}
		else{
			this->elements = NULL;
		}
	};

	VarSet(){
		init();
	}

	VarSet(int n){
		init(n);
	}

	VarSet(const int n, const VarSet* elements){
		init(n, elements);
	}

	VarSet(const VarSet& original){
		init(original.n, original.elements);
	}

	VarSet& operator= (const VarSet& original){
		init(original.n, original.elements);
		return *this;
	}

	bool operator < (const VarSet& that) const{
		if(this->n < that.n){
			return true;
		}
		else if(this->n > that.n){
			return false;
		}
		else if(this->n > 0){
			// Check if one is empty
			if(this->elements == NULL && that.elements == NULL){
				return false;
			}
			else if(this->elements == NULL) {
				std::cout << "Not on the same level" << std::endl;
				return false;
			}
			else if(that.elements == NULL) {
				std::cout << "Not on the same level" << std::endl;
				return false;
			}

			// Compare each element
			for(int i=0; i<n; i++){
				if(this->elements[i] < that.elements[i]){
					return true;
				}
				else if(this->elements[i] > that.elements[i]){
					return false;
				}
			}
			return false;
		}
		else{
			return false;
		}
	}

	bool operator > (const VarSet& that) const{
		if(this->n > that.n){
			return true;
		}
		else if(this->n < that.n){
			return false;
		}
		else if(this->n > 0){
			// Check if one is empty
			if(this->elements == NULL && that.elements == NULL){
				return false;
			}
			else if(this->elements == NULL) {
				std::cout << "Not on the same level" << std::endl;
				return false;
			}
			else if(that.elements == NULL) {
				std::cout << "Not on the same level" << std::endl;
				return false;
			}

			for(int i=0; i<n; i++){
				if(this->elements[i] > that.elements[i]){
					return true;
				}
				else if(this->elements[i] < that.elements[i]){
					return false;
				}
			}
			return false;
		}
		else{
			return false;
		}
	}

	~VarSet(){
		delete[] elements;
		elements = NULL;
	}

	void print(int level = 0){
		if(elements != NULL){
			std::cout << std::endl;
			for(int i=0; i<level; i++){
				std::cout << "   ";
			}
			std::cout << "l"<< level << " n" << n << " : ";
			for(int i=0; i<n; i++){
				elements[i].print(level+1);
			}
		}
		else{
			std::cout << n << ", ";
		}
		if(level == 0){
			std::cout << std::endl;
		}
	}

	void sorted(){
		if(elements != NULL){
			for(int i=0; i<n; i++){
				elements[i].sorted();
			}
			std::sort(elements, elements+n);
		}
	}
};

template<class T>
class GeneralSet{
	friend class MonSet;
	int n;
	T* elements;

	void init(const int n, const T* elements){
		this->n = n;
		this->elements = new T[n];
		for(int i=0; i<n; i++){
			this->elements[i] = elements[i];
		}
	}

public:
	GeneralSet(){
		n = 0;
		elements = NULL;
	};

	GeneralSet(const int n, const T* elements){
		init(n, elements);
	}

	GeneralSet& operator= (const GeneralSet& original){
		init(original.n, original.elements);
		return *this;
	}

	~GeneralSet(){
		delete[] elements;
		elements = NULL;
	}

	int get_n(){
		return n;
	}

	int get_element(int i){
		return elements[i];
	}

	void update(int n, T* elements){
		this->n = n;
		this->elements = elements;
	}

	bool operator < (const GeneralSet& that) const{
		if(this->n < that.n){
			return true;
		}
		else if(this->n > that.n){
			return false;
		}
		else if(this->n > 0){
			// Check if one is empty
			if(this->elements == NULL || that.elements == NULL){
				std::cout << "At least, one set is empty" << std::endl;
				return false;
			}

			// Compare each element
			for(int i=0; i<n; i++){
				if(this->elements[i] < that.elements[i]){
					return true;
				}
				else if(this->elements[i] > that.elements[i]){
					return false;
				}
			}

			return false;
		}
		else{
			return false;
		}
	}

	bool operator > (const GeneralSet& that) const{
		if(this->n > that.n){
			return true;
		}
		else if(this->n < that.n){
			return false;
		}
		else if(this->n > 0){
			// Check if one is empty
			if(this->elements == NULL || that.elements == NULL){
				std::cout << "At least, one set is empty" << std::endl;
				return false;
			}

			for(int i=0; i<n; i++){
				if(this->elements[i] > that.elements[i]){
					return true;
				}
				else if(this->elements[i] < that.elements[i]){
					return false;
				}
			}
			return false;
		}
		else{
			return false;
		}
	}

	bool operator == (const GeneralSet& that) const{
		if(this->n == that.n){
			for(int i=0; i<n; i++){
				if(this->elements[i] != that.elements[i]){
					return false;
				}
			}
		}
		else{
			return false;
		}
		return true;
	}

	void print(){
		std::cout << "n = " << n << ": ";
		for(int i=0; i<n; i++){
			std::cout << elements[i] << "  ";
		}
		std::cout << std::endl;
	}

	void sorted(){
		std::sort(elements, elements+n);
	}

	friend std::ostream& operator << (std::ostream& o, const GeneralSet<T>& c){
		o << "n = " << c.n << ": ";
		for(int i=0; i<c.n; i++){
			o << c.elements[i] << "  ";
		}
		o << std::endl;
		return o;
	}
};


class MonIdxSet{
	friend class MonSet;

	IntSet pos;
	int eq_idx;
	int mon_idx;
	bool sys_idx;
	CT coef;

	void init(const int n, const int* elements){
		pos = IntSet(n, elements);
		eq_idx = 0;
		mon_idx = 0;
		sys_idx = 0;
		coef = 0.0;
	}

	void init(const int n, const int* elements, int eq_idx, int mon_idx, bool sys_idx){
		pos = IntSet(n, elements);
		this->eq_idx = eq_idx;
		this->mon_idx = mon_idx;
		this->sys_idx = sys_idx;
		this->coef = 0.0;
	}

	void init(const int n, const int* elements, int eq_idx, int mon_idx, bool sys_idx, const CT& coef){
		pos = IntSet(n, elements);
		this->eq_idx = eq_idx;
		this->mon_idx = mon_idx;
		this->sys_idx = sys_idx;
		this->coef = coef;
	}

	void init(const IntSet& pos, int eq_idx, int mon_idx, bool sys_idx, const CT& coef){
		this->pos = pos;
		this->eq_idx = eq_idx;
		this->mon_idx = mon_idx;
		this->sys_idx = sys_idx;
		this->coef = coef;
	}

public:
	MonIdxSet(){
		eq_idx = 0;
		mon_idx = 0;
		sys_idx = 0;
		coef = 0.0;
	};

	MonIdxSet(const int n, const int* elements){
		init(n, elements);
	}

	MonIdxSet(const int n, const int* elements, int eq_idx, int mon_idx, bool sys_idx){
		init(n, elements, eq_idx, mon_idx, sys_idx);
	}

	MonIdxSet(const int n, const int* elements, int eq_idx, int mon_idx, bool sys_idx, const CT& coef){
		init(n, elements, eq_idx, mon_idx, sys_idx, coef);
	}

	MonIdxSet(const MonIdxSet& original){
		init(original.pos, original.eq_idx, original.mon_idx, original.sys_idx, original.coef);
	}

	MonIdxSet& operator= (const MonIdxSet& original){
		init(original.pos, original.eq_idx, original.mon_idx, original.sys_idx, original.coef);
		return *this;
	}

	~MonIdxSet(){
		//std::cout << "delete MonIdxSet " << this << endl;
	}

	// Basic functions

	IntSet get_pos(){
		return pos;
	}

	int get_eq_idx(){
		return eq_idx;
	}

	int get_mon_idx(){
		return mon_idx;
	}

	bool get_sys_idx(){
		return sys_idx;
	}

	CT get_coef(){
		return coef;
	}

	// comparison operators
	bool operator < (const MonIdxSet& that) const{
		if(this->pos < that.pos){
			return true;
		}
		else if(this->pos == that.pos){
			// Compare equation index
			if(this->eq_idx < that.eq_idx){
				return true;
			}
			else if(this->eq_idx > that.eq_idx){
				return false;
			}

			// Compare system index
			if(this->sys_idx < that.sys_idx){
				return true;
			}
			else if(this->sys_idx > that.sys_idx){
				return false;
			}
			return false;
		}
		else{
			return false;
		}
	}

	bool operator > (const MonIdxSet& that) const{
		if(this->pos > that.pos){
			return true;
		}
		else if(this->pos == that.pos){
			// Compare equation index
			if(this->eq_idx > that.eq_idx){
				return true;
			}
			else if(this->eq_idx < that.eq_idx){
				return false;
			}

			// Compare system index
			if(this->sys_idx > that.sys_idx){
				return true;
			}
			else if(this->sys_idx < that.sys_idx){
				return false;
			}

			return false;
		}
		else{
			return false;
		}
	}

	bool operator == (const MonIdxSet& that) const{
		if(this->pos == that.pos){
			return true;
		}
		else{
			return false;
		}
	}

	void print(){
		pos.print();
		std::cout << "s" << sys_idx
				  << " e" << eq_idx
				  << " m" << mon_idx << std::endl;
		std::cout << coef;
		std::cout << std::endl;
	}

	void sorted(){
		pos.sorted();
	}
};

class EqIdxCoef{
	int eq_idx;
	CT coef[2];

	void init(const int eq_idx, const CT& coef0, const CT& coef1){
		this->eq_idx = eq_idx;
		coef[0] = coef0;
		coef[1] = coef1;
	}
public:
	EqIdxCoef(){
		eq_idx = 0;
	}
	EqIdxCoef(const int eq_idx, const CT& coef0, const CT& coef1){
		init(eq_idx, coef0, coef1);
	}

	EqIdxCoef& operator= (const EqIdxCoef& original){
		init(original.eq_idx, original.coef[0], original.coef[1]);
		return *this;
	}

	EqIdxCoef(const int eq_idx, const CT& coef, const bool sys){
		this->eq_idx = eq_idx;
		this->coef[sys] = coef;
		this->coef[(!sys)] = 0.0;
	}

	void print(){
		std::cout << "Eq " << eq_idx << std::endl;
		std::cout << coef[0];
		std::cout << coef[1];
	}

	void write_coef(CT*& tmp_coef){
		*tmp_coef++ = coef[0];
		*tmp_coef++ = coef[1];
	}

	int get_eq_idx(){
		return eq_idx;
	}

	friend std::ostream& operator << (std::ostream& o, const EqIdxCoef& c){
	   return o << "Eq " << c.eq_idx << std::endl << c.coef[0] << c.coef[1];
	}
};

class MonSet{
	IntSet pos;
	EqSet eq_idx;

	void init(const int n, const int* elements){
		pos = IntSet(n, elements);
	}

	void init(const MonSet& original){
		pos = original.pos;
		eq_idx = original.eq_idx;
	}

public:
	MonSet(){
	};

	void copy_pos(const MonIdxSet& original){
		this->pos = original.pos;
	}

	void update_eq_idx(int n, EqIdxCoef* eq_idx_elements){
		eq_idx.update(n,eq_idx_elements);
	}

	MonSet(const int n, const int* elements){
		init(n, elements);
	}

	MonSet(const MonSet& original){
		init(original);
	}

	MonSet& operator= (const MonSet& original){
		init(original);
		return *this;
	}

	~MonSet(){
	}

	bool operator < (const MonSet& that) const{
		if(this->pos < that.pos){
			return true;
		}
		else{
			return false;
		}
	}

	bool operator > (const MonSet& that) const{
		if(this->pos > that.pos){
			return true;
		}
		else{
			return false;
		}
	}

	bool operator == (const MonSet& that) const{
		if(this->pos == that.pos){
			return true;
		}
		else{
			return false;
		}
	}

	void sorted(){
		pos.sorted();
	}

	int get_n(){
		return pos.get_n();
	}

	int get_n_mon(){
		return eq_idx.get_n();
	}

	int get_pos(int pos_idx){
		return pos.get_element(pos_idx);
	}

	int get_eq_idx(int idx){
		return eq_idx.elements[idx].get_eq_idx();
	}

	void write_coef(CT*& tmp_coef){
		for(int i=0; i<eq_idx.n; i++){
			eq_idx.elements[i].write_coef(tmp_coef);
		}
	}

	void write_coef(CT*& tmp_coef, int mon_idx){
		eq_idx.elements[mon_idx].write_coef(tmp_coef);
	}

	void write_pos(int*& tmp_pos){
		for(int i=0; i<pos.n; i++){
			(*tmp_pos++) = pos.elements[i];
		}
	}

	void write_pos(unsigned short*& tmp_pos){
		for(int i=0; i<pos.n; i++){
			(*tmp_pos++) = pos.elements[i];
		}
	}

	friend std::ostream& operator << (std::ostream& o, const MonSet& c){
	   return o << c.pos << c.eq_idx << std::endl;
	}
};

#endif /* VARSET_H_ */
