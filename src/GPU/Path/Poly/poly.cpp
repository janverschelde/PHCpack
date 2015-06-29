#include "poly.h"

void PolyMon::read(const string& mon_string, VarDict& pos_dict, int* max_deg)
{
    int end = mon_string.length();
    int loc = 0;
    CT coef = get_coef_complex(mon_string, loc);
    read(mon_string, pos_dict, loc, end, coef, max_deg);
}

void PolyMon::read(const string& mon_string, VarDict& pos_dict, CT coef, int* max_deg)
{
    int end = mon_string.length();
    read(mon_string, pos_dict, 0, end, coef, max_deg);
}


inline int var_to_pos(string var){
    string var1 = var.substr(1,var.length()-1);
	//return atoi(var1.c_str());
	//return atoi(var1.c_str())-1;
	return (atoi(var1.c_str())+6)%8;
}

bool compare_sol(PolySol* sol1, PolySol* sol2){
	if(*sol1 < *sol2){
		return true;
	}
	else{
		return false;
	}
}

void PolyMon::read(const string& eq_string, VarDict& pos_dict, int start, int end, CT coef, int* max_deg)
{
    n_var = 1;

    this->coef = coef;

    for(int i = start; i<end; ++i){
        if(eq_string[i] == '*'){
            if(eq_string[i+1] == '*'){
                i++;
            }
            else{
                n_var++;
            }
        }
    }

    pos = new int[n_var];
    exp = new int[n_var];

    string var;
    int cur_type = 0, pos_ind = 0;

    for(int i=start, next_type=1; i<end; ++i){

        char c = eq_string[i];

        int new_var = 0;

        // handle symbols
        if(c == '*'){
            if(eq_string[i+1] == '*'){
                next_type = 2;
                i++;
            }
            else{
                next_type = 1;
            }
            new_var++;
        }
        else if(c=='^'){
            next_type = 2;
            new_var++;
        }
        else if(c == ' ' || c ==';'){
            new_var++;
        }

        // handle vars
        if(new_var == 1){
            if(cur_type == 1){
            	//pos[pos_ind] = var_to_pos(var);
                //string var1 = var.substr(1,var.length()-1);
                //cout << var << " " << var1 << " " << atoi(var1.c_str()) << endl;
                //pos[pos_ind] = atoi(var1.c_str())-1;
                pos[pos_ind] = pos_dict.get(var);
                exp[pos_ind] = 1;
                pos_ind++;
                cur_type = 0;
            }
            else if(cur_type == 2){
                int tmp_exp =  atoi(var.c_str());
                exp[pos_ind-1] = tmp_exp;
                int tmp_pos = pos[pos_ind-1];
                if(tmp_exp > max_deg[tmp_pos]){
                    max_deg[tmp_pos] = tmp_exp;
                }
                cur_type = 0;
            }
        }
        else{
            if(cur_type == 0){
                var = c;
                cur_type = next_type;            
            }
            else{
                var += c;
            }
        }        
    }

    if(cur_type == 1){
        string var1 = var.substr(1,var.length()-1);
        //cout << var << " " << var1 << " " << atoi(var1.c_str()) << endl;
    	//pos[pos_ind] = var_to_pos(var);
        //pos[pos_ind] = atoi(var1.c_str())-1;
        pos[pos_ind] = pos_dict.get(var);
        exp[pos_ind] = 1;
    }
    else if(cur_type == 2){
        exp[pos_ind-1] = atoi(var.c_str()); 
    }
}

/*
 //speel without coefficient
CT PolyMon::speel(const CT* x_val, CT* deri){
    deri[0].init(1.0,0.0);
    deri[1] = x_val[pos[0]];

    for(int i=1; i<n_var-1; i++){
        deri[i+1] = deri[i]*x_val[pos[i]];
    }

    CT tmp = x_val[pos[n_var-1]];
    for(int i=n_var-2; i>-1; i--){
        deri[i] *= tmp;
        tmp *= x_val[pos[i]];
    }
    return tmp;
}*/

CT PolyMon::speel(const CT* x_val, CT* deri){
    deri[1] = x_val[pos[0]];

    for(int i=1; i<n_var-1; i++){
        deri[i+1] = deri[i]*x_val[pos[i]];
    }

    CT tmp = coef;
    for(int i=n_var-1; i>0; i--){
        deri[i] *= tmp;
        tmp *= x_val[pos[i]];
    }

    deri[0] = tmp;

    return tmp*x_val[pos[0]];
}

CT pow(const CT x, int exp){
    CT tmp = x;
    for(int i=1; i<exp; i++){
        tmp *= x;
    }
    return tmp;
}

CT PolyMon::eval(const CT* x_val){
    CT val = coef;
    for(int i=0; i<n_var; i++){
        val *= pow(x_val[pos[i]],exp[i]);
    }
    return val;
}

CT PolyMon::eval_base(const CT* x_val){
    CT val = coef;
    for(int i=0; i<n_var; i++){
        if(exp[i] > 1){
            val *= pow(x_val[pos[i]],exp[i]-1);
        }
    }
    return val;
}

CT PolyMon::eval(const CT* x_val, CT* deri){
    CT val = speel(x_val, deri);
    /*CT base = eval_base(x_val);
    val *= base;

    for(int i=0; i<n_var; i++){
        deri[i] *= base*exp[i];
    }*/
    return val;
}


CT PolyEq::eval(const CT* x_val){
    CT val = constant;
    for(int i=0; i<n_mon; i++){
        val += mon[i]->eval(x_val);
    }
    return val;
}

CT PolyEq::eval(const CT* x_val, CT* deri){
    for(int i=0; i<dim; i++){
        deri[i].init(0.0,0.0);
    }
    CT val = constant;
    //std::cout << constant << std::endl;
    
    CT* mon_deri = new CT[dim];

    for(int i=0; i<n_mon; i++){
        PolyMon* m = mon[i];
        val += m->eval(x_val, mon_deri);

        for(int j=0; j<m->n_var; j++){
            deri[m->pos[j]] += mon_deri[j];
        }
    }

    delete [] mon_deri;

    return val;
}

CT* PolySys::eval(const CT* x_val){
    CT* val = new CT[n_eq];
    for(int i=0; i<n_eq; i++){
        val[i] = eq[i]->eval(x_val);
    }
    return val;
}

CT* PolySys::eval(const CT* x_val, CT** deri_val){
    CT* val = new CT[n_eq];
    for(int i=0; i<n_eq; i++){
        val[i] = eq[i]->eval(x_val, deri_val[i]);
    }
    return val;
}

void PolySys::eval(const CT* x_val, CT* f_val, CT** deri_val){
    for(int i=0; i<n_eq; i++){
        //cout << "eq " << i << endl;
        f_val[i] = eq[i]->eval(x_val, deri_val[i]);
    }
}

void PolyMon::print(const string* pos_var)
{
    // print coef
	print_coef_complex(coef);
	//cout << endl;
    cout << pos_var[pos[0]];
    if(exp[0]!= 1){
        cout << '^' << exp[0];
    }

    for(int i =1; i< n_var; i++){
        cout << " * " << pos_var[pos[i]];
        if(exp[i]!= 1){
            cout << '^' << exp[i];
        }
    }
}

void PolyEq::read(const string& eq_string, VarDict& pos_dict, int* max_deg){
    int l = eq_string.length();
    read(eq_string, pos_dict, 0, l, max_deg);
}

void PolyEq::read(const string& eq_string, VarDict& pos_dict, int start, int end, int* max_deg){
    n_mon = 0;

    // Get the starting position of first monomial
    int mon_first = start;
    for(int i=start; i< end; i++){
        char c = eq_string[i];
        if(c != ' '){
            mon_first = i;
            break;
        }
    }

    // Count the number of monomials
    // Generate monomials link list, include coef, start and end pos in string

    int mon_start = mon_first;
    int i=mon_first;
    int parentheses_open = 0;
    while(i<end){
        char c = eq_string[i];
        if(c == '('){
        	parentheses_open = 1;

        }
        else if (c == ')'){
        	parentheses_open = 0;
        }
        // If a new monomial appears, record the starting and ending position of last one.
        if(((c== '+' || c == '-') && eq_string[i-1]!='e' && eq_string[i-1]!='E' \
           && eq_string[i-1]!='(') && (i != mon_first) && parentheses_open == 0){
            int tmp_start = mon_start;
            CT tmp_coef = get_coef_complex(eq_string, mon_start);
        	parentheses_open = 0;

            if(mon_start < i){
            	/*for(int k=mon_start; k<end; k++){
            		std::cout << eq_string[k];
            	}
            	std::cout << std::endl;*/
                n_mon++;
                //std::cout << "n_mon = " << n_mon << std::endl;
                PolyMon* mm = new PolyMon;
                mm->read(eq_string, pos_dict, mon_start, i, tmp_coef, max_deg);
                mon.push_back(mm);
                mon_start = i;
            }
            else{
            	//std::cout << tmp_coef;
                constant += tmp_coef;
            }
        }
        i++;
    }

    if (mon_start < end){
        CT tmp_coef = get_coef_complex(eq_string, mon_start);
        if(mon_start < end){
            n_mon++;
            PolyMon* mm = new PolyMon;
            mm->read(eq_string, pos_dict, mon_start, end, tmp_coef, max_deg);
            mon.push_back(mm);
        }
        else{
            constant += tmp_coef;
        }
    }

    dim = pos_dict.n_job;
}

void PolyEq::print(const string* pos_var){
	//std::cout << "n_mon = " << n_mon << std::endl;
    for(int i=0; i<n_mon; i++){
        mon[i]->print(pos_var);
    }

    // print plus and minus
    //if(constant.real > 0){
    //    cout << " + ";
    //}

    // print plus and minus
    print_number_complex(constant);
    //cout << constant;
    cout << endl;
}

void PolySys::print(){
    cout << "dim = " << dim << ", n_eq = " << n_eq << endl;
    for(int i=0; i<n_eq; i++){
        cout << "f"<< i << "=" << endl;
        eq[i]->print(pos_var);
    }
    /*cout << "max deg" << endl;
    for(int i=0; i<dim; i++){
        cout << pos_var[i] << " " << max_deg[i] << endl;
    }*/
}

void PolySys::read(const string& sys_string, VarDict& pos_dict)
{
    int l = sys_string.length();
    n_eq = 1;

    int eq_first = 0;

    LinkList<int>* eq_list = new LinkList<int>();

    eq_list->append(eq_first);

    int ii = eq_first;
    while(ii<l){
        if(sys_string[ii] == ';'){
            ii++;
            while(ii<l){
                if(sys_string[ii] != ' '){
                    n_eq++;
                    eq_list->append(ii);
                    break;
                }
                ii++;
            }
        }
        ii++;
    }
    eq_list->append(l);

    //eq = new PolyEq[n_eq];

    LinkNode<int>* tmp = eq_list->header();

    int eq_start = tmp->data;
    for(int i=0; i<n_eq; i++){
        PolyEq* new_eq = new PolyEq();
        tmp = tmp->next;
        int eq_end = tmp->data;
        //cout << n_eq << " "<< eq_start <<" " << eq_end<< endl;
        //string new_eq = sys_string.substr(eq_start, eq_end-eq_start);
        //cout << new_eq << endl;
        new_eq->read(sys_string, pos_dict, eq_start, eq_end, max_deg);
        eq.push_back(new_eq);
        eq_start = eq_end;
    }

    eq_list->destroy();
    delete eq_list;

    dim = pos_dict.n_job;
    pos_var = pos_dict.reverse();
    //dim = pos_dict.n_job;
}


void PolySys::read(const string* sys_string, int n_eq, VarDict& pos_dict){
    eq.reserve(n_eq);
    this->n_eq = n_eq;
    dim = n_eq;
    max_deg = new int[dim];
    for(int i=0; i<dim; i++){
        max_deg[i] = 0;
    }

    // ***Here is confusion, I should use either vector or array, but not to mix them


    eq_space = new PolyEq[n_eq];
    PolyEq* tmp_eq_space = eq_space;
    for(int i=0; i< n_eq; i++){
        //cout << "eq" << i << endl;
    	tmp_eq_space->read(sys_string[i], pos_dict, max_deg);
    	tmp_eq_space->dim = dim;
        eq.push_back(tmp_eq_space);
        tmp_eq_space++;
    }
}


void PolySys::read_file(const string& file_name){
    VarDict pos_dict;
    ifstream myfile (file_name.c_str());
    if (myfile.is_open()){
        read_file(myfile, pos_dict);
        myfile.close();
    }
    else{
        cout << "There is no such file." << endl;
    }
    //dim = pos_dict.n_job;
    pos_var = pos_dict.reverse();
    std::cout << "dim = " << dim << std::endl;
    for(int var_idx=0; var_idx<dim; var_idx++){
    	std::cout << var_idx << " " << pos_var[var_idx] << std::endl;
    }
}

int not_empty_line(const string& line){
    int not_empty = -1;
    for(unsigned int j=0; j < line.length(); j++){
        char c = line[j];
        if(c!=' '&&c!=';'&&c!='\n'){
            not_empty = j;
        }
    }
    return not_empty;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    std::string empy_str = "";
    while (std::getline(ss, item, delim)) {
    	if(item != empy_str){
    		elems.push_back(item);
    	}
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void PolySys::read_file(ifstream& myfile, VarDict& pos_dict){
    string line;
    getline(myfile,line);

    while(true){
        if(not_empty_line(line)!=-1){
            break;
        }
        getline(myfile,line);
    }

    std::vector<std::string> line_parts = split(line, ' ');
    if(line_parts.size() == 2){
        dim = atoi(line_parts[1].c_str());
    	n_eq = atoi(line_parts[0].c_str());
    }
    else{
        dim = atoi(line_parts[0].c_str());
    	n_eq = dim;
    }

    max_deg = new int[dim];
    for(int i=0; i<dim; i++){
        max_deg[i] = 0;
    }
    
    eq.reserve(n_eq);
    eq_space = new PolyEq[n_eq];
    PolyEq* tmp_eq_space = eq_space;
    /*for(int i=0; i< n_eq; i++){
        eq.push_back(tmp_eq_space);
    	tmp_eq_space->dim = dim;
        tmp_eq_space++;
    }*/

    string tmp_line;
    line = "";
    int n_eq_file = 0;
    while( getline(myfile,tmp_line) ){
        int not_empty = not_empty_line(tmp_line);

        if(not_empty != -1){
            int finish = 0;
            for(unsigned int j=not_empty; j<tmp_line.length(); j++){
                char c =tmp_line[j];
                if(c ==';'){
                    finish = 1;
                }
            }
            line += tmp_line;
            if(finish){
            	//std::cout << line << std::endl << std::endl;
            	tmp_eq_space->read(line, pos_dict, max_deg);
            	tmp_eq_space->dim = dim;
                eq.push_back(tmp_eq_space);
                tmp_eq_space++;
                n_eq_file++;
                //std::cout << "n_eq_file = " << n_eq_file << std::endl;
                if(n_eq_file >= n_eq){
                    cout << "Error: too many lines" << endl;
                    break;
                }
                line = "";
            }
        }
    }
    pos_var = pos_dict.reverse();
}

/*void PolyMon::print_level(){
    for(int i=0; i<n_level; i++){
        cout << "    level " << i  << ": n_job = "<< n_pos_level[i]<< endl;
        cout << "        pos: ";
        for(int j=0; j<n_pos_level[i]; j++){
            cout << pos_level[i][j] << ", ";
        }
        cout << " last l = " << pos_level_last[i] << endl;
    }
}

void PolyEq::print_level(){
    for(int i=0; i<n_mon; i++){
        cout << "  Monomial " << i << endl;
        mon[i] -> print_level();
    }
}

void PolySys::print_level(){
    for(int i=0; i<n_eq; i++){
        cout << "Equation " << i << endl; 
        eq[i] -> print_level();
    }
}

int PolyEq::workspace_size_block(int start_level, int factor_size){
    int level = 0;
    for(int i=0; i<n_mon; i++){
        if(level < mon[i]->n_level){
            level = mon[i]->n_level;
        }
    }
    if(start_level >= level){
        return 0;
    }

    int* level_start = new int[level+1];
    int* base_size = new int[level];

    for(int i=0; i<start_level; i++){
        level_start[i] = 0;
    }

    int tmp_size = 1;
    for(int i=1; i<start_level; i++){
        tmp_size *= factor_size;
        base_size[i-1] = tmp_size + 1;
    }
    for(int i=start_level; i<level+1; i++){
        //cout << "i = " << i << endl;
        tmp_size *= factor_size;
        base_size[i-1] = tmp_size + 1;
        level_start[i] = 0;
    }

    for(int j=0; j<n_mon; j++){
        for(int i=start_level; i<mon[j]->n_level; i++){
            level_start[i] += mon[j]->n_pos_level[i]*base_size[i-1];
        }
    }

    for(int i=start_level; i<level+1; i++){
        level_start[i] += level_start[i-1];
    }
    int workspace_size = level_start[level];

    delete[] base_size;
    delete[] level_start;

    return workspace_size;
}

int PolySys::workspace_size_block(int start_level, int factor_size){
    // max workspace for single equation
    int m_size = 0;
    long eq_workspace_max = 0;
    for(int j=0; j<n_eq; j++){
        int new_m_size = eq[j]->memory_size(factor_size);
        if(m_size < new_m_size){
            m_size = new_m_size;
        }
        int eq_workspace_size = eq[j]->workspace_size_block(start_level, factor_size);
        if(eq_workspace_max < eq_workspace_size){
            eq_workspace_max = eq_workspace_size;
        }
    }
    //cout << "eq_workspace_max = " << eq_workspace_max << endl;
    return eq_workspace_max;
}

int PolyMon::memory_size(int factor_size){
    return n_var + (n_var + factor_size - 1) / factor_size;
}

int PolyEq::memory_size(int factor_size){
    int m_size = 0;
    for(int i=0; i< n_mon; i++){
        m_size += mon[i]->memory_size(factor_size);
    }
    return m_size;
}

int PolySys::memory_size(int factor_size){
    int m_size = 0;
    for(int i=0; i < n_eq; i++){
        m_size += eq[i]->memory_size(factor_size);
    }
    return m_size;
}

void PolySys::job_number(){
    for(int i=0; i<level; i++){
        cout << i << " " << job_number_level[i] << endl;
    }
}

int PolyMon::job_number_block(int start_level){
    int n_job_block_est = 0;
    if(start_level <= n_level){
        for(int i=start_level; i<n_level; i++){
            //cout << i << " " <<  n_pos_level[i] << endl;
            n_job_block_est += n_pos_level[i];
        }
    }
    return n_job_block_est;
}

int PolyEq::job_number_block(int start_level){
    int n_job_block_est = 0;
    //cout << "n_mon = " << n_mon << endl;
    for(int i=0; i<n_mon; i++){
        n_job_block_est += mon[i] -> job_number_block(start_level);
    }
    return n_job_block_est;
}

int PolySys::job_number_block(int start_level){
    cout << endl;
    int n_job_block_est = 0;
    for(int i=0; i<n_eq; i++){
        //cout << "eq    " << i << endl;
        int tmp = eq[i] -> job_number_block(start_level);
        n_job_block_est += tmp;
        //cout << "eq total  " << tmp << endl;
    }
    return n_job_block_est;
}*/

void PolySys::gpu_mon(int& dim, int& level, int& workspace_size, int*& workspace_level,
                      int*& n_mon_level, int& pos_size, unsigned short*& pos, int*& pos_level,
                      int& sum_level, int*& n_sum_level,
                      int& total_n_sum, int& sum_array_size, int*& sum_start, int*& sum_array){
    dim = this->dim;                  
    
    int max_n_var = 0;
    int total_n_mon = 0;

    // Count the largest n_var and the total n_mon
    for(int i=0; i<n_eq; i++){
        PolyEq* eq_tmp = eq[i];
        for(int j=0; j<eq_tmp->n_mon; j++){
            int new_n_var = eq_tmp->mon[j]->n_var;
            if( new_n_var > max_n_var){
                max_n_var = new_n_var;
            }
            total_n_mon++;
        }
    }
    //cout << "max_n_var = " << max_n_var << endl;
    
    // Put monomials index into vector ordered by n_var
    vector<int_idx>* n_var_list = new vector<int_idx>[max_n_var+1];
    for(int i=0; i<n_eq; i++){
        PolyEq* eq_tmp = eq[i];
        for(int j=0; j<eq_tmp->n_mon; j++){
            int new_n_var = eq_tmp->mon[j]->n_var;
            n_var_list[new_n_var].push_back(int_idx(i,j));
            //cout << new_n_var << " " << i << " " << j << endl;
        }
    }
    
    /*for(int i=0; i< max_n_var+1; i++){
        cout << i << " size = " << n_var_list[i].size() << endl;
        for (vector<int_idx>::iterator it = n_var_list[i].begin(); it!=n_var_list[i].end(); ++it) {
            cout << (*it).eq_idx << " " << (*it).mon_idx << endl;
        }
    }*/
    
    // compute level
    level = 0;
    for(int i=1; i<max_n_var; i<<=1){
        level++;
        //cout << "i = " << i << " level = " << level << endl;
    }
    cout << " level = " << level << endl;
    
    
    // Count monomials on each level and count n_vars
    n_mon_level = new int[level + 1];
    
    for(int i=0; i<=level; i++){
        n_mon_level[i] = 0;
        int m_size = 1<<i;
        for(int j= m_size/2+1; j<= min(max_n_var, m_size); j++){
            int size_tmp = n_var_list[j].size();
            n_mon_level[i] += size_tmp;
        }
    }
    
    // Count pos on each level and total pos
    pos_level = new int[level + 1];
    pos_size = 0;
    workspace_level = new int[level + 1];
    workspace_size = 0;
    pos_level[0] = 0;
    workspace_level[0] = 0;
    int pos_level_size = 0;
    for(int i=0; i<=level; i++){
        pos_level[i] = pos_level_size;
        workspace_level[i] = pos_level_size;
        pos_level_size = n_mon_level[i] * ((1<<i)+1);
        //cout << i << " pos_level_size = " << pos_level_size 
        //          << " level_start = " << pos_level[i] << endl;
        pos_size += pos_level_size;
        workspace_size += pos_level_size;
    }
    //cout << "pos_size = " << pos_size << " workspace_size = " << workspace_size << endl;
    
    // Put data into pos
    pos = new unsigned short[pos_size];    
    unsigned short* pos_tmp = pos;
    for(int i=0; i<=level; i++){
        int m_size = 1<<i;
        //cout << "m_size = " << m_size << endl;
        for(int j= m_size/2+1; j<= min(max_n_var, m_size); j++){
            //cout << "size = " << j << endl;
            for (vector<int_idx>::iterator it = n_var_list[j].begin(); it!=n_var_list[j].end(); ++it) {
                int eq_idx = (*it).eq_idx;
                int mon_idx = (*it).mon_idx;
                int* mon_pos = eq[eq_idx] -> mon[mon_idx] -> pos;
                *pos_tmp++ = j;
                for(int k=0; k<j; k++){
                    pos_tmp[k] = mon_pos[k];
                    //cout << pos_tmp[k] << " ";
                }
                //cout << endl;
                pos_tmp += m_size;
                //cout << m_size << endl;
            }
        }
    }
    
    // print pos
    pos_tmp = pos;
    for(int i=1; i<=level; i++){
        int m_size = 1<<i;
        //cout << "level = " << i << " m_size = " << m_size << " level_start = " << pos_level[i] << endl;
        pos_tmp += pos_level[i];
        for(int j=0; j<n_mon_level[i]; j++){
            //cout << m_size << endl;
            unsigned short* pos_start = pos_tmp + (m_size+1)*j;
            int n_var_tmp = *pos_start++;
            /*cout << "size = " << n_var_tmp << "   ";
            for(int k=0; k<n_var_tmp; k++){
                cout << pos_start[k] << " ";
            }
            cout << endl;*/
        }
    }
    
    //Get location of each monomials in pos
    //cout << "n_mon = " << n_mon << endl;
    int* mon_pos_data = new int[total_n_mon];
    int** mon_pos = new int*[n_eq];
    mon_pos[0] = mon_pos_data;
    for(int i=1; i<n_eq; i++){
        mon_pos[i] = mon_pos[i-1] + eq[i-1] -> n_mon;
    }
    
    int mon_loc = 0;
    for(int i=0; i<=level; i++){
        int m_size = 1<<i;
        //cout << "m_size = " << m_size << endl;
        for(int j= m_size/2+1; j<= min(max_n_var, m_size); j++){
            //cout << "size = " << j << endl;
            for (vector<int_idx>::iterator it = n_var_list[j].begin(); it!=n_var_list[j].end(); ++it) {
                int eq_idx = (*it).eq_idx;
                int mon_idx = (*it).mon_idx;
                mon_pos[eq_idx][mon_idx] = mon_loc;
                mon_loc += m_size + 1;
            }
        }
    }
    
    /*for(int i=0; i<n_eq; i++){
        cout << "eq " << i << endl;
        for(int j=0; j<eq[i] -> n_mon; j++){
            cout << mon_pos[i][j] << " "; 
        }
        cout << endl;
    }*/
    
    int* sum_size_data = new int[n_eq*(dim+1)];
    for(int i=0; i<n_eq*(dim+1); i++){
        sum_size_data[i] = 0;
    }
    int** sum_size = new int*[n_eq];
    sum_size[0] = sum_size_data;
    for(int i=1; i<n_eq; i++){
        sum_size[i] = sum_size[i-1] + dim + 1;
    }
    
    //cout << "Equation size:";
    /*for(int i=0; i<n_eq; i++){
        //cout << sum_size[n_eq][i] << " ";
    }
    //cout << endl;*/
    
    // Equation Sum Size and Derivative Sum Size
    for(int eq_idx=0; eq_idx<n_eq; eq_idx++){
        sum_size[eq_idx][dim] = eq[eq_idx] -> n_mon; // +1 for constant
        //cout << "eq " << eq_idx << " eq_size " << eq[eq_idx] -> n_mon << endl;
        for(int mon_idx=0; mon_idx<eq[eq_idx] -> n_mon; mon_idx++){
            int mon_start = mon_pos[eq_idx][mon_idx];
            int mon_size = pos[mon_start];
            //cout << "eq " << eq_idx << " mon_idx " << mon_idx << " mon_size " << mon_size << endl;
            for(int i=1; i<mon_size+1; i++){
                (sum_size[eq_idx][pos[mon_start+i]])++;
            }
        }
    }
    
    // Largest number in Sum
    int max_sum_size = 0;
    total_n_sum = 0;
    int total_n_sum_pos = 0;
    for(int eq_idx=0; eq_idx<n_eq; eq_idx++){
        //cout << "eq " << eq_idx << "   ";
        for(int x_idx=0; x_idx < dim+1; x_idx++){
            int tmp_sum_size = sum_size[eq_idx][x_idx];
            if(tmp_sum_size>max_sum_size){
                max_sum_size = tmp_sum_size;
            }
            if(tmp_sum_size > 0){
                total_n_sum++;
                total_n_sum_pos += tmp_sum_size;
            }
            //cout << sum_size[eq_idx][x_idx] << " ";
        }
        //cout << endl;
    }
    cout << "max_sum_size    = " << max_sum_size << endl;
    cout << "total_n_sum     = " << total_n_sum  << endl;
    cout << "total_n_sum_pos = " << total_n_sum_pos  << endl;
    
    // Put index into Sum Vector
    vector<int_idx>* sum_size_list = new vector<int_idx>[max_sum_size+1];
    for(int eq_idx=0; eq_idx<n_eq; eq_idx++){
        for(int x_idx=0; x_idx<dim+1; x_idx++){
            sum_size_list[sum_size[eq_idx][x_idx]].push_back(int_idx(eq_idx, x_idx));
        }
    }

    // Generate n_sum_level
    sum_level = 0;
    for(int i=1; i<max_sum_size; i<<=1){
        sum_level ++;
        //cout << "i = " << i << " level = " << level << endl;
    }
    cout << "sum_level = " << sum_level << endl;
    
    
    // Count n_sum on each level
    n_sum_level = new int[sum_level + 1];
    
    for(int i=0; i<=sum_level; i++){
        n_sum_level[i] = 0;
        int m_size = 1<<i;
        for(int j= m_size/2+1; j<= min(max_sum_size, m_size); j++){
            int size_tmp = sum_size_list[j].size();
            n_sum_level[i] += size_tmp;
        }
    }

    /*for(int i=0; i<=sum_level; i++){
        cout << i << " " << n_sum_level[i] << endl;
    }*/

    // Generate Sum Start
    int* sum_start_loc_data = new int[n_eq*(dim+1)];
    for(int i=0; i<n_eq*(dim+1); i++){
        sum_start_loc_data[i] = 0;
    }
    
    int** sum_start_loc = new int*[n_eq+1];
    sum_start_loc[0] = sum_start_loc_data;
    for(int i=1; i<n_eq; i++){
        sum_start_loc[i] = sum_start_loc[i-1] + dim + 1;
    }
    
    sum_array_size = total_n_sum*2+total_n_sum_pos;
    sum_array = new int[sum_array_size]; // need to return
    for(int i=0; i<total_n_sum*2+total_n_sum_pos; i++){
        sum_array[i] = 0;
    }

    sum_start = new int[total_n_sum+1]; // need to be return
    int* sum_start_tmp = sum_start;
    *sum_start_tmp = 0;
    
    for(int i=0; i<=max_sum_size; i++){
        //cout << "size = " << i+2  << " n = " << sum_size_list[i].size() << endl;
        for (vector<int_idx>::iterator it = sum_size_list[i].begin(); it!=sum_size_list[i].end(); ++it) {
            //cout << "   " << (*it).eq_idx << " " << (*it).mon_idx << endl;
            int eq_idx = (*it).eq_idx;
            int x_idx = (*it).mon_idx;
            sum_start_loc[eq_idx][x_idx] = *sum_start_tmp;
            sum_array[*sum_start_tmp+1] = eq_idx + x_idx*dim;
            *(sum_start_tmp+1) = (*sum_start_tmp) + i+2;
            sum_start_tmp++;
        }
    }
    
    /*for(int i=0; i<total_n_sum; i++){
        cout << i << " " << sum_start[i] << endl;
    }*/
    
    for(int eq_idx=0; eq_idx<n_eq; eq_idx++){
        //cout << "eq " << eq_idx << endl;
        for(int mon_idx=0; mon_idx < eq[eq_idx]->n_mon; mon_idx++){
            int mon_start = mon_pos[eq_idx][mon_idx];
            int mon_size = pos[mon_start];
            //cout << "eq " << eq_idx << " mon_idx " << mon_idx << " mon_size " << mon_size << endl;
            for(int i=1; i<mon_size+1; i++){
                int x_idx = pos[mon_start+i];
                int deri_start = sum_start_loc[eq_idx][x_idx];
                int current_deri_size = sum_array[deri_start];
                sum_array[deri_start+current_deri_size+2] = mon_start+i;
                sum_array[deri_start]++;
            }
        }
    }
    
    for(int eq_idx=0; eq_idx<n_eq; eq_idx++){
        //cout << "eq[eq_idx]->n_mon = " << eq[eq_idx]->n_mon << endl;
        int deri_start = sum_start_loc[eq_idx][dim];
        //cout << "deri_start = " << deri_start << endl;
        for(int mon_idx=0; mon_idx < eq[eq_idx]->n_mon; mon_idx++){
            int mon_start = mon_pos[eq_idx][mon_idx];
            int current_deri_size = sum_array[deri_start];
            //cout << "current_deri_size = " << current_deri_size << endl;
            sum_array[deri_start+current_deri_size+2] = mon_start;
            sum_array[deri_start]++;
        }
    }
    
    
    /*cout << "sum array" << endl;
    for(int i=0; i<total_n_sum; i++){
        int sum_start_tmp = sum_start[i];
        for(int j=0; j<sum_array[sum_start_tmp]+2; j++){
            cout << sum_array[sum_start_tmp+j] << " ";
        }
        cout << endl;
    }*/
            
    
    delete[] n_var_list;
    delete[] sum_size;
    delete[] sum_size_data;
    delete[] mon_pos_data;
    delete[] mon_pos;
    delete[] sum_size_list;
    delete[] sum_start_loc;
    delete[] sum_start_loc_data;
}


void PolySol::init(CT* sol, int dim, T1 max_residual, T1 max_delta_x, int path_idx, string path_info){
	this->dim = dim;
	this->sol = new CT[dim];
	for(int i=0; i<dim; i++){
		this->sol[i] = sol[i];
	}
	this->path_idx = path_idx;
	idx = 0;
	m = 0;
	t = CT(0.0,0.0);
	err = max_delta_x;
	rco = 0.0;
	res = max_residual;
	info = path_info;
}

void PolySol::init(int dim, T1 t_real, T1 t_imag, T1* sol, \
		T1 max_delta_x, T1 rco, T1 max_residual, \
		int m, int path_idx, string path_info){
	this->dim = dim;
	this->sol = new CT[dim];
	for(int i=0; i<dim; i++){
		this->sol[i].real = sol[2*i];
		this->sol[i].imag = sol[2*i+1];
	}
	this->path_idx = path_idx;
	idx = 0;
	this->m = m;
	t = CT(t_real,t_imag);
	err = max_delta_x;
	this->rco = rco;
	res = max_residual;
	info = path_info;
}

void PolySol::init(ifstream& myfile, int dim){
	this->dim = dim;
	string tmp_line;

    /*getline(myfile,tmp_line);
	std::cout << tmp_line << std::endl;
	std::cout << "dim = " << dim << std::endl;*/

    // Get idx
    getline(myfile,tmp_line, ' ');
    //std::cout << tmp_line << std::endl;
    myfile >> idx;
    //std::cout << "idx = " << idx << std::endl;
    getline(myfile,tmp_line);

    int tmp_line_l = tmp_line.size();
    if(tmp_line_l > 2){
    	unsigned found = tmp_line.find_last_of(" ");
    	tmp_line = tmp_line.substr(found+1);
    	tmp_line_l = tmp_line.size();
    	if(tmp_line[tmp_line_l-1]<'a' \
    			or tmp_line[tmp_line_l-1]>'z'){
    		tmp_line = tmp_line.substr(0,tmp_line_l-1);
    	}
    	info = tmp_line;
    }

    // Get t
    getline(myfile,tmp_line, ':');
    //std::cout << tmp_line << std::endl;
    t = get_complex_number(myfile);
    //std::cout << t << std::endl;
    getline(myfile,tmp_line);

    // Get m
    getline(myfile,tmp_line, ':');
    //std::cout << tmp_line << std::endl;
    int m;
    myfile >> m;
    //std::cout << "m = " << m << std::endl;
    getline(myfile,tmp_line);

    getline(myfile,tmp_line);
    //std::cout << tmp_line << std::endl;

    sol = new CT[dim];

    for(int i=0; i<dim; i++){
        getline(myfile,tmp_line, ':');
        sol[i] = get_complex_number(myfile);
    	//std::cout<< sol[i];
        getline(myfile,tmp_line);
    }


    // Get error
    getline(myfile,tmp_line, ':');
    myfile >> err;
    //std::cout << "error = " << err << endl;

    getline(myfile,tmp_line, ':');
    myfile >> rco;
    //std::cout << "rco = " << rco << endl;

    getline(myfile,tmp_line, ':');
    myfile >> res;
    //std::cout << "res = " << res << endl;
    getline(myfile,tmp_line);
}

void PolySol::init(ifstream& myfile, int dim, VarDict& pos_dict){
	this->dim = dim;
	string tmp_line;

    /*getline(myfile,tmp_line);
	std::cout << tmp_line << std::endl;
	std::cout << "dim = " << dim << std::endl;*/

    // Get idx
    getline(myfile,tmp_line, ' ');
    //std::cout << tmp_line << std::endl;
    myfile >> idx;
    //std::cout << "idx = " << idx << std::endl;
    getline(myfile,tmp_line);

    int tmp_line_l = tmp_line.size();
    if(tmp_line_l > 2){
    	unsigned found = tmp_line.find_last_of(" ");
    	tmp_line = tmp_line.substr(found+1);
    	tmp_line_l = tmp_line.size();
    	if(tmp_line[tmp_line_l-1]<'a' \
    			or tmp_line[tmp_line_l-1]>'z'){
    		tmp_line = tmp_line.substr(0,tmp_line_l-1);
    	}
    	info = tmp_line;
    }

    // Get t
    getline(myfile,tmp_line, ':');
    //std::cout << tmp_line << std::endl;
    t = get_complex_number(myfile);
    //std::cout << t << std::endl;
    getline(myfile,tmp_line);

    // Get m
    getline(myfile,tmp_line, ':');
    //std::cout << tmp_line << std::endl;
    int m;
    myfile >> m;
    //std::cout << "m = " << m << std::endl;
    getline(myfile,tmp_line);

    getline(myfile,tmp_line);
    //std::cout << tmp_line << std::endl;

    sol = new CT[dim];

    //pos_dict.print();
    for(int i=0; i<dim; i++){
        getline(myfile,tmp_line, ':');
        //std::cout << tmp_line << tmp_line.substr(1,tmp_line.length()-2);
        int var_idx = pos_dict.get(tmp_line.substr(1,tmp_line.length()-2));
        //std::cout << var_idx << std::endl;
        sol[var_idx] = get_complex_number(myfile);
    	//std::cout<< sol[i];
        getline(myfile,tmp_line);
    }


    // Get error
    getline(myfile,tmp_line, ':');
    myfile >> err;
    //std::cout << "error = " << err << endl;

    getline(myfile,tmp_line, ':');
    myfile >> rco;
    //std::cout << "rco = " << rco << endl;

    getline(myfile,tmp_line, ':');
    myfile >> res;
    //std::cout << "res = " << res << endl;
    getline(myfile,tmp_line);
}

bool PolySol::operator == (const PolySol& that){
	T1 dis = 0;
	for(int i=0; i<dim; i++){
		T1 tmp_dis = abs(sol[i].real - that.sol[i].real);
		if(tmp_dis > dis){
			dis = tmp_dis;
		}
		tmp_dis = abs(sol[i].imag - that.sol[i].imag);
		if(tmp_dis > dis){
			dis = tmp_dis;
		}
	}
	if(dis > 1E-6){
		return 0;
	}
	return 1;
}

void PolySol::print(){
	std::cout << "dim = " << dim << std::endl;
	for(int i=0; i<dim; i++){
		std::cout << i << " " << sol[i];
	}
	std::cout << std::endl;
}

void PolySol::print_info(){
	std::cout << "Solution " << idx << ":" << std::endl;
	std::cout << "t : " << t;
	for(int i=0; i<dim; i++){
		std::cout << i << " : " << sol[i];
	}
	std::cout << "== err : " << err \
		      << " = rco : " << rco \
		      << " = res : " << res \
		      << " ==" << std::endl;
}

void PolySol::print_info(string* pos_var){
	std::cout << "Solution " << idx << ":" << std::endl;
	std::cout << "t : " << t;
	for(int i=0; i<dim; i++){
		std::cout << pos_var[i] << " : " << sol[i];
	}
	std::cout << "== err : " << err \
		      << " = rco : " << rco \
		      << " = res : " << res \
		      << " ==" << std::endl;
}

CT* PolySol::get_sol(){
	CT* sol_tmp = new CT[dim];
	for(int i=0; i<dim; i++){
		sol_tmp[i] = sol[i];
	}
	return sol_tmp;
}

void PolySol::print_short(){
	std::cout << std::scientific \
			  << " " << err  << " " << res << " " << " x[0] = " << sol[0];
}

int to_int(double a){
	int b = a;
	return b;
}

/*int to_int(dd_real a){
	int b = a.x[0];
	return b;
}

int to_int(qd_real a){
	int b = a.x[0];
	return b;
}*/

bool PolySol::operator<(PolySol& that){
	if(dim < that.dim){
		return true;
	}
	if(dim > that.dim){
		return false;
	}

	CT* this_sol = new CT[dim];
	CT* that_sol = new CT[dim];

	for(int i=0; i<dim; i++){
		this_sol[i] = sol[i];
		that_sol[i] = that.sol[i];
	}

	int digits_per_check = 2;

	T1 digits_multiplier = 1;
	T1 err_roundoff = 1E-4;

	for(int i=0; i<digits_per_check; i++){
		digits_multiplier *= 10;
	}

	int digits_checked = 0;

	while(digits_checked < 8){
		for(int i=0; i<dim; i++){
			/*std::cout << i << std::endl;
			std::cout << this_sol[i].real << std::endl;
			std::cout << that_sol[i].real << std::endl;*/
			int this_digit;
			int that_digit;
			int tmp_number;

			// real check
			this_sol[i].real *= digits_multiplier;
			that_sol[i].real *= digits_multiplier;
			/*std::cout << digits_multiplier << std::endl;
			std::cout << this_sol[i].real << std::endl;
			std::cout << that_sol[i].real << std::endl;*/

			this_digit = to_int(this_sol[i].real);
			tmp_number = to_int(this_sol[i].real+err_roundoff);
			if(this_digit != tmp_number){
				this_digit = tmp_number;
			}
			tmp_number = to_int(this_sol[i].real-err_roundoff);
			if(this_digit != tmp_number){
				this_digit = tmp_number;
			}

			that_digit = to_int(that_sol[i].real);
			tmp_number = to_int(that_sol[i].real+err_roundoff);
			if(that_digit != tmp_number){
				that_digit = tmp_number;
			}
			tmp_number = to_int(that_sol[i].real-err_roundoff);
			if(that_digit != tmp_number){
				that_digit = tmp_number;
			}

			//std::cout << this_digit << " " << that_digit<< std::endl;

			if(this_digit < that_digit){
				return true;
			}
			else if(this_digit > that_digit){
				return false;
			}
			this_sol[i].real -= this_digit;
			that_sol[i].real -= that_digit;

			// imag check
			this_sol[i].imag *= digits_multiplier;
			that_sol[i].imag *= digits_multiplier;
			this_digit = to_int(this_sol[i].imag);
			tmp_number = to_int(this_sol[i].imag+err_roundoff);
			if(this_digit != tmp_number){
				this_digit = tmp_number;
			}
			tmp_number = to_int(this_sol[i].imag-err_roundoff);
			if(this_digit != tmp_number){
				this_digit = tmp_number;
			}
			that_digit = to_int(that_sol[i].imag);
			tmp_number = to_int(that_sol[i].imag+err_roundoff);
			if(that_digit != tmp_number){
				that_digit = tmp_number;
			}
			tmp_number = to_int(that_sol[i].imag-err_roundoff);
			if(that_digit != tmp_number){
				that_digit = tmp_number;
			}
			//std::cout << this_digit << " " << that_digit<< std::endl;
			if(this_digit < that_digit){
				return true;
			}
			else if(this_digit > that_digit){
				return false;
			}
			this_sol[i].imag -= this_digit;
			that_sol[i].imag -= that_digit;
		}
		digits_checked += digits_per_check;
	}

	delete[] this_sol;
	delete[] that_sol;

	return false;
}

void PolySolSet::init(ifstream& myfile){
	n_sol = 0;
	dim = 0;

	string prefix_two[2] = {"START SOLUTIONS : ", "THE SOLUTIONS :"};
	read_until_line(myfile, prefix_two, 2);

	myfile >> n_sol;
	myfile >> dim;

	string prefix = "======";
	read_until_line(myfile, prefix);

	sols.reserve(n_sol);

	//std::cout << "n_sol = " << n_sol << std::endl;
	for(int i=0; i<n_sol; i++){
		//std::cout << i << std::endl;
		PolySol* tmp_sol= new PolySol(myfile, dim);
		sols.push_back(tmp_sol);
	}
}

void PolySolSet::init(ifstream& myfile, VarDict& pos_dict){
	n_sol = 0;
	dim = 0;

	string prefix_two[2] = {"START SOLUTIONS : ", "THE SOLUTIONS :"};
	read_until_line(myfile, prefix_two, 2);

	myfile >> n_sol;
	myfile >> dim;

	string prefix = "======";
	read_until_line(myfile, prefix);

	sols.reserve(n_sol);

	//std::cout << "n_sol = " << n_sol << std::endl;
	for(int i=0; i<n_sol; i++){
		//std::cout << i << std::endl;
		PolySol* tmp_sol= new PolySol(myfile, dim, pos_dict);
		sols.push_back(tmp_sol);
	}
}

bool PolySolSet::find_same_sol(PolySol* tmp_sol){
	for(int i=0; i<n_sol; i++){
		if(*tmp_sol == *sols[i]){
			return true;
		}
	}
	return false;
}

int PolySolSet::count_same_sol(PolySol* tmp_sol){
	int n_same_sol = 0;
	for(int i=0; i<n_sol; i++){
		if(*tmp_sol == *sols[i]){
			n_same_sol++;
		}
	}
	return n_same_sol;
}

bool PolySolSet::add_diff_sol(CT* new_sol){
	PolySol* tmp_sol= new PolySol(new_sol, dim);
	if(find_same_sol(tmp_sol)==true){
		return false;
	}
	std::cout << "Add New Solution" << std::endl;
	sols.push_back(tmp_sol);
	n_sol++;
	return true;
}

void PolySolSet::add_sol(PolySol* tmp_sol){
	sols.push_back(tmp_sol);
	n_sol++;
}

void PolySolSet::add_sol(CT* new_sol, T1 max_residual, T1 max_delta_x, int path_idx, string path_info){
	PolySol* tmp_sol= new PolySol(new_sol, dim, max_residual, max_delta_x, path_idx, path_info);
	add_sol(tmp_sol);
}


void PolySolSet::print(){
	std::cout << "dim   = " << dim << std::endl\
			  << "n_sol = " << n_sol << std::endl;

	for(int i=0; i<n_sol; i++){
		sols[i]->print();
	}

}

void PolySolSet::print_info(string* pos_var){
	std::cout << "dim   = " << dim << std::endl\
			  << "n_sol = " << n_sol << std::endl;

	for(int i=0; i<n_sol; i++){
		sols[i]->print_info(pos_var);
	}

}

void PolySolSet::print_short(){
	std::cout << "n_sol = " << n_sol << std::endl;
	for(int i=0; i<n_sol; i++){
		std::cout << i << " " << sols[i]->info << " x[0] = " << sols[i]->sol[0];
	}
}

CT* PolySolSet::get_sol(int idx){
	return sols[idx]->get_sol();
}

void PolySolSet::sort_set(){
	sort(sols.begin(), sols.end(), compare_sol);
}

void PolySolSet::compare(PolySolSet& that){
	vector<PolySol*> that_sols(that.sols);
	vector<PolySol*> this_sols(this->sols);

	sort(this_sols.begin(), this_sols.end(), compare_sol);
	sort(that_sols.begin(), that_sols.end(), compare_sol);

	vector<PolySol*>::iterator this_pointer = this_sols.begin();
	vector<PolySol*>::iterator that_pointer = that_sols.begin();


	vector<PolySol*> this_sols_only;
	vector<PolySol*> that_sols_only;

	int same_sol = 0;
	int this_sol = 0;
	int that_sol = 0;

	while(true){
		if((**this_pointer)==(**that_pointer)){
			this_pointer++;
			that_pointer++;
			same_sol++;
		}
		else{
			if((**this_pointer)<(**that_pointer)){
				this_sols_only.push_back(*this_pointer);
				this_pointer++;
				this_sol++;
			}
			else{
				/*std::cout << this_pointer - this_sols.begin() << std::endl;
				std::cout << that_pointer - that_sols.begin() << std::endl;
				(*this_pointer)->print();
				(*that_pointer)->print();*/
				that_sols_only.push_back(*that_pointer);
				that_pointer++;
				that_sol++;
			}
		}
		if(this_pointer==this_sols.end() || that_pointer==that_sols.end()){
			break;
		}
	}

	if(this_pointer != this_sols.end()){
		this_sol += this_sols.end() - this_pointer;
		while(this_pointer < this_sols.end()){
			this_sols_only.push_back(*this_pointer);
			this_pointer++;
		}
	}

	if(that_pointer != that_sols.end()){
		that_sol += that_sols.end() - that_pointer;
		while(that_pointer < that_sols.end()){
			that_sols_only.push_back(*that_pointer);
			that_pointer++;
		}
	}
	std::cout << "same_sol = " << same_sol << std::endl;

	if(this_sol > 0){
		std::cout << "this_sol = " << this_sol << std::endl;
	}
	if(that_sol > 0){
		std::cout << "that_sol = " << that_sol << std::endl;
	}

	that_pointer = that_sols_only.begin();
	int i=0;
	while(that_pointer < that_sols_only.end()){
		bool find_same_this = find_same_sol(*that_pointer);
		int that_n_same_sol = that.count_same_sol(*that_pointer);
		std::cout << std::setw(4) << i++ << " " \
				  << std::setw(8) << (*that_pointer)->path_idx << " " \
				  << std::setw(4) << find_same_this << " " \
				  << std::setw(4) << that_n_same_sol;
		(*that_pointer)->print_short();
		that_pointer++;
	}
}
