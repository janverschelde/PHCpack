#include "families.h"

string* x_var(string x, int dim){
	string* var = new string[dim];
	for(int i=0; i< dim; i++){
        ostringstream ss;
		ss << x << i;
		var[i] = ss.str();
	}
	return var;
}

string* string_cyclic(int dim){
    string* poly_sys = new string[dim];

    for(int i=1; i<dim; i++){
        ostringstream ss;
        ss << "x0";
        for(int k=1; k<i; k++){
            ss << "*x" << k;
        }

        for(int j=1; j<dim; j++){
            ss << " + x" << j;
            for(int k=j+1; k<i+j; k++){
                ss << "*x" << k%dim;
            }
        }
        ss << ';';

        poly_sys[i-1] = ss.str();
    }

    ostringstream ss;
    ss << "x0";
    for(int i=1; i<dim; i++){
        ss << "*x" << i;
    }
    ss << "-1;";
    poly_sys[dim-1] = ss.str();

    return poly_sys;
}

void file_cyclic(string f_name, int dim){
    ofstream f;
    f.open(f_name.c_str());
    f << dim << '\n';

    for(int i=1; i<dim; i++){
        f << "x0";
        for(int k=1; k<i; k++){
            f << "*x" << k;
        }

        for(int j=1; j<dim; j++){
            f << " + x" << j ;
            for(int k=j+1; k<i+j; k++){
                f << "*x" << k%dim;
            }
        }
        f << ";\n";
    }

    f << "x0";
    for(int i=1; i<dim; i++){
        f << "*x" << i ;
    }
    f << "-1;" << endl;
    f.close();
}

string* string_two(int dim){
    string* poly_sys = new string[2];

    ostringstream ss;
    ss << "x0";
    for(int i=1; i<dim; i++){
        ss << "*x" << i;
    }
    poly_sys[0] = ss.str();

    ostringstream ss1;
    ss1 << "x2";
    for(int i=3; i<dim; i++){
        ss1 << "*x" << i;
    }
    poly_sys[1] = ss1.str();

    return poly_sys;
}



/*string* string_cyclic(int dim){
    string* poly_sys = new string[dim];

    for(int i=1; i<dim; i++){
        string new_poly = "x0";
        for(int k=1; k<i; k++){
            new_poly += "*x" + to_string(k);
        }

        for(int j=1; j<dim; j++){
            new_poly += " + x" + to_string(j) ;
            for(int k=j+1; k<i+j; k++){
                new_poly += "*x" + to_string(k%dim);
            }
        }
        new_poly += ';';

        poly_sys[i-1] = new_poly;
    }

    string new_poly = "x0";
    for(int i=1; i<dim; i++){
        new_poly += "*x" + to_string(i);
    }
    new_poly += "-1;";
    poly_sys[dim-1] = new_poly;

    return poly_sys;
}*/
