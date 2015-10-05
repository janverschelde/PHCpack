class PolySol():
    def __init__(self):
        self.dim = 0
        self.ind = 0
        self.m = 0
        self.t_real = 0.0
        self.t_imag = 0.0
        self.sols = []
        self.err = 0.0
        self.rco = 0.0
        self.res = 0.0

    def read_sol(self, f, f_line):
        self.read_sol_ind(f_line)

        self.read_var_value(f.readline())
        
        self.read_sol_mul(f.readline())
        
        f.readline()

        dim = 0
        while True:
            f_line = f.readline()
            if f_line[0] != '=':
                self.read_var_value(f_line)
                dim += 1
            else:
                break;
        self.dim = dim

        self.info_var_value(f_line)

    def read_sol_ind(self, line):
        line = line[8:].split(':')
        self.ind = int(line[0])

    def read_sol_mul(self, line):
        line = line.split(':')
        self.m = int(line[1])
        
    def read_var_value(self,var_string):
        l = len(var_string)
        i = 0
        var_name = ''
        while True:
            c = var_string[i]
            i = i + 1
            if c == ':':
                break
            elif c != ' ':
                var_name += c


        new_element = ' '
        var_real = 0.0
        var_imag = 0.0
        var_end = 1

        while True:
            if var_string[i] != ' ':
                break
            i += 1

        real_start  = i

        while True:
            if var_string[i] == ' ':
                break;
            i += 1

        var_real = float(var_string[real_start:i])
        var_imag = float(var_string[i:])

        if var_name != 't':
            self.sols.append([var_name, var_real, var_imag])
        else:
            self.t_real = var_real
            self.t_imag = var_imag

        
    def info_var_value(self,info_string):
        j = 0
        i = 6
        start_i = 6
        var = [0.0,0.0,0.0]
        for j in range(3):
            var_start = 0
            while True:
                if var_start:
                    if info_string[i] == '=':
                        var[j] = float(info_string[start_i:i])
                        break;
                else:
                    if info_string[i] == ':':
                        var_start = 1
                        start_i = i+1
                i += 1

        self.err = var[0]
        self.rco = var[1]
        self.res = var[2]

    def print_sol(self):
        print 'dim = ', self.dim
        print 'ind = ', self.ind
        print 'm = ', self.m
        print 't_real = ', self.t_real
        print 't_m = ', self.t_imag
        print 'sol = ', self.sols
        print 'err = ', self.err
        print 'rco = ', self.rco
        print 'res = ', self.res

    def print_sol_table(self):
        print """<tr style ="border:1px solid black;"><td style="width:20%%"><b>SOL %d</b></td><td colspan="2", align = "right">  Mult: %d err: %.2e   rco:%.2e   res:%.2e</td></tr>"""%(self.ind,self.m, self.err,self.rco,self.res)
        for i in xrange(self.dim):
            print """<tr class="sol_high"><td>%s</td><td align="right">%.14e</td>\
                         <td align="right">%.14e</td></tr>"""\
                    %(self.sols[i][0],self.sols[i][1],self.sols[i][2])

            print """<tr class="sol_low"><td>%s</td><td align="right">%.5e</td>\
                         <td align="right">%.5e</td></tr>"""\
                    %(self.sols[i][0],self.sols[i][1],self.sols[i][2])


def read_sol_file(filename, sol_page):
    from phc_config import n_sol_page
    with open(filename) as f:
        while True:
            a = f.readline()
            if a[:5] == '=====':
                break
        print """<table style ="border-collapse: collapse;border: 1px solid black;">"""
        sol_start = sol_page*n_sol_page
        sol_end = (sol_page+1)*n_sol_page
        sols = []
        n_sol = 0
        while True:
            f_line = f.readline()
            if f_line =='':
                break
            s = PolySol()
            # this part need to be changed, don't read sols, skip sols before sol_start
            s.read_sol(f,f_line)
            if n_sol < sol_start:
                sols.append('')
            else:
                sols.append(s)
            n_sol += 1
            if(len(sols) > sol_end):
                break

        if (sol_start + n_sol_page > len(sols)):
            sol_end = len(sols)
        else:
            sol_end = sol_start + n_sol_page

        for i in range(sol_start, sol_end):
            sols[i].print_sol_table()

        print '</table>'

def read_sol_number_file(filename):
    with open(filename) as f:
        sol_symbol = '========='
        tmp_line = ''
        while True:
            a = f.readline()
            if ((len(a)>len(sol_symbol)  and a[:len(sol_symbol)] == sol_symbol)) or (a == ''):
                n_sol_line = last_line
                break
            last_line = a
        if n_sol_line != '':
            n_sol_line = n_sol_line.split(' ')
            return int(n_sol_line[0])
        return 0
    
