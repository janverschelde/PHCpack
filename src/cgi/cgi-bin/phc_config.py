phc_debug = 0

# PHC_Server host
hostname     = ''  # on same host
hostadd      = ''
creator_name = ''

# SQL Host
sql_host=""
sql_user=""
sql_passwd=""
sql_datebase=""


# PHC_Server
portnumber  = 12091          # same port number
server_address = (hostname, portnumber)
usersfolder = "%s@%s:public_html/users/"%(creator_name, hostname)

file_recover = '../Server/phcserver.recover' # Server Recovery File
file_ext_all = [".pid", ".start", ".cmd",".sta", ".poly",".sol",".phc",".outfile",".maple",".python"]


# Web Interface

header= "Solve Polynomial Systems by PHCpack(beta)"
title = "PHCpack Web Interface(beta)"
solve_sleep_time = 1
small_screen = 400
n_sol_page = 10

# Web Style

style_login     = ['login-small', 'login-large', 'login', 'common']
style_demo     = ['login-small', 'login-large', 'login', 'common']
style_reg = ['login-small', 'login-large', 'login', 'common']
style_activate  = ['login-small', 'login-large', 'login', 'common']
style_contact   = ['login-small', 'login-large', 'login', 'common']
style_forgetpwd = ['login-small', 'login-large', 'login', 'common']
style_about = ['login-small', 'login-large', 'login', 'common']

style_web       = ['style-small', 'style-large', 'style', 'common']


# client message

client_id_len=40
Sid_len=5
buffer = 3072         # Message size


# PHC Solver Client

poly_small_threads= 1 # local threads to solve small polynomial system
poly_small = 25        # small polynomial degree
phc = '../PHC/phc64'
phc_cmd = [phc,'-b']
phc_check = phc + " -g "  # check poly before solving
phc_python = phc + " -x " # print python
phc_maple = phc + " -z "  # print maple
phc_hom = [phc,'-p']    # homotopy by 1 parameter
Solver_Remote_Folder="../TMP" # Remote Solver Temprary Solution Path
sleeptime = 5

# Gmail STMP
phcstmp = ''
phcmail  = ''  
phcmailps= ''

