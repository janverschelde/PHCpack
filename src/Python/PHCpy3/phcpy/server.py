"""
The module server.py exports routines to send and receive lists
of strings through sockets.  These strings represent either polynomials
or solutions as data interchanges between a client and a server.
A simple illustration of the use of server is to solve many
polynomial systems over a client/server network.
The interactive main program starts up a multithreaded server.
The server generates a list of random polynomial systems which
are distributed among clients in a static scheme.  The clients
solve the polynomial systems and send solutions to the server.
"""

BUFSIZE = 2048
PORTNUM = 11267

# server and clients exchange lists of strings

def send_strings(sock, bufsize, items):
    """
    A list of strings in items will be sent
    via the socket sock using buffer size bufsize.
    First the number of strings is sent,
    followed by the strings in items.
    """
    dim = len(items)
    strdim = str(dim)
    data = strdim + (bufsize-len(strdim))*' '
    sock.send(data)
    for i in range(0, dim):
        data = items[i] + (bufsize-len(items[i]))*' '
        sock.send(data)

def recv_strings(sock, bufsize):
    """
    A list of strings is received
    via the socket sock using buffer size bufsize.
    First the number of strings is received,
    before the strings.
    On return is the list of strings.
    """
    result = []
    data = sock.recv(bufsize).strip()
    try:
        dim = int(data)
    except ValueError:
        print('no integer received')
        return result
    for i in range(0, dim):
        data = sock.recv(bufsize).strip()
        result.append(data)
    return result

def solve_system(pols):
    """
    Calls the black box solver of PHCpack
    for valid input lists pols.
    Returns the solution found.
    """
    if len(pols) > 0:
        from phcpy import solver
        return solver.solve(pols)
    else:
        return []

# code for the server :

from socket import socket, AF_INET, SOCK_STREAM

def server_connect(nbclients, bufsize, portnum):
    """
    Connects a server to listen to nbclients clients,
    using buffer size bufsize and port number p.
    Returns the socket to connect clients.
    """
    hostname = ''      # to use any address
    number = portnum   # number for the port
    # buffer = bufsize   # size of the buffer
    server_address = (hostname, number)
    sock = socket(AF_INET, SOCK_STREAM)
    sock.bind(server_address)
    sock.listen(nbclients)
    return sock

from threading import Thread

class ServerHandler(Thread):
    """
    Defines the action of the handler threads,
    using static workload balancing.
    """
    def __init__(self, n, s, m, problems):
        """
        The data attributes of each handler are
        n : the name of the thread;
        s : the server socket;
        m : the number of handler threads;
        problems : the list of problems to solve.
        """
        print('initialization of thread ', n)
        Thread.__init__(self, name=n)
        self.svs = s
        self.nht = m
        self.sys = problems 

    def run(self):
        """
        Each handler accepts a connection from a client.
        Thread k send those problems in the list L
        whose index modulo the number of threads is k.
        """
        tname = self.getName()
        print('running of thread ', tname)
        server = self.svs
        client, address = server.accept()
        print(tname + ' accepted request from ', address)
        index = int(tname)
        while index < len(self.sys):
            send_strings(client, BUFSIZE, self.sys[index])
            sols = recv_strings(client, BUFSIZE)
            print('received %d solutions for problem %d' % (len(sols), index))
            # for i in range(0, len(sols)): print sols[i]
            index = index + self.nht
        send_strings(client, BUFSIZE, [])

def start_server(systems, nbclients):
    """
    The server has a list of systems in systems
    for nbclients clients to solve.
    """
    sevsck = server_connect(nbclients, BUFSIZE, PORTNUM)
    print('server is ready for %d clients' % nbclients)
    servers = [ServerHandler(str(i), sevsck, nbclients, systems) \
               for i in range(0, nbclients)]
    print('server starts %d threads' % len(servers))
    for aserver in servers:
        aserver.start()
    # sevsck.close()

# code for the client

def client():
    """
    The client connects to a server and
    solves a system.
    """
    sock = socket(AF_INET, SOCK_STREAM)
    try:
        sock.connect(('localhost', PORTNUM))
    except:
        print('could not connect to server')
        return
    print('connected to server')
    from phcpy import solver
    while True:
        pols = recv_strings(sock, BUFSIZE)
        print('received ', pols)
        if pols == []:
            break
        sols = solver.solve(pols)
        send_strings(sock, BUFSIZE, sols)
    sock.close()

def main():
    """
    Launches one server and clients on the
    same computer.
    """
    clorsv = input("client or server ? (type c or s) ")
    if clorsv == 'c':
        client()
    else:
        from phcpy import solver
        if clorsv == 's':
            nbsys = eval(input('-> how many systems ? '))
            probs = [solver.random_trinomials() for i in range(0, nbsys)]
            nbclt = eval(input('-> how many clients ? '))
            start_server(probs, nbclt)
        else:
            print('sorry, expected c or s, please try again')

if __name__ == "__main__":
    main()
