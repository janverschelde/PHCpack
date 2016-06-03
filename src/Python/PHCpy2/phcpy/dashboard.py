"""
This module protypes a graphical user interface to phcpy.
"""

from Tkinter import Tk, Message, Button, W, E, END, Label
from Tkinter import Entry, INSERT, Text
from tkFont import Font
from Tkinter import StringVar

class Scroller(object):
    """
    Scrolls through a solution list.
    """
    def __init__(self, wdw, sols):
        """
        Stores the list of solutions in sols
        and defines the layout of the GUI.
        """
        wdw.title('solutions scroller')
        self.sols = sols
        self.cursor = 0
        self.lbl = Label(wdw, text="solution : ")
        self.lbl.grid(row=0, column=0, sticky=E) 
        self.ent = Entry(wdw)
        self.ent.grid(row=0, column=1, stick=W)
        self.ent.insert(INSERT, "0 of %d" % len(sols))
        self.myft = Font(family="Courier New", size=12, weight="normal")
        self.mlen = self.myft.measure("M")
        lines = sols[0].split('\n')
        self.width = max([len(line) for line in lines])
        self.display = StringVar()
        self.display.set(self.sols[0])
        self.mess = Message(wdw, textvariable=self.display, \
            font=self.myft, width=self.width*self.mlen, background='white')
        self.mess.grid(row=1, column=0, columnspan=2)
        self.btnext = Button(wdw, command=self.next, text='next')
        self.btnext.grid(row=2, column=1, sticky=W+E)
        self.btprev = Button(wdw, command=self.previous, text='previous')
        self.btprev.grid(row=2, column=0, sticky=W+E)

    def show(self):
        """
        Shows the solution at position self.cursor
        in the message widget and updates the entry widget.
        """
        self.display.set(self.sols[self.cursor])
        self.ent.delete(0, END)
        self.ent.insert(INSERT, '%d of %d' % (self.cursor, len(self.sols)))

    def next(self):
        """
        Increases the cursor by one if possible.
        """
        if self.cursor < len(self.sols) - 1:
            self.cursor = self.cursor + 1
        self.show()

    def previous(self):
        """
        Decreases the cursor by one if possible.
        """
        if self.cursor > 0:
            self.cursor = self.cursor - 1
        self.show()

def scrollsols(sols):
    """
    Instantiates the window and launches the GUI
    to scroll through the solutions in the lis sols.
    """
    wdw = Tk()
    Scroller(wdw, sols)
    wdw.mainloop()

def testscroller():
    """
    Solves the cyclic 5-roots problems
    and launches the scroller.
    """
    from phcpy.families import cyclic
    from phcpy.solver import solve
    cyc5 = cyclic(5)
    sols = solve(cyc5, silent=True)
    scrollsols(sols)

def pols2str(pols):
    """
    Returns the input string to put into the text input widget 
    for the string representations of polynomials in pols.
    """
    result = str(len(pols)) + '\n'
    for pol in pols:
        result = result + pol + '\n'
    return result

def str2pols(strp):
    """
    Returns the list of string representations of the polynomials
    in the string strp.
    """
    data = strp.split('\n')
    dim = int(data[0])
    result = [data[k] for k in range(1, dim+1)]
    return result

class SolveButton(object):
    """
    Simple input/output text widget and button
    graphical user interface to a solver.
    """
    def __init__(self, wdw, pols=[]):
        """
        The layout consists of two labels,
        two text widgets and one button.
        The optional input is the list of polynomials.
        """
        wdw.title('solve a polynomial system')
        self.myft = Font(family="Courier New", size=12, weight="normal")
        self.inputd = Text(wdw, width=80, height=10, font=self.myft)
        if pols != []:
            self.inputd.insert(END, pols2str(pols))
        self.output = Text(wdw, width=80, height=20, font=self.myft)
        self.submit = Button(wdw, text="solve the system", \
            command=self.solve, font=self.myft)
        self.inlbl = Label(wdw, text="enter a polynomial system:", \
            font=self.myft)
        self.outlbl = Label(wdw, text="scroll solutions with up/down arrows", \
            font=self.myft)
        self.inlbl.grid(row=0, column=0, sticky=W)
        self.inputd.grid(row=1, column=0)
        self.submit.grid(row=2, column=0, sticky=W+E)
        self.output.grid(row=3, column=0)
        self.outlbl.grid(row=4, column=0)

    def solve(self):
        """
        Takes the data from the input text widget,
        places the data into a list of polynomials,
        calls the blackbox solver and displays the
        solutions in the output text widget.
        """
        self.output.delete(1.0, END)
        data = self.inputd.get(1.0, END)
        self.output.insert(END, data)
        psys = str2pols(data)
        from phcpy.solver import solve
        sols = solve(psys, silent=True)
        self.output.insert(END, 'Computed %d solutions.\n\n' % len(sols))
        cnt = 0
        for sol in sols:
            cnt = cnt + 1
            self.output.insert(END, 'solution %d :\n' % cnt)
            self.output.insert(END, sol + '\n')

def launchsolver(pols=[]):
    """
    Instantiates a Tk object
    and launches the event loop.
    """
    wdw = Tk()
    SolveButton(wdw, pols)
    wdw.mainloop()

def testsolvebutton():
    """
    Solves the cyclic 5-roots problems
    and launches the solve button.
    """
    from phcpy.families import cyclic
    cyc5 = cyclic(5)
    launchsolver(cyc5)

if __name__ == "__main__":
    # testscroller()
    testsolvebutton()
