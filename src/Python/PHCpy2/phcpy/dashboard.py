"""
This module develops a graphical user interface to phcpy.
"""

from Tkinter import Tk, Message, Button, W, E, END, Label
from Tkinter import Entry, INSERT
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

if __name__ == "__main__":
    testscroller()
