"""
Python 3 script with matplotlib and numpy to illustrate
a motivating example for the Pade approximation:
x(t) = sqrt((1 + 1/2*t)/(1 + 2*t)).
The [2/2] Pade approximant is
(9.53125*t^2 + 2.125*t + 1)/(1.890625*t^2 + 2.875*t + 1)
but that is too close to x(t) for a good plot.
Instead, use the [1/1] Pade approximant
(0.875*t + 1)/( + 1.625*t + 1)
with corresponding series
1 - 0.75*t + 1.21875*t^2.
"""
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
axs = fig.add_subplot(111)
axs.set_xlabel('t')
axs.set_ylabel('x')
axs.set_xlim([-0.1, 4.1])
axs.set_ylim([0.50, 1.1])
tvl = np.arange(0, 4.01, 0.01)
fun = lambda x: np.sqrt((1 + 0.5*x)/(1 + 2*x))
# pad = lambda t: (0.953125*t**2 + 2.125*t + 1)/(1.890625*t**2 + 2.875*t + 1)
pad = lambda t: (0.875*t + 1)/(1.625*t + 1)
ser = lambda t: 1.0 - 0.75*t + 1.21875*t**2
yvl = fun(tvl)
zvl = [pad(t) for t in tvl]
svl = [ser(t) for t in tvl]
fp = axs.plot(tvl, yvl, 'b-')
pp = axs.plot(tvl, zvl, 'r-')
sp = axs.plot(tvl, svl, 'g-')
axs.legend((fp[0], sp[0], pp[0]), 
('sqrt((1 + 0.5*x)/(1 + 2*x))', '1.0 - 0.75*t + 1.21875*t**2', \
 '(0.875*t + 1)/(1.625*t + 1)'))
axs.set_title('rational versus series approximation of a function')
plt.show()
