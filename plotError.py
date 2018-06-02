import numpy as np
import matplotlib.pyplot as plt
import pylab

errorData = pylab.loadtxt("error.txt")
plotData = pylab.loadtxt("data.txt")

print(errorData)
print(plotData)

x = []
y = []

for datapoint in plotData:
    x.append(datapoint[0])
    y.append(datapoint[1])

xer = []
yer = []

for datapoint in errorData:
    xer.append(datapoint[0])
    yer.append(datapoint[1])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("Ambient T")
ax.set_ylabel("<E>")
ax.set_title("Averge Energy Given Ambient Temperature for a Quantum Electron")
plt.plot(x,y)
plt.errorbar(x,y,yerr=yer, linestyle="None")

plt.show()