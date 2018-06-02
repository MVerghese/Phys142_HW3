import pylab


data=pylab.loadtxt("out.txt")
pylab.plot(data[:,0],data[:,1])
pylab.legend()
pylab.title("Avg Energy Level per T for Classical Electron")
pylab.xlabel("T")
pylab.ylabel("Avg Energy Level ")
pylab.savefig("output.png")
