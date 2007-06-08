from dolfin import *
from pylab import *

if (size(sys.argv) < 4):
  print 'Usage:', sys.argv[0], '<title> <objective> <plotfile> [plotfile 2 ... plotfile n]'
  print """title - The plot title
           objective - What to measure. Legal values are "time", "edgecut" 
                       and "balance".
           plotfile - Output file from the benchmark program
        """
  sys.exit(1)

plottitle = sys.argv[1]
objective = sys.argv[2]

if objective == "time":
  pos = 1
elif objective == "edgecut":
  pos = 2
elif objective == "balance":
  print 'not implemented'
  sys.exit(1)

for filename in sys.argv[3:]:
  file = open(filename, 'r')
  lines = file.readlines();
  linelabel = lines[0]
  axis = lines[1].split()
  (xname, yname) = (axis[0], axis[pos])
  xarr = []
  yarr = []
  for line in lines[2:]:
    xarr.append(line.split()[0])
    yarr.append(line.split()[pos])

  loglog(xarr, yarr, '-o', label=linelabel)
  legend(loc='best')
  file.close()

title(plottitle)
xlabel(xname)
ylabel(yname)
show()
