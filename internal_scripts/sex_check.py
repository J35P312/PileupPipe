import sys
import numpy

autosomal_cov=[]
x_cov=[]
y_cov=[]
first=True

female=True
chromosome=False
for line in open(sys.argv[1]):
    if first:
        first=False
        continue

    content=line.strip().split()
    if content[0] == "Y" or content[0] == "chrY":
        y_cov.append(float(content[3]))
    elif content[0] == "X" or content[0] == "chrX":
        x_cov.append(float(content[3]))
    else:
        autosomal_cov.append(float(content[3]))

ylen=0
xlen=0
for line in open(sys.argv[2]):
    content=line.strip().split()
    if content[0] == "Y" or content[0] == "chrY":
        ylen=content[1]
    elif content[0] == "X" or content[0] == "chrX":
        xlen=content[1]
    
x_avg_cov=numpy.average(x_cov)
y_avg_cov=numpy.average(y_cov)
autosomal_avg_cov=numpy.average(autosomal_cov)

#i.e the individual is a male if the coverage of X is about half of the average coverage, and if the coverage of Y is about half of the X coverage
if x_avg_cov/autosomal_avg_cov < 0.7 and y_avg_cov/x_avg_cov > 0.3:
    female=False

if female:

    print "Y	1	{}	*	0".format(ylen)
    print "chrY	1	{}	*	0".format(ylen)
    print "*	*	*	*	2"
else:
    print "Y	1	{}	*	1".format(ylen)
    print "chrY	1	{}	*	1".format(ylen)
    print "X	1	{}	*	1".format(xlen)
    print "chrX	1	{}	*	1".format(xlen)
    print "*	*	*	*	2"
