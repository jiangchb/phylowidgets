import math,re,sys,os
fin = open(sys.argv[1], "r")
newline = ""
for line in fin.readlines():
    if line.__len__() > 10:
        tokens = line.split(":")
        for t in tokens:
            #print t
            if t.__contains__(")") and False == t.__contains__(";"):
                #print t
                ts = t.split(")")
                newline += ts[0] + ")"
                alrt = float(ts[1])
                if alrt > 1418:
                    alrt = 1418.0
                alr = math.exp(alrt/2.0)
                #print alrt, alr
                newline += "%.4e"%alr
                newline += ":"
            else:
                newline += t
                newline += ":"
        print newline
fin.close()