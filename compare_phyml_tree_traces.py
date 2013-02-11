#
# Compare two *_trace_ files from PhyML, using different optimization algorithms or different models.
# Did the two searches explore similar regions of tree space?
#

import re,sys,os,pickle
from dendropy import Tree
from argParser import ArgParser
argParser = ArgParser(sys.argv)

t1path = sys.argv[1]
t2path = sys.argv[2]

t1trees = []
t2trees = []

HUGE = 100000000; # approaching infinity

colors = {}
colors["Multimax"] = "royalblue"
colors["Unimax"] = "orangered"
colors[1] = "black"
colors["UMtd"] = colors["Unimax"] 
colors["MMtd"] = colors["Multimax"] 

pch = {}
pch["Multimax"] = "15"
pch["Unimax"] = "15"
pch[1] = "15"
pch["UMtd"] = "17" 
pch["MMtd"] = "16"

#
# points is an array of arrays, points[row] = column
#
def draw_heatmap(points, rowlabels, collabels, outputfilenameseed):
    dataseg = "c("
    for i in points:
        for d in i:
            dataseg += d.__str__() + ","
    dataseg = re.sub(",$", "", dataseg)
    dataseg += ")\n"
    
    rownames = "c("
    for i in rowlabels:
        rownames += i.__str__() + ","
    rownames = re.sub(",$", "", rownames)
    rownames += ")\n"

    colnames = "c("
    for i in collabels:
        colnames += i.__str__() + ","
    colnames = re.sub(",$", "", colnames)
    colnames += ")\n"
                    
    # build a CRAN script
    fout = open("./" + outputfilenameseed + ".cran", "w")
    fout.write("library(RColorBrewer)\n")
    fout.write("tm<-matrix(" + dataseg + ",ncol=" + collabels.__len__().__str__() + ",byrow=TRUE)\n")
    fout.write("rownames(tm)<-" + rownames + "\n")
    fout.write("colnames(tm)<-" + colnames + "\n")
    fout.write("pdf(\"" + outputfilenameseed + ".pdf\")\n")
    #fout.write("h <- heatmap(tm, col=brewer.pal(32,\"YlGnBu\"), scale=\"none\", margins=c(5,10), Rowv=NA, Colv=NA)\n")
    fout.write("h <- heatmap(tm, col=colorRampPalette(brewer.pal(9,\"YlGn\"))(256), scale=\"none\", margins=c(5,10), Rowv=NA, Colv=NA)\n")
    fout.write("dev.off()\n")
    fout.close()
    os.system("r --save < " + "./" + outputfilenameseed + ".cran",)

#
# Write CRAN scripts and plot PDFs.
#
# points[dataset] = {} where key = xvalue, value = yvalue
def plot_in_r(points, output_filename_seed, title, xlab, ylab):
    maxx = None
    maxy = None
    miny = None

    filepaths = points.keys()
    filepaths.sort()

    # the ingredients for our CRAN script:
    plotstring = ""
    pointsstring = ""
    legendstring = "legend(\"bottomright\", c("
    for f in points.keys():
        legendstring += "\"" + f.__str__() + "\","
    legendstring = re.sub(",$", "", legendstring)
    legendstring += "), cex=0.8, col=c("
    for f in points.keys():
        legendstring += "\"" + colors[f].__str__() + "\","
    legendstring = re.sub(",$", "", legendstring)
    legendstring += "), pch=c("
    for f in points.keys():
        legendstring +=  pch[f].__str__() + ","
    legendstring = re.sub(",$", "", legendstring)    
    legendstring += "));\n"

    count = -1
    for f in points.keys():
        count += 1
        x_sorted = points[f].keys()
        x_sorted.sort()
        
        string = "x" + count.__str__() + "<-c("
        for x in x_sorted:
            if maxx == None:
                maxx = x
            elif x > maxx:
                maxx = x
            string += x.__str__() + ","   
        string = re.sub(",$", "", string)
        string += ");\n"
        
        string += "y" + count.__str__() + "<-c("
        for x in x_sorted:
            if miny == None:
                miny = points[f][x]
            elif points[f][x] < miny:
                miny = points[f][x]
            if maxy == None:
                maxy = points[f][x]
            elif points[f][x] > maxy:
                maxy = points[f][x]
            string += points[f][x].__str__() + ","
        string = re.sub(",$", "", string)
        string += ");\n"
        if f.__str__().__contains__("td"):
            string += "points(x" + count.__str__() + ",y" + count.__str__() + ", col='" + colors[f] + "', type='p',pch=" + pch[f].__str__() + ");\n"
        else:
            string += "points(x" + count.__str__() + ",y" + count.__str__() + ", col='" + colors[f] + "', type='l',lwd=3);\n"
        pointsstring += string
    
    plotstring = "x <-c(0.0," + maxx.__str__() + ");\n"
    plotstring += "y <-c(" + miny.__str__() + "," + maxy.__str__() + ");\n"
    plotstring += "plot(x, y, type='n', main='" + title + "', xlab='" + xlab + "', ylab='" + ylab + "');\n"
    
    fout_cran = open(output_filename_seed + ".cran", "w")
    fout_cran.write("pdf('" + output_filename_seed + ".pdf', width=7, height=4);\n")
    fout_cran.write(plotstring)
    fout_cran.write(pointsstring)
    fout_cran.write(legendstring)
    fout_cran.write("dev.off();\n")
    fout_cran.close()
    
    os.system("r --no-save < " + output_filename_seed + ".cran")


#
# Returns [an array of unique trees in the path, an array of lnls for those trees]
#
def return_trees_from_trace(path):
    print "Parsing trace:", path
    trees = []
    lnls = []
    fin = open(path, "r")
    last_tree = None
    last_lnl = 0.0
    count_unique_trees = 0
    for line in fin.xreadlines():
        treestring = ""
        lnlstring = ""
        found_tree = False
        for c in line:
            if found_tree == False and c != "]" and c != "[" and c != "(":
                lnlstring += c
            if c == "(":
                found_tree = True
            if found_tree == True:
                treestring += c
        lnl = float(lnlstring)
        t = Tree()
        t.read_from_string(line, "newick")
        if last_tree != None: #2nd->nth trees in the list
            #sd = last_tree.symmetric_difference(t)
            #sd = t.symmetric_difference(last_tree)
            if last_lnl < lnl:
                trees.append(t)
                lnls.append("%.2f"%lnl)
                count_unique_trees += 1
            else:
                trees[trees.__len__()-1] = t
                lnls[lnls.__len__()-1] = "%.2f"%lnl
        else: #first tree in the list
            trees.append(t)
            lnls.append("%.2f"%lnl)
            count_unique_trees += 1
        last_tree = t
        last_lnl = lnl
        print count_unique_trees, lnl
    trees.append(last_tree)
    lnls.append("%.2f"%lnl)
    fin.close()
    return [trees, lnls]

def build_distance_matrix(t1trees, t2trees):
    points = []
    for i in t1trees:
        p = []
        for j in t2trees:
            #p.append(i.symmetric_difference(j) )
            p.append(j.symmetric_difference(i) )
        points.append(p)
        print p
    return points

def get_my_cost(i, j, d, c):
    n = d.__len__()
    m = d[0].__len__()
    if i == (n-1) and j == (m-1):
        return d[i][j]
    elif i == (n-1):
        return c[i][j+1]
    elif j == (m-1):
        return c[i+1][j]
    else:
        above = c[i+1][j]
        right = c[i][j+1]
        aboveright = c[i+1][j+1]
        #print above, right, aboveright
        # return my cost plus the minimum cost available from my neighbors 
        return d[i][j] + min(min(above,right),aboveright)

# d = distance matrix
def build_cost_matrix(d):
    c = []
    n = d.__len__()
    m = d[0].__len__()
    
    # init the matrix to have very large numbers
    for i in range(0,n):
        this_c = []
        for j in range(0,m):
            this_c.append(HUGE)
        c.append(this_c)
    for i in reversed(range(0,n)):
        for j in reversed(range(0,m)):
            c[i][j] = get_my_cost(i, j, d, c)
    print "Cost matrix = "
    for i in range(0,c.__len__()):
        print c[i]
    return c

# return [i,j] for minimum neighbor
def get_min_neighbor_index(c, i, j):    
    n = c.__len__()
    m = c[0].__len__()
    if i == (n-1) and j == (m-1):
        return None
    elif i == (n-1):
        return [i,j+1]
    elif j == (m-1):
        return [i+1,j]
    else:
        above = c[i+1][j]
        right = c[i][j+1]
        aboveright = c[i+1][j+1]
        print above, right, aboveright
        if above < right and above < aboveright:
            return [i+1,j]
        elif right < aboveright:
            return [i,j+1]
        else:
            return[i+1,j+1]

# c = cost matrix
# returns an array of tuples[(,),(,),(,)...]
def find_min_path(c):
    path = []
    n = c.__len__()
    m = c[0].__len__()
    
    path.append([0,0])
    next = get_min_neighbor_index(c,0,0)
    path.append(next)
    while (next != None):
        next = get_min_neighbor_index(c, next[0], next[1])
        path.append(next)
    print "min path = ", path
    return path

    
def write_trace(trees, lnls, outputpath):
    fout = open(outputpath, "w")
    for i in range(0, trees.__len__()):
        fout.write("[" + lnls[i].__str__() + "]" + trees[i].__str__() + "\n")
    fout.close()

# d = the distance matrix, i.e. "points"
def plot_aligned_distances(min_path, d):
    points = {} # we'll send points to the method plot_in_r
    count = 0
    for p in min_path:
        if p == None:
            continue
        #print p
        points[ count ] = d[ p[0] ][ p[1] ]
        count += 1
    plot_in_r({1:points}, "plot_aligned_distance", "Aligned Tree Distance", "trees", "d")

def plot_aligned_lnls(min_path, t1lnls, t2lnls):
    points = {}
    points["Unimax"] = {}
    points["Multimax"] = {}
    count_i = 0
    for i in min_path:
        if i == None:
            continue
        n = i[0]
        m = i[1]
        points["Unimax"][count_i] = t1lnls[n]
        points["Multimax"][count_i] = t2lnls[m]
        count_i += 1
    plot_in_r(points, "plot_aligned_lnls", "lnL(tree steps)", "topology swap steps", "lnL")       

def plot_lnls(t1lnls, t2lnls,di=None):
    points = {}
    points["Unimax"] = {}
    points["Multimax"] = {}
    for i in range(0, t1lnls.__len__()):
        points["Multimax"][i] = t1lnls[i]  
    for i in range(0, t2lnls.__len__()):
        points["Unimax"][i] = t2lnls[i]
    if di != None:
        if False == points.keys().__contains__("UMtd"):
            points["UMtd"] = {}
        if False == points.keys().__contains__("MMtd"):
            points["MMtd"] = {}
        points["UMtd"][di] = t2lnls[di]
        points["MMtd"][di] = t1lnls[di]
    plot_in_r(points, "plot_lnls", "lnL(T,b,m|d)", "trees", "lnL") 

def write_tree(tree, tpath):
    fout = open(tpath, "w")
    fout.write( tree.__str__() + ";\n")
    fout.close()

def write_all_trees(unitrees, multitrees):
    count = 0
    for t in unitrees: 
        path = "tree.unimax." + count.__str__() + ".tre"
        count += 1
        write_tree( t, path)
    count = 0
    for t in multitrees:
        path = "tree.multimax." + count.__str__() + ".tre"
        write_tree( t, path )
        count += 1

def plot_treelength(min_path, t1trees, t2trees):
    points = {}
    points["Unimax"] = {}
    points["Multimax"] = {}
    count_i = 0
    for i in min_path:
        if i == None:
            continue
        n = i[0]
        m = i[1]
        points["Unimax"][count_i] = t1trees[n].length()
        points["Multimax"][count_i] = t2trees[m].length()
        count_i += 1        
        #os.system("python /Users/victor/Documents/EclipseWorkspace/PhyloWidgets/plot_bl_distribution.py " + t1path + " " + t2path + " --id " + i.__str__())
        
    plot_in_r(points, "plot_treelength", "tree length", "topology swap steps", "tree length")    

#
# Returns td (as a Tree object), and four lnls: UM on td, UM on td+1, MM on td, MM on td+1
#
def get_td(t1trees, t1lnls, t2trees, t2lnls):
    dindex = 0
    td = Tree()
    umtd1 = Tree()
    mmtd1 = Tree()
    l_umtd = None
    l_umtd1 = None
    l_mmtd = None
    l_mmtd1 = None
    for i in range(0, t1trees.__len__()):
        t1 = t1trees[i]
        t2 = t2trees[i]
        d = t1.symmetric_difference(t2)
        d = t2.symmetric_difference(t1)
        print "d", d
        if d == 0:
            td = t1
            l_umtd = t1lnls[i]
            l_mmtd = t2lnls[i]
            dindex = i
            print "ML td", l_umtd, l_mmtd
        else:
            umtd1 = t1
            mmtd1 = t2
            break
    
    return [dindex, td, umtd1, mmtd1]  

################################
#
# main:
#

points = None
if False != argParser.getOptionalArg("--saved"):
    fin = open( argParser.getOptionalArg("--saved"), "r")
    [points, t1lnls, t2lnls, t1trees, t2trees] = pickle.load( fin )
    fin.close()
else:
    [t1trees, t1lnls] = return_trees_from_trace(t1path)
    [t2trees, t2lnls] = return_trees_from_trace(t2path)
    write_all_trees(t1trees, t2trees)
    write_trace(t1trees, t1lnls, t1path + ".unique.txt")
    write_trace(t2trees, t2lnls, t2path + ".unique.txt")    
    points = build_distance_matrix(t1trees, t2trees)
    fout = open("points.txt", "w")
    pickle.dump([points, t1lnls, t2lnls, t1trees, t2trees], fout)
    fout.close()
draw_heatmap(points, t1lnls, t2lnls, "compare_phyml_tree_traces")



print "Finding the divergence tree..."
[dindex, td, umtd1, mmtd1] = get_td(t1trees, t1lnls, t2trees, t2lnls)
write_tree(td, "td.tre")
write_tree(umtd1, "umtd+1.tre")
write_tree(mmtd1, "mmtd+1.tre")
plot_lnls(t1lnls,t2lnls,di=dindex)

print "calculating the cost matrix = "
cost_matrix = build_cost_matrix(points)
draw_heatmap(cost_matrix, t1lnls, t2lnls, "compare_phyml_tree_traces.cost")
min_path = find_min_path(cost_matrix)
plot_aligned_distances(min_path, points)
plot_aligned_lnls(min_path, t1lnls, t2lnls)
plot_treelength(min_path, t1trees, t2trees)