import numpy as np

BAR_LIST = []
with open("/Users/efrainsegarra/work/band/software/bandsoft_tools/include/bar_list.txt") as f:
    for line in f:
        parse = line.strip().split()
        BAR_LIST.append( int(parse[0])*100 + int(parse[1])*10 + int(parse[2]) )


files = ["init_paramsamp","iter0_paramsamp","iter1_paramsamp","iter2_paramsamp"]
parameters = [3,5,2,5]
PARAMETERS = {}
for PMT in ['L','R']:
    for i in range(len(files)):
        fi = files[i]
        orig = {};
        iters = {};
        new = {};
        path = "/Users/efrainsegarra/work/band/documentation/calibrations/spring2019/time_walk/"
        with open(path+fi+PMT+".txt") as f:
            for line in f:
                parse = line.strip().split()
                sector = int(parse[0])
                layer = int(parse[1])
                comp = int(parse[2])
                order = 0 if PMT=='L' else 1

                these_pars = []
                for par in range(parameters[i]):
                    p = float(parse[3+par])
                    these_pars.append(p)
                for empty_par in range(parameters[i],6):
                    these_pars.append(0)
                
                key = sector*1000+layer*100+comp*10+order
                if PARAMETERS.get(key):
                    PARAMETERS[key].append( these_pars )
                else:
                    PARAMETERS[key] = [ these_pars ]

    with open("time_walk_amp_"+PMT+".txt","w") as g:
        for bar in BAR_LIST:
            order = 0 if PMT=='L' else 1
            pars = PARAMETERS[bar*10+order]

            sector = int(bar/100)
            layer = int((bar-sector*100)/10)
            comp = int(bar-sector*100-layer*10)

            line = str(sector) + "\t" + str(layer) + "\t" + str(comp) + "\t"

            for it in range(len(pars)):
                for par in range(len(pars[it])):
                    line += str(pars[it][par]) + "\t" 
            g.write(line+"\n")




                #init[sector*1000+layer*100+comp*10+order] = np.asarray([ p1,p2,p3,p4,p5,p6] )

                #orig[sector*100+layer*10+comp] = np.asarray([ p1, p2, p3, p4, p5, p6 ])
'''
        with open("../../bin/calibration/"+"combo1_"+fi+PMT+".txt","w") as g:
            for key in orig:
                new[key] = orig[key] + iters[key]
                print(new[key])
                sector = str(key)[0]
                layer = str(key)[1]
                comp = str(key)[2]
                g.write(sector+" "+layer+" " +comp+" "+
                        str(new[key][0]) + " " +
                        str(new[key][1]) + " " +
                        str(new[key][2]) + " " +
                        str(new[key][3]) + " " +
                        str(new[key][4]) + " " +
                        str(new[key][5]) + "\n")
'''

