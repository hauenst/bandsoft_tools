import sys
import re

if len(sys.argv) != 3 :
    print('Incorrect number of arguments. Instead use:')
    print('\tpython code.py [input.txt clas12tags] [output.txt for reconstruction]')
    exit(-1)


sector = []
layer = []
component = []
y = []
z = []
with open(sys.argv[1],'r') as f:
    for line in f:
        line = line.strip()
        #print(line)
        if re.match("scintillator", line):
            #split line on | chars
            parse = line.split("|")
            strippedparse = []
            for item in parse:
                #strip unnecessary white space in the list entries
                item = item.strip()
                strippedparse.append(item)
            #copy position (x,y,z) string and bar id string (sector, layer, component)
            barpos_string = strippedparse[3]
            barid_string = strippedparse[17]
            #1. split string with x,y,z pos and get y and z-pos string from splitted list (array entry 1 and 2).
            #2. Remove *mm from y-pos and z-pos string and convert to float
            y.append(float( barpos_string.split(" ")[1].rstrip("*mm")))
            z.append(float( barpos_string.split(" ")[2].rstrip("*mm")))
            #sector is 3rd list entry (array 2)
            sector.append(int( barid_string.split(" ")[2]))
            #layer is 6rd list entry (array 5)
            layer.append(int( barid_string.split(" ")[5]))
            #sector is 9rd list entry (array 8)
            component.append(int( barid_string.split(" ")[8]))
            #prints for debugging
#            print(barpos_string + " " + barid_string)
#            print(y)
#            print(z)
#            print(sector)
#            print(layer)
#            print(component)
        elif re.match("band", line):
            #split line on | chars
            parse = line.split("|")
            motheritems = []
            for item in parse:
                #strip unnecessary white space in the list entries
                item = item.strip()
                motheritems.append(item)
            #copy mother volume position (x,y,z) string
            motherpos = motheritems[3]
            #1. split string with x,y,z pos and get z-pos string from splitted list (array entry 2).
            #2. Remove *mm from z-pos string and convert to float
            motherzpos = float( motherpos.split(" ")[2].rstrip("*mm"))
            print("Mother volume z-pos in mm: " + str(motherzpos))
        else:
            pass

if len(y) == len(z) == len(sector) == len(layer) == len(component):
    print("all lists have same length")
    #open output file
    output_file = open(sys.argv[2],'w')
    for i in range(0,len(y)):
        #add BAND mother volume shift to each z-position for global z value
        shiftedz = motherzpos - z[i]
        #convert into cm
        shiftedz/= 10.
        y[i]/= 10.
        print(str(sector[i]) + " " + str(layer[i]) +  " " + str(component[i]) + " " + str(y[i]) + " " + str(shiftedz))
        output_file.write(str(sector[i]) + " " + str(layer[i]) +  " " + str(component[i]) + " " + str(y[i]) + " " + str(shiftedz) + '\n')

    output_file.close()
