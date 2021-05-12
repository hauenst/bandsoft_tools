import sys
import re

if len(sys.argv) != 3 :
    print('Incorrect number of arguments. Instead use:')
    print('\tpython code.py [input option 1 = spring2019, 2 = fall2019] [output.txt]')
    exit(-1)


if sys.argv[1] == '1' :
    print('Create table for spring2019')
elif sys.argv[1] == '2' :
    print('Create table for fall2019')
else:
    print('Option for table not known')
    exit(-1)

output_file = open(sys.argv[2],'w')

layer = 6
sector = 5
#number of components for each layer and sector. Overall 30 numbers
component = [ [3,7,6,6,2],
    [ 3,7,6,6,2 ],
    [ 3,7,6,6,2 ],
    [ 3,7,6,6,2 ],
    [ 3,7,5,5,0 ],
    [ 3,7,6,6,2 ] ]

#bad PMTS Spring2019
#sector*1000+layer*100+component*10+Order (left = 0 and right = 1)
badspring19 = [3410 , 3411, 4320, 4321, 4251  ]
badfall19 = [3410 , 3411, 4320, 4321, 4251 , 3261, 4161 ]
tempbad = []
if sys.argv[1] == '1' :
    tempbad = badspring19
elif sys.argv[1] == '2' :
    tempbad = badfall19
else:
    print('Option for table not known')
    exit(-1)

output_file.write("#Table for BAND status words. Sector, Layer, Component, Status PMT Left, Status PMT Right"+ '\n')

for ilay in range(1,layer+1) :
    for isec in range(1, sector+1) :
    #    print (str(isec) + " , " + str(ilay) + "  , " + str(component[ilay-1][isec-1]))
        for icomponent in range(1, component[ilay-1][isec-1]+1 ) :
            print (str(isec) + " , " + str(ilay) + "  , " + str(icomponent))
            leftid = isec*1000 + ilay*100 + icomponent*10
            rightid = isec*1000 + ilay*100 + icomponent*10 + 1
            leftstatus = 0
            rightstatus = 0
    #Check if id or PMT is in bad list. Set status from default 0 to 5
            if any([i == leftid for i in tempbad] ) :
                leftstatus = 5
                print("found bad PMT number " + str(leftid))
            if any([i == rightid for i in tempbad] ) :
                rightstatus = 5
                print("found bad PMT number " + str(rightid))
            output_file.write(str(isec) + " " + str(ilay) + " " + str(icomponent) + " " + str(leftstatus) + " " + str(rightstatus) + '\n')

output_file.close()
