import numpy as np
import matplotlib.pyplot as plt

s_TDC = {}
s_FTDC = {}
f_TDC = {}
f_FTDC = {}
with open("../../include/calibrations/spring2019/lr_offsets.txt") as f:
    for line in f:

        parse = line.strip().split()
        sector = int(parse[0])
        layer = int(parse[1])
        comp = int(parse[2])

        tdc_vel = float(parse[4])
        ftdc_vel = float(parse[7])

        barID = sector*100+layer*10+comp
        s_TDC[barID] = tdc_vel
        s_FTDC[barID] = ftdc_vel

with open("../../include/calibrations/fall2019/lr_offsets.txt") as f:
    for line in f:

        parse = line.strip().split()
        sector = int(parse[0])
        layer = int(parse[1])
        comp = int(parse[2])

        tdc_vel = float(parse[4])
        ftdc_vel = float(parse[7])

        barID = sector*100+layer*10+comp
        f_TDC[barID] = tdc_vel
        f_FTDC[barID] = ftdc_vel

for key in s_TDC:
    if int(int(key - int(key/100)*100)/10)==6: continue
    if( s_TDC[key]==0 ): print(key);continue;
    if( f_TDC[key]==0 ): print(key);continue;
    plt.scatter(key,s_TDC[key],color='red')
    plt.scatter(key,s_FTDC[key],color='blue')
    plt.scatter(key,f_TDC[key],color='red',facecolor='none')
    plt.scatter(key,f_FTDC[key],color='blue',facecolor='none')

plt.xlim([75,575])
plt.axvline(x=180,linestyle='--',color='black')
plt.axvline(x=280,linestyle='--',color='black')
plt.axvline(x=380,linestyle='--',color='black')
plt.axvline(x=480,linestyle='--',color='black')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Bar ID',fontsize=16)
plt.ylabel('Effective Velocity [cm/ns]',fontsize=16)
plt.text(100,11.5,'TDC',color='red',fontsize=16)
plt.text(100,11,'FADC',color='blue',fontsize=16)

plt.savefig('effective_vel_compare.pdf')
plt.show()
