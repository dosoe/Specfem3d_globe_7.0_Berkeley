g=open('Crust.asc','r')
lines_in=g.readlines()
g.close()

f=open('Crust.asc','w')
i=0
while i<2:
    f.write(lines_in[i])
    i+=1

while i<len(lines_in):
    f.write(lines_in[i])
    f.write('0.00000    2.61873    5.50377    3.25969   99900.00     300.00    5.50377    3.25169    1.00000\n')
    f.write('5.18019    2.59992    5.39792    3.20202   99900.00     300.00    5.39792    3.19127    1.00000\n')
    f.write('15.00000    2.65398    5.69225    3.36430   99900.00     300.00    5.69225    3.34992    1.00000\n')
    f.write('24.81981    3.01883    7.17236    4.11992   99900.00     300.00    7.17236    4.11882    1.00000\n')
    f.write('30.00000    3.16208    7.62491    4.34962   99900.00     300.00    7.62491    4.35490    1.00000\n')
    i+=6
f.close()

