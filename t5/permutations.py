from itertools import permutations
import cmath
import math
import sys


def updateres(data):
	global res
	a=0
	for re in res:
		if re[0] == data[0] and re[1] == data[1] and re[2] == data[2] and re[3] == data[3]:
			a=1

	if a==0:
		res.append(data[:])

def combinationUtilres(arr, data, start,
                    end, index, r):  

    if (index == r):
        updateres(data)
        return;
 

    i = start;
    while(i <= end and end - i + 1 >= r - index):
        data[index] = arr[i];
        combinationUtilres(arr, data, i + 1,
                        end, index + 1, r);
        i += 1;




def updatecon(data):
	global con
	a=0
	for co in con:
		if co[0] == data[0] and co[1] == data[1]:
			a=1

	if a==0:
		con.append(data[:])

def combinationUtilcon(arr, data, start,
                    end, index, r):  

    if (index == r):
        updatecon(data)
        return;
 

    i = start;
    while(i <= end and end - i + 1 >= r - index):
        data[index] = arr[i];
        combinationUtilcon(arr, data, i + 1,
                        end, index + 1, r);
        i += 1;
 

global res
global con
res=[]
arr = [1000,1000,1000,10000,10000,10000,100000,100000,100000,2000,3000,20000,30000,200000,300000];
r = 4;
n = len(arr);

datares = [0]*r;
combinationUtilres(arr, datares, 0,
                    n - 1, 0, r);


con = []
arrr = [0.000000220,0.000000220,0.000000220,0.000001,0.000001,0.000001,0.000000110,0.0000005];
r1 = 2;
n = len(arrr);
datacon = [0]*r1;
combinationUtilcon(arrr, datacon, 0,
                    n - 1, 0, r1);

sys.stdout = open('t5_file.txt', 'w')



max_merit = 0
min_gain_dev=9999
min_freq_dev=9999

for re in res:
	perm_re = list(permutations(re))
	for per_re in perm_re:
		R1 = per_re[0]
		R2 = per_re[1]
		R3 = per_re[2]
		R4 = per_re[3]

		for co in con:
			perm_co = list(permutations(co))
			for per_co in perm_co:
				C1 = per_co[0]
				C2 = per_co[1]

				wL = 1/(R1*C1)
				wH = 1/(R2*C2)
				wO = math.sqrt(wL*wH)
				f = wO/(2*math.pi)

				gain = abs((R1*C1*wO*complex(0,1))/(1+R1*C1*wO*complex(0,1))*(1+R3/R4)*(1/(1+R2*C2*wO*complex(0,1))))

				if C2 == 0.000000110:
					C2 = C2*4

				if C1 == 0.000000110:
					C1 = C2*4

				Cost = (R1+R2+R3+R4)/1000 + (C1+C2)*1000000 + 13323
				gain_deviation = abs(100-gain)
				frequency_deviation = abs(f-1000)
				Merit = 1/(Cost*(gain_deviation+frequency_deviation+10**(-6)))
				
				if Merit > max_merit:
					max_merit = Merit

				if gain_deviation < min_gain_dev:
					min_gain_dev = gain_deviation

				if frequency_deviation < min_freq_dev:
					min_freq_dev = frequency_deviation


				print('R1 = %f'% R1)
				print('R2 = %f'% R2)
				print('R3 = %f'% R3)
				print('R4 = %f'% R4)
				print('C1 = %.9f'% C1)
				print('C2 = %.9f'% C2)
				print('gain = %f'% gain)
				print('freq = %f'% f)
				print('gain_dev = %f'% gain_deviation)
				print('freq_dev = %f'% frequency_deviation)
				print('Cost = %.9f'% Cost)
				print('Merit = %.9f'% Merit)
				print('\n\n')

print()
print(max_merit)
print(min_freq_dev)
print(min_gain_dev)


sys.stdout.close()
