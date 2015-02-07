## For Chapter six of Bioinfo
import sys

def greedysorting(P):
	'''
        return approxReversalDistance
	'''
	l = len(P)
	f = open('ans.txt', 'w')
	approxReversalDistance = 0
	for k in range(l):
		if P[k] != k+1 and P[k] != -k-1 :
			idx = [abs(x) for x in P].index(k+1)
			P[k:(idx+1)] = [-1*x for x in reversed(P[k:(idx+1)])]
			approxReversalDistance += 1
			f.write('(' +" ".join(['+'+str(x) if x >0 else str(x) for x in P]) + ')\n')

		if P[k] == -1*(k+1):
			P[k] = -1*P[k]
			approxReversalDistance += 1
			f.write('(' +" ".join(['+'+str(x) if x >0 else str(x) for x in P]) + ')\n')

	return(approxReversalDistance)

def main():
	fin = open(sys.argv[1], 'r')
	#fin.readline()
	dat = fin.readline().rstrip()
	print('read in file\n')
	dat = dat[1:(len(dat)-1)]
	dat = dat.split()
	P = [int(x) for x in dat]
	greedysorting(P)
	print('done\n')


if __name__ == '__main__':
	main()
