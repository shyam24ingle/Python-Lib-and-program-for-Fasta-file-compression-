import sys
import string
import timeit
import fastacomp 

start = timeit.default_timer()

fil = open(sys.argv[1])
hinfo = fil.readline()
fastacomp.empty_outf()
lp = 0
while lp == 0:
	s = fil.read(999999)
	if s != "":
	  s = s.replace("\n","")
	  seqline = fastacomp.find_del_N(s)
	  seqline = fastacomp.find_del_repeats(seqline)
	  seqline = fastacomp.Dbase_comp(seqline)
	else:
	  lp = 1


stop = timeit.default_timer()
print  stop - start

fil.close()
