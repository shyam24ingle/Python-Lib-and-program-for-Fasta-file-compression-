#Functions lib for fasta file compression
#/////////////////////////////////////////////

#Functions for finding and deleting N sysmbols from seq string 
#It writes N position and length in Npos.txt file and return without N seq string

def find_del_N(s):
	s = s + "_"
	seq = list(s)
	out = open("Npos.txt","a");i = 0
	ncount =seq.count("N");#print ncount
	while "N" in seq:
		  pos = seq.index("N",i)
		  i = pos
		  if seq[pos] == seq[pos+1]:
			l = 1
			while seq[pos] == seq[pos+1]:
			   l = l+1;pos = pos+1
			out.write(((str(i)+'_')+(str(l)+ ' ')))
			del seq[i:pos+1]
		  else:
			out.write((str(i)+' '))
			del seq[i]
	del seq[len(seq)-1];out.write(('*'+' '));out.close()
	return ''.join(seq)
#funcion end$$$

def add_N(seq):
	IN = open("Npos.txt","r") 
	npos = IN.read()
	arr = npos.split('*');lnNN=0
	for ps in arr:
		arr = ps.split(' ')
		for pos in arr:
			#print pos
			if pos != '':
				if '_' in pos:
					p = pos.index('_')
					psNN = int(pos[:p])
					ln = int(pos[p+1:]);psNN +=lnNN;lnNN+=ln
					#print psNN,lnNN
					seq = seq[:psNN]+(('N'*ln)+seq[psNN:]);#print "oooooooooo",len('N'*20)
					#print seq
				else:
					po = int(pos);#print po
					po += lnNN;lnNN+=1
					seq = seq[:po]+('N'+ seq[po:]);#print seq
	return seq
	


#find all repeats

def find_del_repeats(outfname):
	comf = open(outfname,'rb');lp = 0;comff = open(outfname+'1','w');comff.write('');comff.close()
	comff = open(outfname+'1','ab')
	while lp == 0:
		compseq = comf.read(999999)
		if compseq != '':
			compseq+="_"
			seq = list(compseq)
			for i in range(len(seq)):
				if i < len(seq)-1:
					if seq[i] == seq[i+1]:
						pos = i;l = 1
						while seq[pos] == seq[pos+1]:
							l = l + 1;pos = pos+1
						if l > 3:
							del seq[(i+1):(pos+1)];seq.insert(i+1,('_'+str(l)))		
			del seq[len(seq)-1];comff.write(''.join(seq))
		else:lp =1;
	comf.close();comff.close();
	make_gzip(outfname+'1');comlen=len(seq)
	return comlen
	##END

###return repeates

def repeatsIN(fname):
	#un_gzip(fname)
	comf = open(fname,"rb");lp=0
	while lp == 0:
		comseq = comf.read(999999)
		comseq+="$"
		if comseq[0]== '_':comseq =last+comseq		
		comseq= list(comseq);indx=0
		if '_' in comseq:
			while '_' in comseq:
				indx = comseq.index('_',1);i =1;l=''
				while comseq[indx+i].isdigit():l += comseq[indx+i];i+=1;print l
				del comseq[indx:(indx+1)+int(len(l))]
				comseq.insert(indx,(comseq[indx-1]*int(l)))
			del comseq[-1:];last = comseq[-1:]
			if comseq[-1:] == "_" or comseq[0].isdigit():comseq+=comf.read(6)
		else:lp =1
	comf.close()
	
#fnnv= 'compressnew.txt1'
#repeatsIN(fnnv)


 
###Gzip

def make_gzip(filename):
	import gzip
	IN = open(filename,'rb');s = IN.read();IN.close()
	OUT = gzip.GzipFile(filename +'.gz','wb');OUT.write(s);OUT.close()
###End

def un_gzip(filename):
	import gzip
	IN = gzip.open(filename,'rb');s = IN.read();print s;IN.close()
	#OUT = open(filenamereplace(".gz",""),'wb');OUT.write(s);OUT.close()



##largest pattern

def find_largest(text):
    import re
    largest = ''
    i = 1
    while 1:
        m = re.search("(" + ("\w" * i) + ").*\\1.*\\1", text)
        #print m
        if not m:
            break
        largest = m.group(1)
        i += 1
    print largest 
    return largest




####dictionry base compression

def Dbase_comp(seqline):
 import sys
 import string
 genetic_code={ 
 'AAAT':'\x00',	'TAAT':'\x3f',	'GAAT':'\x7d',	'CAAT':'\xbd',\
 'AAAG':'\x01',	'TAAG':'\x40',	'GAAG':'\x7e',	'CAAG':'\xbe',\
 'AAAC':'\x02',	'GAAC':'\x7f',	'CAAC':'\xbf',\
 'AATA':'\x03',	'TATA':'\x42',	'GATA':'\x80',	'CATA':'\xc0',\
 'AATT':'\x04',	'GATT':'\x81',	'CATT':'\xc1',\
 'AATG':'\x05',	'TATG':'\x44',	'GATG':'\x82',	'CATG':'\xc2',\
 'AATC':'\x06',	'TATC':'\x45',	'GATC':'\x83',	'CATC':'\xc3',\
 'AAGA':'\x07',	'TAGA':'\x46',	'GAGA':'\x84',	'CAGA':'\xc4',\
 'AAGT':'\x08',	'GAGT':'\x85',	'CAGT':'\xc5',\
 'AAGG':'\x09',	'TAGG':'\x48',	'GAGG':'\x86',	'CAGG':'\xc6',\
 'AAGC':'\x0a',	'TAGC':'\x49',	'GAGC':'\x87',	'CAGC':'\xc7',\
 'AACA':'\x0b',	'TACA':'\x4a',  'GACA':'\x88',	'CACA':'\xc8',\
 'AACT':'\x0c',  'TACT':'\x4b',	'GACT':'\x89',	'CACT':'\xc9',\
 'AACG':'\x0d',	'TACG':'\x4c',	'GACG':'\x8a',	'CACG':'\xca',\
 'AACC':'\x0e',	'TACC':'\x4d',	'GACC':'\x8b',	'CACC':'\xcb',\
 'ATAA':'\x0f',	'GTAA':'\x8c',	'CTAA':'\xcc',\
 'ATAT':'\x10',	'TTAT':'\x4f',	'GTAT':'\x8d',	'CTAT':'\xcd',\
 'ATAG':'\x11',	'TTAG':'\x50',	'GTAG':'\x8e',	'CTAG':'\xce',\
 'ATAC':'\x12',	'TTAC':'\x51',	'GTAC':'\x8f',	'CTAC':'\xcf',\
 'ATTA':'\x13',	'TTTA':'\x52',	'GTTA':'\x90',	'CTTA':'\xd0',\
 'ATTT':'\x14',	'TTTT':'\x53',	'GTTT':'\x91',	'CTTT':'\xd1',\
 'ATTG':'\x15',	'GTTG':'\x92',	'CTTG':'\xd2',\
 'ATTC':'\x16',	'TTTC':'\x55',	'GTTC':'\x93',	'CTTC':'\xd3',\
 'ATGA':'\x17',	'TTGA':'\x56',	'GTGA':'\x94',	'CTGA':'\xd4',\
 'ATGT':'\x18',	'TTGT':'\x57',	'GTGT':'\x95',	'CTGT':'\xd5',\
 'ATGG':'\x19',	'TTGG':'\x58',	'GTGG':'\x96',	'CTGG':'\xd6',\
 'ATGC':'\x1a',	'TTGC':'\x59',	'GTGC':'\x97',	'CTGC':'\xd7',\
 'ATCA':'\x1b',	'TTCA':'\x5a',	'GTCA':'\x98',	'CTCA':'\xd8',\
 'ATCT':'\x1c',	'TTCT':'\x5b',	'GTCT':'\x99',	'CTCT':'\xd9',\
 'ATCG':'\x1d',	'TTCG':'\x5c',	'GTCG':'\x9a',	'CTCG':'\xda',\
 'ATCC':'\x1e',	'TTCC':'\x5d',	'GTCC':'\x9b',	'CTCC':'\xdb',\
 'AGAA':'\x1f',	'TGAA':'\x5e',	'GGAA':'\x9c',	'CGAA':'\xdc',\
 'AGAT':'\x20',		'GGAT':'\x9d',	'CGAT':'\xdd',\
 'AGAG':'\x21',	'TGAG':'\x60',	'GGAG':'\x9e',	'CGAG':'\xde',\
 'AGAC':'\x22',	'TGAC':'\x61',	'GGAC':'\x9f',	'CGAC':'\xdf',\
 'AGTA':'\x23',	'TGTA':'\x62',	'GGTA':'\xa0',	'CGTA':'\xe0',\
 'AGTT':'\x24',	'TGTT':'\x63',	'GGTT':'\xa1',	'CGTT':'\xe1',\
 'AGTG':'\x25',	'TGTG':'\x64',  'GGTG':'\xa2',	'CGTG':'\xe2',\
 'AGTC':'\x26',	'TGTC':'\x65',	'GGTC':'\xa3',	'CGTC':'\xe3',\
 'AGGA':'\x27',	'TGGA':'\x66',	'GGGA':'\xa4',	'CGGA':'\xe4',\
 'AGGT':'\x28',	'TGGT':'\x67',	'GGGT':'\xa5',	'CGGT':'\xe5',\
 'AGGG':'\x29',	'TGGG':'\x68',	'GGGG':'\xa6',	'CGGG':'\xe6',\
 'AGGC':'\x2a',	'TGGC':'\x69',	'GGGC':'\xa7',	'CGGC':'\xe7',\
 'AGCA':'\x2b',	'TGCA':'\x60',	'GGCA':'\xa8',	'CGCA':'\xe8',\
 'AGCT':'\x2c',	'TGCT':'\x6a',	'GGCT':'\xa9',	'CGCT':'\xe9',\
 'AGCG':'\x2d',	'TGCG':'\x6b',	'GGCG':'\xaa',	'CGCG':'\xea',\
 'AGCC':'\x2e',	'TGCC':'\x6c',	'GGCC':'\xab',	'CGCC':'\xeb',\
 'ACAA':'\x2f',	'TCAA':'\x6d',	'GCAA':'\xac',	'CCAA':'\xec',\
 'ACAT':'\x30',	'TCAT':'\x6e',	'GCAT':'\xad',	'CCAT':'\xed',\
 'ACAG':'\x31',	'TCAG':'\x6f',	'GCAG':'\xae',	'CCAG':'\xee',\
 'ACAC':'\x32',	'TCAC':'\x70',	'GCAC':'\xaf',	'CCAC':'\xef',\
 'ACTA':'\x33',	'TCTA':'\x71',	'GCTA':'\xb0',	'CCTA':'\xf0',\
 'ACTT':'\x34',	'TCTT':'\x72',	'GCTT':'\xb1',	'CCTT':'\xf1',\
 'ACTG':'\x35',	'TCTG':'\x73',	'GCTG':'\xb2',	'CCTG':'\xf2',\
 'ACTC':'\x36',	'TCTC':'\x74',	'GCTC':'\xb3',	'CCTC':'\xf3',\
 'ACGA':'\x37',	'TCGA':'\x75',	'GCGA':'\xb4',	'CCGA':'\xf4',\
 'ACGT':'\x38',	'TCGT':'\x76',	'GCGT':'\xb5',	'CCGT':'\xf5',\
 'ACGG':'\x39', 'TCGG':'\x77',	'GCGG':'\xb6',	'CCGG':'\xf6',\
 'ACGC':'\x3a',	'TCGC':'\x78',	'GCGC':'\xb7',	'CCGC':'\xf7',\
 'ACCA':'\x3b',	'TCCA':'\x79',	'GCCA':'\xb8',	'CCCA':'\xf8',\
 'ACCT':'\x3c',	'TCCT':'\x7a',	'GCCT':'\xb9',	'CCCT':'\xf9',\
 'ACCG':'\x3d',	'TCCG':'\x7b',	'GCCG':'\xbb',	'CCCG':'\xfa',\
 'ACCC':'\x3e',	'TCCC':'\x7c',  'GCCC':'\xbc',	'CCCC':'\xfb',\
 'AAAA':'\xfc', 'TAAA':'\xfd', 'GAAA':'\xfe', 'CAAA':'\xff'}
 compseq = ""
 outfname = 'compressnew.txt'
 fil = open("compressnew.txt", "ab")
 for i in range(0,len(seqline),4):
    codon = seqline[i:i+4]
    if (genetic_code.has_key(codon)):
	 compseq = genetic_code[codon]
	 fil.write(compseq)
    else:fil.write(codon)
 fil.close()
 return outfname

###

def decomp(outfname):
	genetic_code = {'\x00'	:	'AAAT'	,	'\x3f'	:	'TAAT'	,	'\x7d'	:	'GAAT'	,	'\xbd'	:	'CAAT'	,	
	'\x01'	:	'AAAG'	,	'\x40'	:	'TAAG'	,	'\x7e'	:	'GAAG'	,	'\xbe'	:	'CAAG'	,	
	'\x02'	:	'AAAC'	,	'\x7f'	:	'GAAC'	,	'\xbf'	:	'CAAC'	,	'\xc0'	:	'CATA'	,	
	'\x03'	:	'AATA'	,	'\x42'	:	'TATA'	,	'\x80'	:	'GATA'	,	'\xc2'	:	'CATG'	,	
	'\x04'	:	'AATT'	,	'\x81'	:	'GATT'	,	'\xc1'	:	'CATT'	,	'\xc3'	:	'CATC'	,	
	'\x05'	:	'AATG'	,	'\x44'	:	'TATG'	,	'\x82'	:	'GATG'	,	'\xc4'	:	'CAGA'	,
	'\x06'	:	'AATC'	,	'\x45'	:	'TATC'	,	'\x83'	:	'GATC'	,	'\xc6'	:	'CAGG'	,	
	'\x07'	:	'AAGA'	,	'\x46'	:	'TAGA'	,	'\x84'	:	'GAGA'	,	'\xc7'	:	'CAGC'	,	
	'\x08'	:	'AAGT'	,	'\x85'	:	'GAGT'	,	'\xc5'	:	'CAGT'	,	'\xc8'	:	'CACA'	,	
	'\x09'	:	'AAGG'	,	'\x48'	:	'TAGG'	,	'\x86'	:	'GAGG'	,	'\xc9'	:	'CACT'	,	
	'\x0a'	:	'AAGC'	,	'\x49'	:	'TAGC'	,	'\x87'	:	'GAGC'	,	'\xca'	:	'CACG'	,
	'\x0b'	:	'AACA'	,	'\x4a'	:	'TACA'	,	'\x88'	:	'GACA'	,	'\xcb'	:	'CACC'	,	
	'\x0c'	:	'AACT'	,	'\x4b'	:	'TACT'	,	'\x89'	:	'GACT'	,	'\xcd'	:	'CTAT'	,	
	'\x0d'	:	'AACG'	,	'\x4c'	:	'TACG'	,	'\x8a'	:	'GACG'	,	'\xce'	:	'CTAG'	,	
	'\x0e'	:	'AACC'	,	'\x4d'	:	'TACC'	,	'\x8b'	:	'GACC'	,	'\xcf'	:	'CTAC'	,	
	'\x0f'	:	'ATAA'	,	'\x8c'	:	'GTAA'	,	'\xcc'	:	'CTAA'	,	'\xd0'	:	'CTTA'	,	
	'\x10'	:	'ATAT'	,	'\x4f'	:	'TTAT'	,	'\x8d'	:	'GTAT'	,	'\xd1'	:	'CTTT'	,	
	'\x11'	:	'ATAG'	,	'\x50'	:	'TTAG'	,	'\x8e'	:	'GTAG'	,	'\xd3'	:	'CTTC'	,	
	'\x12'	:	'ATAC'	,	'\x51'	:	'TTAC'	,	'\x8f'	:	'GTAC'	,	'\xd4'	:	'CTGA'	,	
	'\x13'	:	'ATTA'	,	'\x52'	:	'TTTA'	,	'\x90'	:	'GTTA'	,	'\xd5'	:	'CTGT'	,	
	'\x14'	:	'ATTT'	,	'\x53'	:	'TTTT'	,	'\x91'	:	'GTTT'	,	'\xd6'	:	'CTGG'	,	
	'\x15'	:	'ATTG'	,	'\x92'	:	'GTTG'	,	'\xd2'	:	'CTTG'	,	'\xd7'	:	'CTGC'	,	
	'\x16'	:	'ATTC'	,	'\x55'	:	'TTTC'	,	'\x93'	:	'GTTC'	,	'\xd8'	:	'CTCA'	,	
	'\x17'	:	'ATGA'	,	'\x56'	:	'TTGA'	,	'\x94'	:	'GTGA'	,	'\xd9'	:	'CTCT'	,	
	'\x18'	:	'ATGT'	,	'\x57'	:	'TTGT'	,	'\x95'	:	'GTGT'	,	'\xda'	:	'CTCG'	,	
	'\x19'	:	'ATGG'	,	'\x58'	:	'TTGG'	,	'\x96'	:	'GTGG'	,	'\xdb'	:	'CTCC'	,	
	'\x1a'	:	'ATGC'	,	'\x59'	:	'TTGC'	,	'\x97'	:	'GTGC'	,	'\xdc'	:	'CGAA'	,	
	'\x1b'	:	'ATCA'	,	'\x5a'	:	'TTCA'	,	'\x98'	:	'GTCA'	,	'\xdd'	:	'CGAT'	,	
	'\x1c'	:	'ATCT'	,	'\x5b'	:	'TTCT'	,	'\x99'	:	'GTCT'	,	'\xde'	:	'CGAG'	,	
	'\x1d'	:	'ATCG'	,	'\x5c'	:	'TTCG'	,	'\x9a'	:	'GTCG'	,	'\xdf'	:	'CGAC'	,	
	'\x1e'	:	'ATCC'	,	'\x5d'	:	'TTCC'	,	'\x9b'	:	'GTCC'	,	'\xe0'	:	'CGTA'	,	
	'\x1f'	:	'AGAA'	,	'\x5e'	:	'TGAA'	,	'\x9c'	:	'GGAA'	,	'\xe1'	:	'CGTT'	,	
	'\x20'	:	'AGAT'	,		'\x9d'	:	'GGAT'	,	'\xe2'	:	'CGTG'	,	
	'\x21'	:	'AGAG'	,	'\x60'	:	'TGAG'	,	'\x9e'	:	'GGAG'	,	'\xe3'	:	'CGTC'	,	
	'\x22'	:	'AGAC'	,	'\x61'	:	'TGAC'	,	'\x9f'	:	'GGAC'	,	'\xe4'	:	'CGGA'	,	
	'\x23'	:	'AGTA'	,	'\x62'	:	'TGTA'	,	'\xa0'	:	'GGTA'	,	'\xe5'	:	'CGGT'	,	
	'\x24'	:	'AGTT'	,	'\x63'	:	'TGTT'	,	'\xa1'	:	'GGTT'	,	'\xe6'	:	'CGGG'	,	
	'\x25'	:	'AGTG'	,	'\x64'	:	'TGTG'	,	'\xa2'	:	'GGTG'	,	'\xe7'	:	'CGGC'	,	
	'\x26'	:	'AGTC'	,	'\x65'	:	'TGTC'	,	'\xa3'	:	'GGTC'	,	'\xe8'	:	'CGCA'	,	
	'\x27'	:	'AGGA'	,	'\x66'	:	'TGGA'	,	'\xa4'	:	'GGGA'	,	'\xe9'	:	'CGCT'	,	
	'\x28'	:	'AGGT'	,	'\x67'	:	'TGGT'	,	'\xa5'	:	'GGGT'	,	'\xea'	:	'CGCG'	,	
	'\x29'	:	'AGGG'	,	'\x68'	:	'TGGG'	,	'\xa6'	:	'GGGG'	,	'\xeb'	:	'CGCC'	,	
	'\x2a'	:	'AGGC'	,	'\x69'	:	'TGGC'	,	'\xa7'	:	'GGGC'	,	'\xec'	:	'CCAA'	,	
	'\x2b'	:	'AGCA'	,	'\x60'	:	'TGCA'	,	'\xa8'	:	'GGCA'	,	'\xed'	:	'CCAT'	,	
	'\x2c'	:	'AGCT'	,	'\x6a'	:	'TGCT'	,	'\xa9'	:	'GGCT'	,	'\xee'	:	'CCAG'	,	
	'\x2d'	:	'AGCG'	,	'\x6b'	:	'TGCG'	,	'\xaa'	:	'GGCG'	,	'\xef'	:	'CCAC'	,	
	'\x2e'	:	'AGCC'	,	'\x6c'	:	'TGCC'	,	'\xab'	:	'GGCC'	,	'\xf0'	:	'CCTA'	,	
	'\x2f'	:	'ACAA'	,	'\x6d'	:	'TCAA'	,	'\xac'	:	'GCAA'	,	'\xf1'	:	'CCTT'	,	
	'\x30'	:	'ACAT'	,	'\x6e'	:	'TCAT'	,	'\xad'	:	'GCAT'	,	'\xf2'	:	'CCTG'	,	
	'\x31'	:	'ACAG'	,	'\x6f'	:	'TCAG'	,	'\xae'	:	'GCAG'	,	'\xf3'	:	'CCTC'	,	
	'\x32'	:	'ACAC'	,	'\x70'	:	'TCAC'	,	'\xaf'	:	'GCAC'	,	'\xf4'	:	'CCGA'	,	
	'\x33'	:	'ACTA'	,	'\x71'	:	'TCTA'	,	'\xb0'	:	'GCTA'	,	'\xf5'	:	'CCGT'	,	
	'\x34'	:	'ACTT'	,	'\x72'	:	'TCTT'	,	'\xb1'	:	'GCTT'	,	'\xf6'	:	'CCGG'	,	
	'\x35'	:	'ACTG'	,	'\x73'	:	'TCTG'	,	'\xb2'	:	'GCTG'	,	'\xf7'	:	'CCGC'	,	
	'\x36'	:	'ACTC'	,	'\x74'	:	'TCTC'	,	'\xb3'	:	'GCTC'	,	'\xf8'	:	'CCCA'	,	
	'\x37'	:	'ACGA'	,	'\x75'	:	'TCGA'	,	'\xb4'	:	'GCGA'	,	'\xf9'	:	'CCCT'	,	
	'\x38'	:	'ACGT'	,	'\x76'	:	'TCGT'	,	'\xb5'	:	'GCGT'	,	'\xfa'	:	'CCCG'	,	
	'\x39'	:	'ACGG'	,	'\x77'	:	'TCGG'	,	'\xb6'	:	'GCGG'	,	'\xfb'	:	'CCCC'	,	
	'\x3a'	:	'ACGC'	,	'\x78'	:	'TCGC'	,	'\xb7'	:	'GCGC'	,	'\xff'	:	 'CAAA'	,	
	'\x3b'	:	'ACCA'	,	'\x79'	:	'TCCA'	,	'\xb8'	:	'GCCA'	,					
	'\x3c'	:	'ACCT'	,	'\x7a'	:	'TCCT'	,	'\xb9'	:	'GCCT'	,				
	'\x3d'	:	'ACCG'	,	'\x7b'	:	'TCCG'	,	'\xbb'	:	'GCCG'	,					
	'\x3e'	:	'ACCC'	,	'\x7c'	:	'TCCC'	,	'\xbc'	:	'GCCC'	,					
	'\xfc'	:	'AAAA'	,	'\xfd'	:	'TAAA'	,	'\xfe'	:	'GAAA'}					
	comf = open(outfname,'rb')
	fil = open("orginal_seq.txt","a");lp = 0
	while lp == 0:
		compseq = comf.read(999999)
		if compseq != '':
			for s in compseq:
				if (genetic_code.has_key(s)):seq = genetic_code[s];fil.write(seq)
				else:fil.write(s)
		else:lp = 1
	fil.close()




#empty all out file

def empty_outf():
  fil = open("orginal_seq.txt","w");fil.write('');fil.close()
  comf = open("compressnew.txt", "w");comf.write('');comf.close() 
  allrepf = open("allrepeats.txt","w");allrepf.write('');allrepf.close()
  Npos = open("Npos.txt","w");Npos.write('');Npos.close()
  print "***all output file are empty***\n compressing seq....."
