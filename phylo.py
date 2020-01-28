
class phylo:
	
		
	def __init__(self,argv):

		import os
		import sys, getopt
		mp_params = {}
		
		def read_parameters():
			for i in range(2,len(sys.argv),2):
				mp_params[sys.argv[i]] = sys.argv[i+1]

			if ('-par' not in mp_params ):
				print "Error: indicate -par argument (parameter file)"
				exit()

			fi = open(mp_params['-par'], 'r')
			line = fi.readline()
			while line != '':
				a = line.split()
				if a[0] not in mp_params.keys() and a[0][0] == '-':
					if len(a) > 2:
						mp_params[a[0]] = (a[1],a[2])
					else:
						mp_params[a[0]] = a[1]
				line = fi.readline()

			if "-mode" not in mp_params:
				print "please indicate -mode parameters"
				exit()
			else:
				self.act = mp_params["-mode"]
		
		read_parameters()

		self.prefix		     						= 		mp_params.get('-outputPrefixName'	,"ranbow")
		self.outputFolder 							= 		mp_params['-outputFolderBase']
		self.phyloFolder								= 		self.outputFolder+'phylo/'
		self.ploidy									=		int(mp_params['-ploidy'])
		self.hapFile 									= 		self.outputFolder + self.prefix + '.single.hap'
		#self.hapFile 									= 		'/scratch/cluster/moeinzadeh/nt_plant_review/hap/FINAL_Single/ranbow.single.hap'
		self.c1File									=		self.hapFile+".c1"


		if self.act == "run":
			self.vcf 									= 		mp_params['-vcfFile']
			self.vcf_inx									= 		self.vcf+'.index'
			self.noCpu									=		int(mp_params['-noProcessor'])
			self.cpuInx									=		int(mp_params['-processorIndex'])
			self.ref 									= 		mp_params['-refFile']
			self.ref_fai									= 		self.ref + '.fai'
			if not os.path.exists(self.ref_fai):
				print 'No index file: ',self.ref + '.fai'
				print 'Run the "ranbow.py index" first'
				exit()
			
		if self.act == "collect":
			self.noCpu									=		int(mp_params['-noProcessor'])
			self.mnHapLen									=		int(mp_params.get('-minHapLen',3))

	
	def load_inx_file(self,f):
		mp = {}
		for i in open(f):
			a = i.split()
			mp[a[0]] = int (a[1])
		return mp


	def make_global_alignment_based_on_cigar(self, l_alleles, l_cigar):
		l_cigar = [str(len(l_alleles[0]))+'M'] + l_cigar
		l_cigar_extended = []
		for i in l_cigar:
			import re
			k = re.split('(\d+)',i)	
			tmp = ''
			for t in range(1 , len(k) , 2):
				for tt in range(int(k[t])):
					tmp += k[t+1]
			l_cigar_extended.append(tmp)
		#print l_cigar_extended
		pointer_cigar = [0] * len (l_cigar_extended)
		pointer_allele = [0] * len (l_cigar_extended)
		new_seq = [''] * len (l_cigar_extended)
		no_allele = len(l_cigar)
		
		
		bo_continiue = True
		
		
		bo_test = False
		
		while bo_continiue:
			no_d = 0
			no_i = 0
			no_m = 0
			l_d = []
			l_i = []
			for j in range(no_allele):
				i = l_cigar_extended[j][pointer_cigar[j]]
				if i == 'M' or i == 'X':
					no_m += 1
				if i == 'I':
					no_i += 1
					l_i.append(j)
				if i == 'D':
					no_d += 1
					l_d.append(j)
					
			if no_i > 0:
				bo_test = True
				for i in range(no_allele):
					if i not in l_i:
						new_seq[i] += '-'
					else:
						new_seq[i] += l_alleles[i][pointer_allele[i]]
						pointer_cigar[i] += 1
						pointer_allele[i] += 1

			elif no_d > 0:
				bo_test = True
				for i in range(no_allele):
					if i in l_d:
						new_seq[i] += '-'
						pointer_cigar[i] += 1
					else:
						if l_cigar_extended[i][pointer_cigar[i]] == 'M' or l_cigar_extended[i][pointer_cigar[i]] == 'X':
							new_seq[i] += l_alleles[i][pointer_allele[i]]
							pointer_cigar[i] += 1
							pointer_allele[i] += 1
						
						
			elif no_allele == no_m:
				for i in range(no_allele):
					new_seq[i] += l_alleles[i][pointer_allele[i]]
					pointer_cigar[i] += 1
					pointer_allele[i] += 1

			bo_continiue = True
			for i in range(no_allele):
				if pointer_cigar[i] >= len(l_cigar_extended[i]):
					bo_continiue = False
		
		return new_seq, l_cigar_extended
		if bo_test:
			
			
			print '---------------------------------'
			for i in new_seq		:
				print i
			print 
			print l_cigar_extended
			print new_seq		
			print l_alleles
			cnt_indel += 1
			print '---------------------------------'
			print 
			


		
	def indexingVCF_aligned(self, scName, fi):
		ploidy_no = self.ploidy
		arrInx2pos = []
		arrInxVariation = []
		arrVarPattern = []
		
		line = fi.readline()
		a = line.split('\t')
		
		while a[0] == scName:
			#print a
			allele_pos = int(a[1])
			end_ref_allele = allele_pos + len(a[3])
			l_allele = [a[3]]+ a[4].split(',')
			#print l_allele
			#exit()
			
			
			tmp1 = a[7].split(';')
			l_cigar = tmp1[6][6:].split(',')
			
			l_allele_new, l_cigar_new = self.make_global_alignment_based_on_cigar(l_allele,l_cigar)
			
			
			arrInxVariation.append(l_allele_new[:ploidy_no])
			arrVarPattern.append(l_cigar_new[:ploidy_no])
			arrInx2pos.append(allele_pos)
				
				
			line = fi.readline()
			a = line.split('\t')
			# for those consecutive SNPs that one covers the other
			while a[0] == scName and end_ref_allele > int(a[1]):
				line = fi.readline()
				a = line.split('\t')
			
		return (arrInx2pos,arrInxVariation,arrVarPattern)




	def extract_dna_of_a_fragment(self,hap,hap_str,hap_name,scf,vcf_info):
		l_allele__index_2_realPos 	= 	vcf_info[0]
		l_index_2_alleles			=	vcf_info[1]
		#l_index_2_cigar			=	vcf_info[2]
		name 					=	hap_name
		frag_start	 			=	hap_str
		frag						=	hap
		frag_end				= 	frag_start + len(hap) - 1
		interval_no_tail_both_side = ( l_allele__index_2_realPos[frag_start] , l_allele__index_2_realPos[frag_end] + len(l_index_2_alleles[frag_end][0]) )

		ref_allele = 0
		if frag_start == 0:
			frag_nc = scf[:l_allele__index_2_realPos[0]] 
			frag_nc_start_pos = 0
			#frag_cigar = str(l_allele__index_2_realPos[0])+'M'
		else:
			frag_nc = scf[l_allele__index_2_realPos[frag_start - 1] + len (l_index_2_alleles[frag_start - 1 ][ref_allele]) :l_allele__index_2_realPos[frag_start]]
			frag_nc_start_pos = l_allele__index_2_realPos[frag_start - 1] + len (l_index_2_alleles[frag_start - 1 ][0])
			#frag_cigar = str( l_allele__index_2_realPos[frag_start] - l_allele__index_2_realPos[frag_start - 1] - len (l_index_2_alleles[frag_start - 1 ][ref_allele]))+'M' 
		
		shared_region_lenght = len(frag_nc)
		
		end_pos_on_noneallelic_region = start_pos_on_noneallelic_region = -1
		for i in range(len(frag)): 
			if frag[i] == '0':
				allele = l_index_2_alleles[frag_start + i][int(frag[i])]
				frag_nc += allele
				#frag_cigar += str(len(allele))+	'M'
			
			elif frag[i] == '-':
				alleleSize = len (l_index_2_alleles[frag_start + i][ref_allele])
				frag_nc += alleleSize * '?'
				#frag_cigar += str(alleleSize)+	'M'
				
			else:
				#frag_cigar += l_index_2_cigar[frag_start + i][int(frag[i]) - 1]
				frag_nc += l_index_2_alleles[frag_start + i][int(frag[i])]



			start_pos_on_noneallelic_region = l_allele__index_2_realPos[frag_start + i] + len (l_index_2_alleles[frag_start + i][ref_allele])
			
			if frag_start + i + 1 < len (l_index_2_alleles):
				end_pos_on_noneallelic_region = l_allele__index_2_realPos[frag_start + i + 1]
			elif frag_start + i + 1 == len (l_index_2_alleles):
				end_pos_on_noneallelic_region = len(scf)
			else:
				print 'error: fragment passed the boundary of scaffold'

			frag_nc += scf[start_pos_on_noneallelic_region:end_pos_on_noneallelic_region]
			shared_region_lenght += len(scf[start_pos_on_noneallelic_region:end_pos_on_noneallelic_region])
			#frag_cigar += str(end_pos_on_noneallelic_region - start_pos_on_noneallelic_region) + 'M'
		return frag_nc_start_pos, frag_nc, shared_region_lenght , interval_no_tail_both_side #, self.combined_cigar(frag_cigar)


	def trim_tails_new(self,l):
		import sys
		mx_l = -1
		mn_r = mn_l = sys.maxint
		hap_b = []
		for i in l:
			st = i[1]
			en = i[1] + len(i[2])
			mx_l = mx_l if st < mx_l else st
			mn_r = mn_r if en > mn_r else en
			mn_l = mn_l if st > mn_l else st
		if mx_l >= mn_r:
			return True,[], -1
		new_l = []
		for i in l:
			st = i[1]
			en = i[1] + len(i[2])
			
			if en - mn_r == 0:
				new_l += [ (i[2][ mx_l - st: ] , i[0] )]
			else:
				new_l += [ (i[2][ mx_l - st: -(en - mn_r) ] , i[0])]
		return False, new_l, mx_l


		
	def pre_run (self):
		import numpy as np
		import os
		import subprocess
		import pickle

		ploidy = self.ploidy
		hapFile = self.hapFile
		l_bk = []
		tmp = 0
		for line in open(hapFile):
			a = line.split()
			if line.startswith(">>>"):
				scfName 	= 	a[1]
				scfLen	=	int(a[2])
				l = []
				cnt = 0
				continue
			cnt +=1
			l+=[(a[1], int(a[2]), a[3] )]
			
			if cnt < ploidy:
				continue
			bo_empty, haps, start = self.trim_tails_new(l)
			if not bo_empty:
				l_bk += [(haps, start,scfName)]
			cnt = 0
			l = []
			
			#tmp += 1
			#if tmp %1000 == 0:
			#	break
			#	pass
		
		
		l_bk.sort(key=lambda x: x[2])
		pickle.dump( l_bk, open( "c1.p", "wb" ) )
		#return l_bk

	def make_pickle (self):
		import numpy as np
		import os
		import subprocess
		import pickle

		ploidy = self.ploidy
		hapFile = self.hapFile
		l_bk = []
		tmp = 0
		for line in open(hapFile):
			a = line.split()
			if line.startswith(">>>"):
				scfName 	= 	a[1]
				scfLen	=	int(a[2])
				l = []
				cnt = 0
				continue
			cnt +=1
			l+=[(a[1], int(a[2]), a[3] )]
			
			if cnt < ploidy:
				continue
			bo_empty, haps, start = self.trim_tails_new(l)
			if not bo_empty:
				l_bk += [(haps, start,scfName)]
			cnt = 0
			l = []
			
			tmp += 1
			if tmp %1000 == 0:
				break
				pass
		
		
		l_bk.sort(key=lambda x: x[2])
		pickle.dump( l_bk, open( "c1.p", "wb" ) )
		#pickle.dump( l_bk, open( "c1.p_tmp", "wb" ) )


			
	def find_pairs  (self, nwk, mp_name_hap, vcf_info, shared_region_len,hap_st):
		if nwk.find(".") == -1 :
			print 'empty'
			return False, '', []
			
		s = nwk
		d = [',',';',':','(',')']
		for i in d:
			s = s.replace(i,' ')
		s = s.split()	
		top = ''
		for c in nwk:
			if c in ("(",")"):
				top += c
		
		if top == "((()())())" or top == "(()(()()))":
			ln_hap = len(mp_name_hap[mp_name_hap.keys()[0]])
			l_sub = []
			for i in xrange(hap_st,hap_st + ln_hap):
				sw = True
				for j in vcf_info[1][i]:
					if "-" in j:
						sw = False
						break
				if sw :
					l_sub += [i - hap_st]
			if len(l_sub) == 0:
				return False, '', []

			cnt = 0
			l = []
			l_mutation_rate = []
			for i in s:
				if not ('0' <= i[0] <= '9'):
					l += [i]
					cnt += 1
					if cnt % 2 == 0:
						hapName1 = l[0]
						hapName2 = l[1]
						mutation = 0
						for i in xrange(ln_hap):
							if mp_name_hap[hapName1][i] != mp_name_hap[hapName2][i] and i in l_sub:
								mutation += 1
						l_mutation_rate += [mutation / float(shared_region_len)]		
						l =	[]
			return True, top, l_mutation_rate
		return True, top, []
				
			
		
		
	def run (self):
		import numpy as np
		import os
		import subprocess
		import pickle

		ploidy = self.ploidy
		hapFile = self.hapFile
		l_bk = []
		noCpu = self.noCpu
		#noCpu = 100
		cpuInx = self.cpuInx
		
		l_bk = pickle.load( open( "c1.p", "rb" ) )
			
		no_bk = len(l_bk)
		
		
		mp_vcf 		= 	self.load_inx_file (self.vcf_inx)
		fo_vcf 		= 	open(self.vcf)

		import pysam
		fastafile 	= 	pysam.FastaFile(self.ref)
		

		l_1cpu = np.array_split(range(no_bk),noCpu)
		scfPrev = ""
		
		if not os.path.exists(self.phyloFolder): 
			os.makedirs(self.phyloFolder)
			
		if not os.path.exists(self.phyloFolder+'log'): 
			os.makedirs(self.phyloFolder+ 'log')

		
		file_hap_len 	= self.phyloFolder+'haplen.'+str(cpuInx)+'.txt'
		fo_hap_len 		= open(file_hap_len,'w')
		
		l_this_cpu = l_bk[l_1cpu[cpuInx][0]:l_1cpu[cpuInx][-1]+ 1]
		del l_bk
		for i in l_this_cpu:
			
			file_meg_name_in 	= self.phyloFolder+'tmp'+str(cpuInx)+'.meg'
			file_meg_name_out	= self.phyloFolder+'tmp'+str(cpuInx)+'.nwk'	
			fi_mega = open(file_meg_name_in, 'w')
			fi_mega.write ( "#mega\n!Title ;\n!Format DataType=DNA indel=-;\n\n" )
			
			haps 		= 	i[0]
			hap_st		=	i[1]
			scf 			= 	i[2]
			refSeq = fastafile.fetch(scf)
			
			
			if scfPrev != scf:
				fo_vcf.seek(mp_vcf[scf])
				vcf_info = self.indexingVCF_aligned(scf,fo_vcf)
				scfPrev = scf
			
			mp_name_hap = {}
			start_pos_dna = -1 
			end_pos_dna = -1 
			extract_dna_of_a_fragment = (-1, -1)
			
			for hap,hap_name in haps:
				start_pos_dna, seq, shared_region_len , extract_dna_of_a_fragment = self.extract_dna_of_a_fragment(hap,hap_st,hap_name,refSeq,vcf_info)
				
				mp_name_hap [hap_name] = hap
				fi_mega.write("#" + hap_name + '\n' + seq + '\n\n')
			fi_mega.close()

			a = subprocess.check_output(['bash','./phylo.sh',file_meg_name_in, file_meg_name_out])
			
			#print a+ 'asdfasdfasdfasd'

			for nwk in open(file_meg_name_out):
				is_correct, top, mutRate = self.find_pairs(nwk,mp_name_hap,vcf_info ,shared_region_len,hap_st)
				#print is_correct, top, mutRate
				if is_correct:
					fo_hap_len.write ( str(len(haps[0][0])) + '\t')
					fo_hap_len.write ( top + '\t' + nwk.rstrip() + '\t' + str(extract_dna_of_a_fragment[0]) + '\t' + str(extract_dna_of_a_fragment[1]) + '\t' +str(mutRate) +'\n')
				break # just the first line is used
			

	def tree_char(self, infile):
		fi = open(infile, 'r')
		line = fi.readline()
		cnt = 0
		l_all =	[0.0] * 10
		l_len_all_10 = [[],[],[],[],[],[],[],[],[],[]]
		l_a1 = []
		l_a2_12 = []
		while line!= '':
			a = line.split()
			if "(()(()()))" == a[1]:
				cnt += 1
				b = a[2].split(':')
				l = []
				for i in range(1, 11):
					s= ''
					for j in b[i]:
						if '0' <= j <='9' or j == '.':
							s += j
						else:
							break
					l.append(float(s))
				l_4_2 = l[3:]+l[:3]
				for i in range(len(l)):
					l_all[i] += l_4_2[i]
					l_len_all_10[i].append(l_4_2[i])
				
				a1 = float(a[5][1:-1])
				a2_1 = float(a[6][	:-1])
				a2_2 = float(a[7][1:-1])
				l_a1	   += [a1]
				l_a2_12 += [a2_1, a2_2]
				
			if "((()())())" == a[1]:
				cnt += 1
				b = a[2].split(':')
				l = []
				for i in range(1, 11):
					s= ''
					for j in b[i]:
						if '0' <= j <='9' or j == '.':
							s += j
						else:
							break
					l.append(float(s))
				for i in range(len(l)):
					l_all[i] += l[i]
					l_len_all_10[i].append(l[i])
				a2_1 = float(a[5][1:-1])
				a2_2 = float(a[6][:-1])
				a1 = float(a[7][1:-1])

				l_a1	   += [a1]
				l_a2_12 += [a2_1, a2_2]
		

			line = fi.readline()

		
		print "Tree topology:		(((A-B)I (C-D)J ) M (E-F) N)"
		print "Branch Name\tMean relative distance\tSD relative distance"
		edge_name = "ABICDJMEFN"
		if cnt == 0:
			'no such a (((A-B)I (C-D)J ) M (E-F) N) topology'
			exit()
		l_avr = [0.0]*10
		for i in range(len(l_all)):
			l_avr[i] = l_all[i] / float(cnt)
	
		import numpy
				
		for i in range(10):
			print edge_name[i]+':\t\t', numpy.mean(l_len_all_10[i]),'\t\t', numpy.std(l_len_all_10[i]) 
		print 
		print 'Mutation rate for A-B and C-D branches:' , numpy.mean(l_a2_12)
		print 'Mutation rate for E-F branch:' , numpy.mean(l_a1)
