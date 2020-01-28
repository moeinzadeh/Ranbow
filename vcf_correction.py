class vcf_correction:
	def __init__(self,hap):
		self.hapFile = hap.outputFolder +hap.prefix + '.single.hap'
		self.modVCFlog = hap.outputFolder + hap.prefix + '.modVCF.log'
		self.vcf = hap.vcf
		self.ploidy = hap.ploidy

		return
	def run (self):
		def trim_tails():
			import sys
			mx_l = -1
			mn_r = sys.maxint
			hap_b = []
			for i in l:
				a = i.split()
				mx_l = mx_l if int(a[2]) < mx_l else int(a[2])
				mn_r = mn_r if int(a[2])  + len(a[3]) > mn_r else int(a[2])  + len(a[3])
				hap_b.append((a[1],int(a[2]), a[3]))


			cnt_block_len = 0
			l_haps = [''] * ploidy
			for j in range(mx_l, mn_r):
				l_snp = []
				for ii in range(ploidy):
					pos = hap_b[ii][1]
					hap = hap_b[ii][2]
					l_snp.append(hap[j - pos])
				
				
				if '-' in l_snp:
					for k in range(ploidy):
						l_haps[k] += '-'
				else:
					for k in range(ploidy):
						l_haps[k] += l_snp[k]
						
			h =  []	
			hname  = []		
			if len(l_haps[0]) != 0 :
				for i in range(ploidy):
					h += [l_haps[i]]
					hname += [hap_b[i][0]]
				return True, h,hname,mx_l
			return False,[],[],-1


		def snp_modification():
			
			for i in range(len(haps[0])):
				mp_snp = {}
				for j in range(ploidy):
					al = haps[j][i]
					if al in mp_snp:
						mp_snp[al] += 1
					else:
						mp_snp[al] = 1

				if  '-' in mp_snp:
					continue
				
				if len(mp_snp) == 1:
					mp[scfName][start+i] = 'D'
				else:
					l_t = []
					cnt_e = 0
					for k in mp_snp:
						l_t.append((mp_snp[k], int(k)))
						cnt_e += 1

					snps = []
					for k in sorted(l_t, key=lambda x:(-x[0],x[1])):
						snps += [(k[1],k[0])]
					mp[scfName][start+i] = ('M',snps)
						
					
		self.deletion = 0
		self.modified = 0
		ploidy = self.ploidy
		
		withPickle = False
		if withPickle:
			import pickle
			if False:
			
				hapFile = self.hapFile
				mp = {}
				for line in open(hapFile):
					if line.startswith(">>>"):
						a = line.split()
						scfName 	= 	a[1]
						scfLen	=	int(a[2])
						mp[scfName]	= {}
						l = []
						cnt = 0
						continue
					cnt +=1
					l+=[line.rstrip()]
					if cnt < ploidy:
						continue
					bo_empty, haps, hapsName, start = trim_tails()
					snp_modification()
					l = []
					cnt = 0

			
				pickle.dump( mp, open( "save.p", "wb" ) )
			else:
			
				mp = pickle.load( open( "save.p", "rb" ) )

		else:
			hapFile = self.hapFile
			mp = {}
			for line in open(hapFile):
				if line.startswith(">>>"):
					a = line.split()
					scfName 	= 	a[1]
					scfLen	=	int(a[2])
					mp[scfName]	= {}
					l = []
					cnt = 0
					continue
				cnt +=1
				l+=[line.rstrip()]
				if cnt < ploidy:
					continue
				bo_empty, haps, hapsName, start = trim_tails()
				snp_modification()
				l = []
				cnt = 0
			



		
		vcfFile = self.vcf
		modVCFFile = self.vcf[:-3] + 'mod.vcf'
		fo = open(modVCFFile,'w')
		fo_mod = open(self.modVCFlog,'w')
		print 'modified VCF file:\t'+modVCFFile
		print "log of modified SNPs:\t"+self.modVCFlog
		
		scf = ''
		fo.write("###Modified vcf file with ranbow haplotyper\n")
		for line in open(vcfFile):
			if line[0] == "#":
				fo.write(line)
				continue
			b = line.split()
			if scf != b[0]:
				cnt = 0
				scf = b[0]
			else:
				if end_of_prev_ref_allele > int(b[1]): # the overlapping SNPs in the code are ignored. the first SNP is maintained and rest of them are removed.
					continue
				cnt += 1
			end_of_prev_ref_allele = 	int(b[1]) + len(b[3])
			b_1 = b
			scf_mp = mp.get(scf, -1)
			if scf_mp == -1:
				fo.write(line)
			else:
				pos_scf_mp = scf_mp.get(cnt,-1)
				if pos_scf_mp == -1:
					fo.write(line)
				else:
					if pos_scf_mp == 'D':
						self.deletion += 1
						fo_mod.write("D\t"+b[0]+'\t'+b[1]+'\n')
						#print "D\t"+b[0]+'\t'+b[1]
						continue
					elif pos_scf_mp[0] == "M":
						self.modified += 1
						new_all = pos_scf_mp[1]
						ans = ''
						ans  = b[0] + '\t' + b[1] + '\t'+ b[2]+'\t'
						b1 = b[4].split(',')
						alleles = [b[3]] + b1
						mp_alleles = {} 
						for i,j in enumerate(alleles):
							mp_alleles[i]  = j
						

						ans += mp_alleles[new_all[0][0]]+'\t'
						GT	= "0/"* new_all[0][1]


						fo_mod.write("M\t"+b[0]+'\t'+b[1]+'\t'+str(pos_scf_mp)+'\n')
						#print "M\t"+b[0]+'\t'+b[1]#+str(pos_scf_mp)
						
						for i in xrange(1,len(new_all)):
							ans += mp_alleles[new_all[i][0]]+','
							GT	+= (str(i)+'/')*new_all[i][1]
						
						ans = ans[:-1] + '\t20\t.\t'
						
						freq_ref = new_all[0][1]
						ans += 'AC=' + str(freq_ref) + ";"
						ans += 'Af=' + str(freq_ref/float(ploidy)) + ";AN=6\tGT\t"+GT[:-1]+"\n"
						fo.write(ans)

		print "Modified SNPs:",self.modified, "\tDeleted SNPs:", self.deletion
