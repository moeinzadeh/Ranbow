import pysam

class IGV:

    def __init__(self, argv,outputFolder, prefix):
        self.hap        = outputFolder+"/"+prefix+'.single.hap'
        self.fas        = outputFolder+"/"+prefix+'.single.hap.igv.fa'
        self.bam        = outputFolder+"/"+prefix+'.single.hap.igv.bam'
        self.sam        = outputFolder+"/"+prefix+'.single.hap.igv.sam'
        self.samtype    = outputFolder+"/"+prefix+'.single.hap.igvtype.sam'
        self.bamSort    = outputFolder+"/"+prefix+'.single.hap.igv.sort'
        import glob, os
        for f in glob.glob(outputFolder+prefix+"*sam.sai"):
            os.remove(f)
    def combined_cigar(self, frag_cigar):
        import re
        l = re.split('(\d+)', frag_cigar)
        i = 2
        ans = ''
        while i < len(l):
            ch = l[i]
            no = int(l[i - 1])
            i += 2
            while i < len(l) and ch == l[i]:
                no += int(l[i - 1])
                i += 2
            if no != 0:
                ans += str(no) + ch
        return ans



    def makeBam_ATCG(self):
        l_scf = []
        for line in open(self.hap):
            a = line.split()
            if a[0] == ">>>":
                l_scf += [{'LN': int(a[2]),'SN':a[1]}]

        fo = open(self.fas,'w')

        header = {'HD': {'VN': '1.0'},
                  'SQ': l_scf}
        outfile = pysam.AlignmentFile(self.bam,'wb', header=header)

        ref_id = -1
        for line in open(self.hap):
            a = line.split()
            if a[0] == ">>>":
                fo.write(">"+a[1]+'\n'+int(a[2])*'A'+'\n')
                ref_id += 1
            else:
                align = pysam.AlignedSegment()
                align.query_name = a[1]
                align.query_sequence = a[3]
                align.reference_start = int(a[2]) + 1
                align.reference_id = ref_id

                prev = '0'
                cnt = 1
                l_cigar = []
                for i in a[3][1:]:
                    if i == '-':
                        if prev == '-':
                            cnt += 1
                        else:
                            l_cigar += [(0, cnt)]
                            cnt = 1
                            prev = '-'
                    else:
                        if prev == '-':
                            l_cigar += [(3,cnt)]
                            cnt = 1
                            prev = '0'
                        else:
                            cnt += 1
                l_cigar += [(0,cnt)]
                align.cigar = tuple(l_cigar)
                outfile.write(align )
        outfile.close()

        pysam.sort('-o', self.bamSort, self.bam)
        pysam.index(self.bamSort)



    def write_to_sam_file(self,l_haps,scf_name):
        for i in sorted(l_haps, key=lambda x: x[0]):
            a = i[1].split('\t')
            sub_name = a[1]
            sub_pos = a[2]
            sub_seq_tmp = a[3]
            sub_cigar_tmp = ""
            sub_seq = ''
            for i in sub_seq_tmp:
                if i != '-':
                    sub_cigar_tmp += '1M'
                    sub_seq += i
                else:
                    sub_cigar_tmp += '1N'
            sub_cigar = self.combined_cigar(sub_cigar_tmp)
            self.fosam.write(sub_name + '\t0\t' + scf_name + '\t' + str(
                int(sub_pos) + 1) + '\t30\t' + sub_cigar + '\t*\t0\t0\t' + sub_seq + '\t*\n')

    def write_to_samtype_file(self,l_type,scf_name):

        for i in sorted(l_type, key=lambda x: x[0]):
            type_poses = i[1]
            cigar = ''
            for inx, v in enumerate(type_poses[:-1]):
                d = type_poses[inx + 1] - v - 1
                if d > 0:
                    cigar += '1M'+str(d )+'N'
                else:
                    cigar += '1M'
            cigar += '1M'
            self.fotype.write(i[3] + '\t0\t' + scf_name + '\t' + str(i[0] + 1) + '\t30\t' + self.combined_cigar(cigar) + '\t*\t0\t0\t' + i[2] + '\t*\n')


    def makeSam(self):
        self.fosam = open(self.sam, 'w')
        self.fotype = open(self.samtype,'w')

        self.fofas = open(self.fas, 'w')

        is_first = True

        for line in open(self.hap, 'r'):
            a = line.split()
            if a[0] == ">>>" and is_first:
                scf_name = a[1]
                scf_len = int(a[2])
                l_haps = []
                l_type = []
                is_first = False
                continue


            if a[0] == ">>>":
                self.write_to_sam_file(l_haps,scf_name)
                self.write_to_samtype_file(l_type,scf_name)
                self.fofas.write(">"+scf_name+'\n'+scf_len*'0'+'\n')
                scf_name = a[1]
                scf_len = int(a[2])
                l_haps = []
                l_type = []
                continue
            type_poses = list(set([int(pos) for pos in a[-2].split('-') if pos!='']))
            type_poses.sort()
            l_haps.append((int(a[2]), line))
            l_type.append((type_poses[0], type_poses, a[-3],a[1]))

        self.write_to_sam_file(l_haps, scf_name)
        self.write_to_samtype_file(l_type, scf_name)
        self.fosam.close()

