
import pysam
import datetime
import os
import sys
import pysam as ps
import re
import ast
from collections import Counter
from bisect import bisect_left
from bisect import bisect_right


class haplotyper:
    class _outfiles:
        def __init__(self, mp_params):

            basefolder = mp_params['-outputFolderBase']
            folder = mp_params['-processorIndex']
            prefix = mp_params['-outputPrefixName']
            if not os.path.exists(basefolder):
                os.makedirs(basefolder)

            outputFolder = basefolder + '/' + folder
            if not os.path.exists(outputFolder):
                os.makedirs(outputFolder)

            self.outputFilePath = outputFolder + '/' + prefix

            self.file_hap_single = self.outputFilePath + '.single.hap'
            self.fo_hap_single = open(self.file_hap_single, 'w')

            # self.file_hap_pair = outputFilePath + '.pair.hap'
            # self.fo_hap_pair = open(self.file_hap_pair, 'w')

            self.file_hap_single_sam = self.outputFilePath + '.single.hap.sam'
            self.fo_hap_single_sam = open(self.file_hap_single_sam, 'w')

            # self.file_hap_pair_sam = outputFilePath + '.pair.hap.sam'
            # self.fo_hap_pair_sam = open(self.file_hap_pair_sam, 'w')

            self.file_hap_single_fasta = self.outputFilePath + '.single.hap.fasta'
            self.fo_hap_single_fasta = open(self.file_hap_single_fasta, 'w')

            # self.file_hap_pair_fasta = outputFilePath + '.pair.hap.fasta'
            # self.fo_hap_pair_fasta = open(self.file_hap_pair_fasta, 'w')

            # self.file_tmp_PE_reads = outputFilePath + '.PE_frag'
            # self.fo_file_tmp_PE_reads = open(self.file_tmp_PE_reads, 'w')
            self.folder_index = folder



    class _outfiles_eval:
        def __init__(self, mp_params):

            basefolder = mp_params['-outputFolderBase'] + 'eval/'
            folder = mp_params['-processorIndex']
            prefix = mp_params['-outputPrefixName']
            if not os.path.exists(basefolder):
                os.makedirs(basefolder)

            outputFolder = basefolder + '/' + folder
            if not os.path.exists(outputFolder):
                os.makedirs(outputFolder)

            outputFilePath = outputFolder + '/' + prefix

            fileEval = outputFilePath + '.single.eval'
            self.fo_eval_singel = open(fileEval, 'w')

            # fileEval = outputFilePath + '.pair.eval'
            # self.fo_eval_pair = open(fileEval, 'w')

            # fileEvalNoBlock = outputFilePath + '.pair_noBlock.eval'
            # self.fo_eval_pair_NoBlock = open(fileEvalNoBlock, 'w')

            self.outputFilePath = outputFilePath
            self.folder_index = folder



    class _outfiles_addSMRT:
        def __init__(self, mp_params):

            basefolder = mp_params['-outputFolderBase'] + 'addSMRT/'
            folder = mp_params['-processorIndex']
            prefix = mp_params['-outputPrefixName']
            if not os.path.exists(basefolder):
                os.makedirs(basefolder)

            outputFolder = basefolder + folder
            if not os.path.exists(outputFolder):
                os.makedirs(outputFolder)

            outputFilePath = outputFolder + '/' + prefix

            fileaddSMRT = outputFilePath + '.single.addSMRT.hap'
            print fileaddSMRT
            self.fo_addSMRT_singel = open(fileaddSMRT, 'w')

            # fileEval = outputFilePath + '.pair.eval'
            # self.fo_eval_pair = open(fileEval, 'w')

            # fileEvalNoBlock = outputFilePath + '.pair_noBlock.eval'
            # self.fo_eval_pair_NoBlock = open(fileEvalNoBlock, 'w')

            self.outputFilePath = outputFilePath
            self.folder_index = folder


    class Struct:
        def __init__(self, *argv, **argd):
            if len(argd):
                self.__dict__.update(argd)
            else:
                attrs = filter(lambda x: x[0:4] != '__', dir(self))
                for n in range(len(argv)):
                    setattr(self, attrs[n], argv[n])

    class readd(Struct):
        pos = -1
        cigar = ''
        read = ''
        name = ''
        leftRightSingle = 's'
        quality = ''

    def check_all_params(self, mp):
        l = ['-ploidy',
             '-MaxTypeLength',
             '-method',
             '-noProcessor',
             '-bamFile',
             '-refFile',
             '-vcfFile',
             '-selectedScf',
             '-outputFolderBase',
             '-mode',
             '-par',
             '-processorIndex']
        s1 = set(l) - set(mp.keys())
        if len(s1) != 0:
            print 'missing parameters:'
            print list(s1)
            return False
        else:
            self.okey_params = True
            return True

    def make_index_vcf_file(self, file_name_in):
        file_name_out = file_name_in + '.index'
        fiVCF = open(file_name_in, 'r')
        foVCF = open(file_name_out, 'w')
        lineVCF = fiVCF.readline()
        startPos = fiVCF.tell()
        while lineVCF[0] == '#' :
            startPos = fiVCF.tell()
            lineVCF = fiVCF.readline()
        cnt = 0
        while lineVCF != '':
            aVCF = lineVCF.split()
            scfName = aVCF[0]
            foVCF.write(scfName + '\t' + str(startPos) + '\n')
            # if cnt %5000 == 0:
                # print scfName+ '\t' + str(startPos) + '\t' + str(cnt)
            while lineVCF != '' and scfName == aVCF[0]:
                startPos = fiVCF.tell()
                lineVCF = fiVCF.readline()
                aVCF = lineVCF.split()
            cnt += 1

    def make_index_hap_file(self, hapFile):
        fi = open(hapFile, 'r')
        fo = open(hapFile + '.index', 'w')

        # line = fi.readline()
        startPos = fi.tell()
        line = fi.readline()
        a = line.split()
        while line != '':
            scfName = a[1]
            fo.write( scfName + '\t' + str(startPos) + '\n')
            line = fi.readline()
            while line != '' and line.strip()[:3] != '>>>':
                startPos = fi.tell()
                line = fi.readline()
            a = line.split()

    def __init__(self, argv, main_act):
        self.okey_params = False
        mp_params = {}

        def read_parameters():
            for i in range(2, len(sys.argv), 2):
                mp_params[sys.argv[i]] = sys.argv[i + 1]

            if '-par' not in mp_params:
                print 'Error: indicate -par argument (parameter file)'
                exit()

            fi = open(mp_params['-par'], 'r')
            line = fi.readline()
            while line != '':
                a = line.split()
                if a[0] not in mp_params.keys() and a[0][0] == '-':
                    if len(a) > 2:
                        mp_params[a[0]] = (a[1], a[2])
                    else:
                        mp_params[a[0]] = a[1]
                line = fi.readline()

            if '-mode' not in mp_params:
                print 'please indicate -mode parameters'
                exit()
            else:
                self.act = mp_params['-mode']

        read_parameters()
        if self.act == 'hap' and main_act == 'hap':
            self.is_eval = False
            mp_params['-MaxTypeLength'] = int(mp_params.get('-MaxTypeLength', 2)) + 1
            mp_params['-method'] = mp_params.get('-method', 'kruskal')
            mp_params['-outputPrefixName'] = mp_params.get('-outputPrefixName', 'ranbow')
            mp_params['-noProcessor'] = int(mp_params.get('-noProcessor', '1'))

            if mp_params['-noProcessor'] == 1:
                mp_params['-processorIndex'] = '0'

            if self.check_all_params(mp_params):
                self.ref = mp_params['-refFile']
                self.ref_fai = self.ref + '.fai'
                if not os.path.exists(self. ref_fai):
                    print 'No index file: ' + self.ref + '.fai'
                    print 'Run the "ranbow.py index" first'
                    exit()

                self.bam = mp_params['-bamFile']
                self.bam_bai = self.bam + '.bai'
                if not os. path.exists(self.bam_bai):
                    print 'No index file: ',self.bam + '.bai'
                    print 'Run the "ranbow.py index" first'
                    exit()

                self.vcf = open(mp_params['-vcfFile'])
                if not os.path.exists(mp_params['-vcfFile'] + '.index'):
                    print 'No index file: ' + mp_params['-vcfFile'] + '.index'
                    print 'Run the "ranbow.py index" first'
                    exit()

                self.vcf_inx = open(mp_params['-vcfFile']+'.index')
                self.no_processor = mp_params['-noProcessor']
                self.processorNumber = int( mp_params['-processorIndex'])
                self.l_selected_scfs = [i.rstrip() for i in open(mp_params['-selectedScf'])]
                self.outfiles = self. _outfiles(mp_params)
                self.algType = mp_params['-method']
                self.ploidy = int(mp_params['-ploidy'])
                self.max_type_len = int(mp_params['-MaxTypeLength'])


        if self.act == 'index' and main_act == 'hap':
            self. ref = mp_params['-refFile']
            self.ref_fai = self.ref + '.fai'
            if not os.path.exists(self.ref_fai):
                ps.faidx(self.ref)

            self.bam = mp_params['-bamFile']
            self.bam_bai = self.bam + '.bai'
            if not os.path.exists(self.bam_bai):
                ps.index(self.bam)

            self.vcf = open(mp_params['-vcfFile'])
            if not os.path.exists(mp_params['-vcfFile'] + '.index'):
                self.make_index_vcf_file(mp_params['-vcfFile'])
            self. vcf_inx = open(mp_params['-vcfFile'] + '.index')

            if not os.path.exists(mp_params['-outputFolderBase']):
                os.makedirs(mp_params['-outputFolderBase'])
            exit()

        if self.act == 'collect' and main_act == 'hap':
            self.prefix = mp_params.get('-outputPrefixName', 'ranbow')
            self.no_processor = int(mp_params.get('-noProcessor', '1'))
            self.outputFolder = mp_params[ '-outputFolderBase']
            self.ref = mp_params['-refFile']
            self.ref_fai = self.ref + '.fai'
            if not os.path.exists(self.ref_fai):
                print 'No index file: ' + self.ref + '.fai'
                print 'Run the "ranbow.py index" first'
                exit()

        if self.act == 'modVCF' and main_act == 'hap':
            self.prefix = mp_params.get('-outputPrefixName', 'ranbow')
            self.outputFolder = mp_params['-outputFolderBase']
            self.vcf = mp_params['-vcfFile']
            self.ploidy = int(mp_params['-ploidy'])

        if self.act == 'index' and main_act == 'eval':
            mp_params['-outputPrefixName'] = mp_params.get('-outputPrefixName', 'ranbow')

            self.bam = mp_params['-bamFileEval']
            self.bam_bai = self.bam + '.bai'
            if not os.path.exists(self.bam_bai):
                ps.index(self.bam)

            self.vcf = open(mp_params['-vcfFile'])
            if not os.path.exists(mp_params[ '-vcfFile'] + '.index'):
                self.make_index_vcf_file(mp_params['-vcfFile'])

            self.vcf_inx = open(mp_params['-vcfFile'] + '.index')

            basefolder = mp_params['-outputFolderBase']
            out_fo = basefolder + 'eval/'

            if not os.path.exists(out_fo):
                os.makedirs(out_fo)

            self.hapSing = basefolder+'/' + mp_params['-outputPrefixName'] + '.single.hap'
            #self.hapPair = basefolder + mp_params['-outputPrefixName'] + '.pair.hap'

            self.hapSingInx = self.hapSing + '.index'
            if not os.path.exists(self.hapSingInx):
                self.make_index_hap_file(self.hapSing)

            # self.hapPairInx = self.hapPair + '.index'
            # if not os.path. exists(self.hapPairInx):
            #     self. make_index_hap_file(self.hapPair)

        if self.act == 'collect' and main_act == 'eval':
            # self.prefix = mp_params.get('-outputPrefixName', 'ranbow')
            self.no_processor = int(mp_params.get('-noProcessor', '1'))
            self.evalOutput = mp_params['-outputFolderBase'] + 'eval/'



        if self.act == 'run' and main_act == 'eval':
            self.ploidy = int(mp_params['-ploidy'])
            self.is_eval = True
            self.is_addSMRT = False
            mp_params['-outputPrefixName'] = mp_params.get('-outputPrefixName', 'ranbow')
            mp_params['-noProcessor'] = int(mp_params.get('-noProcessor', '1'))

            if mp_params['-noProcessor'] == 1:
                mp_params['-processorIndex'] = '0'

            self.ref = mp_params['-refFile']
            self.ref_fai = self.ref + '.fai'
            if not os.path. exists(self. ref_fai):
                print 'No index file: ' + self.ref + '.fai'
                print 'Run the "ranbow.py index" first'
                exit()

            self.bam = mp_params['-bamFileEval']
            self.bam_bai = self.bam + '.bai'
            if not os.path.exists(self.bam_bai):
                print 'No index file: ' + self.bam + '.bai'
                print 'Run the "ranbow.py index" first'
                exit()

            self.vcf = open(mp_params['-vcfFile'])
            if not os.path. exists(mp_params['-vcfFile'] + '.index'):
                print 'No index file: ' + mp_params['-vcfFile'] + '.index'
                print 'Run the "ranbow.py index" first'
                exit()

            basefolder = mp_params['-outputFolderBase']
            out_fo = basefolder + 'eval/'

            self.hapSing = basefolder + '/'+  mp_params ['-outputPrefixName'] + '.single.hap'
            # self.hapPair = basefolder + mp_params['-outputPrefixName'] + '.pair.hap'

            hapSingInx = self.hapSing + '.index'
            print hapSingInx
            if not os.path.exists(hapSingInx):
                print 'no Index'
                exit()
            else:
                self.mpHapSing = {}
                for i in open(hapSingInx):
                    a = i.split()
                    self.mpHapSing[a[0]] = self.mpHapSing.get(a[0], []) + [int(a[1])]

            # hapPairInx = self.hapPair + '.index'
            # if not os.path.exists(hapPairInx):
            #     print 'no Index'
            #     exit()
            # else:
            #     self.mpHapPair = {}
            #     for i in open(hapPairInx):
            #         a = i.split()
            #         self.mpHapPair[a[0]] = self.mpHapPair.get(a[0], []) + [int(a[1])]

            self.vcf_inx = open(mp_params['-vcfFile'] + '.index')
            self.no_processor = mp_params['-noProcessor']
            self.processorNumber = int( mp_params['-processorIndex'])
            self.l_selected_scfs = [i.rstrip() for i in open(mp_params['-selectedScf'])]
            self.outfiles = self._outfiles_eval(mp_params)


        if self.act == 'run' and main_act == 'addSMRT':
            self.ploidy = int(mp_params['-ploidy'])
            self.is_addSMRT = True
            self.is_eval= False

            mp_params['-outputPrefixName'] = mp_params.get('-outputPrefixName', 'ranbow')
            mp_params['-noProcessor'] = int(mp_params.get('-noProcessor', '1'))

            if mp_params['-noProcessor'] == 1:
                mp_params['-processorIndex'] = '0'

            self.ref = mp_params['-refFile']
            self.ref_fai = self.ref + '.fai'
            if not os.path. exists(self. ref_fai):
                print 'No index file: ' + self.ref + '.fai'
                print 'Run the "ranbow.py index" first'
                exit()

            self.bam = mp_params['-bamFileSMRT']
            self.bam_bai = self.bam + '.bai'
            if not os.path.exists(self.bam_bai):
                print 'No index file: ' + self.bam + '.bai'
                print 'Run the "ranbow.py index" first'
                exit()

            self.vcf = open(mp_params['-vcfFile'])
            if not os.path. exists(mp_params['-vcfFile'] + '.index'):
                print 'No index file: ' + mp_params['-vcfFile'] + '.index'
                print 'Run the "ranbow.py index" first'
                exit()

            basefolder = mp_params['-outputFolderBase']
            out_fo = basefolder + 'addSMRT/'

            self.hapSing = basefolder +'/' +mp_params ['-outputPrefixName'] + '.single.hap'
            self.hapNanoCompare = basefolder + mp_params ['-outputPrefixName'] + '.hap_vs_SMRT'
            # self.hapPair = basefolder + mp_params['-outputPrefixName'] + '.pair.hap'
            print self.hapSing
            hapSingInx = self.hapSing + '.index'
            if not os.path.exists(hapSingInx):
                print 'no Index'
                exit()
            else:
                self.mpHapSing = {}
                for i in open(hapSingInx):
                    a = i.split()
                    self.mpHapSing[a[0]] = self.mpHapSing.get(a[0], []) + [int(a[1])]

            # hapPairInx = self.hapPair + '.index'
            # if not os.path.exists(hapPairInx):
            #     print 'no Index'
            #     exit()
            # else:
            #     self.mpHapPair = {}
            #     for i in open(hapPairInx):
            #         a = i.split()
            #         self.mpHapPair[a[0]] = self.mpHapPair.get(a[0], []) + [int(a[1])]

            self.vcf_inx = open(mp_params['-vcfFile'] + '.index')
            self.no_processor = mp_params['-noProcessor']
            self.processorNumber = int( mp_params['-processorIndex'])
            self.l_selected_scfs = [i.rstrip() for i in open(mp_params['-selectedScf'])]
            self.outfiles = self._outfiles_addSMRT(mp_params)

        if self.act == 'index' and main_act == 'addSMRT':
            mp_params['-outputPrefixName'] = mp_params.get('-outputPrefixName', 'ranbow')

            self.bam = mp_params['-bamFileSMRT']
            self.bam_bai = self.bam + '.bai'
            if not os.path.exists(self.bam_bai):
                ps.index(self.bam)

            self.vcf = open(mp_params['-vcfFile'])
            if not os.path.exists(mp_params[ '-vcfFile'] + '.index'):
                self.make_index_vcf_file(mp_params['-vcfFile'])

            self.vcf_inx = open(mp_params['-vcfFile'] + '.index')

            basefolder = mp_params['-outputFolderBase']
            out_fo = basefolder + 'addSMRT/'

            if not os.path.exists(out_fo):
                os.makedirs(out_fo)

            self.hapSing = basefolder +'/'+ mp_params['-outputPrefixName'] + '.single.hap'
            #self.hapPair = basefolder + mp_params['-outputPrefixName'] + '.pair.hap'

            self.hapSingInx = self.hapSing + '.index'
            if not os.path.exists(self.hapSingInx):
                self.make_index_hap_file(self.hapSing)

            # self.hapPairInx = self.hapPair + '.index'
            # if not os.path. exists(self.hapPairInx):
            #     self. make_index_hap_file(self.hapPair)

    def make_mp_scfLen(self):
        mp = {}
        for i in open(self.ref_fai):
            a = i.split()
            mp[a[0]] = int(a[1])
        return mp

    def distribute_selected_scaffolds(self):
        mp = self.make_mp_scfLen()

        l_scf = []
        su = 0
        for i in self.l_selected_scfs:
            a = i.split()
            if len(a) == 4:
                l_scf += [(a[0], int(a[1]), int(a[2]), a[3])]
            else:
                l_scf += [(a[0], 1, int(mp[a[0]]), a[0])]
            su += l_scf[-1][2] - l_scf[-1][1] + 1

        #print self.l_selected_scfs
        #print su
        #print self.no_processor
        avr = su/self.no_processor

        l_l_selected_scf = []
        l_l_sum = []
        inx = 0
        for i in xrange(min(self.no_processor, len(l_scf))):
            l_l_selected_scf.append([l_scf[inx]])
            l_l_sum.append(l_scf[inx][2] - l_scf[inx][1] + 1)
            inx += 1
        while inx < len(l_scf):
            ch = False
            for i in range (self.no_processor):
                if l_l_sum[i] > avr:
                    continue
                else:
                    ch = True
                    l_l_selected_scf[i].append(l_scf[inx])
                    l_l_sum[i] += l_scf[inx][2] - l_scf[inx][1] + 1
                    inx += 1
                if inx == len(l_scf):
                    break

            if not ch:
                for j in range(inx, len(l_scf)):
                    inx_min = l_l_sum.index(min(l_l_sum))
                    l_l_sum[inx_min] += l_scf[inx][2] - l_scf[inx][1] + 1
                    l_l_selected_scf[inx_min].append(l_scf[inx])
                break
        return l_l_selected_scf

    def load_scaffolds_length(self):
        l = []
        su = 0
        for i in open(self.ref_fai):
            a = i.split()
            l += [(a[0], int(a[1]))]
            su += int(a[1])
        return l, su/self.no_processor

    def distribute_scaffolds(self, l_scf_len, avr):
        print l_scf_len[0]
        exit()
        l_l_selected_scf = []
        l_l_sum = []
        inx = 0
        for i in xrange(self.no_processor):
            l_l_selected_scf.append([(l_scf_len[inx][0], l_scf_len[inx][1])])
            l_l_sum.append(l_scf_len[inx][1])
            inx += 1

        while inx < len(l_scf_len):
            ch = False
            for i in range(self.no_processor):
                if l_l_sum[i] > avr:
                    continue
                else:
                    ch = True
                    l_l_selected_scf[i].append((l_scf_len[inx][0], l_scf_len[inx][1]))
                    l_l_sum[i] += l_scf_len[inx][1]
                    inx += 1
                if inx == len(l_scf_len):
                    break
            if not ch:
                for j in range(inx, len(l_scf_len)):
                    inx_min = l_l_sum.index(min(l_l_sum))
                    l_l_sum[inx_min] += l_scf_len[j][1]
                    l_l_selected_scf[inx_min].append((l_scf_len[inx][0], l_scf_len[inx][1]))
                break

        return l_l_selected_scf

    def load_inx_file(self, f):
        mp = {}
        for i in f:
            a = i. split()
            mp[a[0]] = int(a[1])
        return mp

    def indexingVCF(self, scName, fi):
        ploidy_no = self.ploidy
        arrInx2pos = []
        arrInxVariation = []
        arrVarPattern = []

        line = fi.readline()
        a = line.split('\t')

        while a[0] == scName:
            allele_pos = int(a[1])
            end_ref_allele = allele_pos + len(a[3])
            l_allele = [a[3]] + a[4].split(',')

            arrInx2pos.append(allele_pos)
            arrInxVariation.append(l_allele[:ploidy_no])

            tmp1 = a[7]. split(';')
            l_cigar = tmp1[6][6:]. split(',')
            arrVarPattern.append(l_cigar[:ploidy_no-1])

            line = fi.readline()
            a = line.split('\t')
            # for those consecutive SNPs that one covers the other
            while a[0] == scName and end_ref_allele > int(a[1]):
                line = fi.readline()
                a = line.split('\t')

        return arrInx2pos, arrInxVariation, arrVarPattern

    def makeMapRefandRead(self, r):
        mp_read_ref = {}
        mp_ref_read = {}
        arr_cigar = re.split('(\d+)', r.cigar)
        new_cigar = ''
        inxRef = r.pos - 1
        inxRead = -1
        indel_tag = False
        for i in range(1,len(arr_cigar), 2):
            if arr_cigar[i + 1] in ['M','X','=']:
                new_cigar += arr_cigar[i] + 'M'
                no_match = int(arr_cigar[i])

                if indel_tag:
                    no_match -= 1
                    indel_tag = True

                for j in range(no_match):
                    inxRef += 1
                    inxRead += 1
                    mp_read_ref[inxRead] = inxRef

            elif arr_cigar[i+1] == 'D':
                new_cigar += arr_cigar[i] + 'D'
                inxRef += int(arr_cigar[i]) + 1
                inxRead += 1
                mp_read_ref[inxRead] = inxRef
                indel_tag = True

            elif arr_cigar[i+1] == 'I':
                new_cigar += arr_cigar[i] + 'I'
                inxRead += int(arr_cigar[i]) + 1
                inxRef += 1
                mp_read_ref [inxRead] = inxRef
                indel_tag = True

            elif arr_cigar[i+1] == 'S':
                # inxRef += int(arr_cigar[i])
                inxRead += int(arr_cigar[i])

        delList = []

        for i in mp_read_ref:
            if i >= len(r.read):
                delList.append(i)
        for i in delList:
            del mp_read_ref[i]

        for i in mp_read_ref:
            mp_ref_read[mp_read_ref[i]] = i

        r.cigar = new_cigar
        return mp_read_ref, mp_ref_read

    def cigar2length(self, c):
        #		M	0*  alignment match (can be a sequence match or mismatch)
        #		I	1*  insertion to the reference
        #		D	2*  deletion from the reference
        #		N	3-  skipped region from the reference
        #		S	4-  soft clipping (clipped sequences present in	SEQ)
        #		H	5-  hard clipping (clipped sequences NOT present in SEQ)
        #		P	6-  padding (silent deletion from padded reference)
        #		=	7*  sequence match
        #		X	8*  sequence mismatch


        ans = 0
        arr_cigar = re.split('(\d+)', c)
        for i in range(1, len(arr_cigar), 2):
            if arr_cigar[i + 1] == 'M':
                ans += int(arr_cigar[i])

            elif arr_cigar[i + 1] == 'D':
                ans += int(arr_cigar[i])

            elif arr_cigar[i+1] == 'I':
                pass

            elif arr_cigar[i+1] == 'S':
                pass
                # ans += int(arr_cigar[i])

            elif arr_cigar[i + 1] == 'X':
                ans += int(arr_cigar[i])

            elif arr_cigar[i + 1] == '=':
                ans += int(arr_cigar[i])

            elif arr_cigar[i + 1] == 'H':
                pass

            elif arr_cigar[i + 1] == 'N':
                ans += int(arr_cigar[i])

            else:
                print 'unexpected cigar' + c
                exit()
        return ans

    def make_segments_for_read(self, r, vcf_info):
        import numpy as np



        QualityThreshold = 10


        arr_pos = vcf_info [0]
        arr_var = vcf_info [1]

        mp_r_M, mp_M_r = self.makeMapRefandRead(r)

        #print mp_r_M, mp_M_r

        segList = {}
        segQual = {}
        seqLen = self.cigar2length(r.cigar)

        left_pos = bisect_left(arr_pos, r.pos, 0, len(arr_pos)) - 1
        right_pos = bisect_right(arr_pos, r.pos + seqLen, 0, len(arr_pos)) + 1

        for i in range(max(0, left_pos), min(len(arr_pos), right_pos)):
            var_s_M = arr_pos[i] - 1  # starts before the variant
            var_e_M = arr_pos[i] + len(arr_var[i][0])  # ends after the variant

            bo_starts_in = True if (var_s_M in mp_M_r) else False
            bo_ends_in = True if (var_e_M in mp_M_r) else False

            cnt_matched = 0
            if bo_starts_in  and bo_ends_in:
                pos_s_r = mp_M_r[var_s_M]
                pos_e_r = mp_M_r[var_e_M]
                seg = r.read[pos_s_r+1 : pos_e_r]
                seg_q  = 0 if len(r.qual[pos_s_r+1 : pos_e_r]) == 0 else int(np.mean(r.qual[pos_s_r+1 : pos_e_r]))


                if seg_q < QualityThreshold:
                    continue


                for inxV in range(len(arr_var[i])):
                    v = arr_var[i][inxV]
                    if seg == v:
                        cnt_matched += 1
                        tmp_inxV = inxV
                        tmp_Vq = seg_q

                if cnt_matched == 1:
                    segList[i] = tmp_inxV
                    segQual[i] = tmp_Vq

            elif bo_starts_in and not bo_ends_in:
                pos_s_r = mp_M_r[var_s_M]
                seg = r.read[pos_s_r+1:]
                seg_q  = 0 if len(r.qual[pos_s_r+1:]) == 0 else int(np.mean(r.qual[pos_s_r+1:]))

                if seg_q < QualityThreshold:
                    continue

                for inxV in range(len(arr_var[i])):
                    v = arr_var[i][inxV]
                    if v.startswith(seg):
                        cnt_matched += 1
                        tmp_inxV = inxV
                        tmp_Vq = seg_q

                if cnt_matched == 1:
                    segList[i] = tmp_inxV
                    segQual[i] = tmp_Vq

            elif not bo_starts_in and bo_ends_in:
                pos_e_r = mp_M_r[var_e_M]
                seg = r.read[:pos_e_r]
                seg_q  = 0 if len(r.qual[:pos_e_r]) == 0 else int(np.mean(r.qual[:pos_e_r]))
                if seg_q < QualityThreshold:
                    continue

                for inxV in range(len(arr_var[i])):
                    v = arr_var[i][inxV]
                    if v[:: -1].startswith(seg[::-1]):
                        cnt_matched += 1
                        tmp_inxV = inxV
                        tmp_Vq = seg_q


                if cnt_matched == 1:
                    segList[i] = tmp_inxV
                    segQual[i] = tmp_Vq



        return segList, segQual

    def segment2frag(self, mp_ReadSegments):
        mn = min(mp_ReadSegments.keys())
        mx = max(mp_ReadSegments.keys())
        fr_t = (mx - mn + 1) * ['-']
        for i in mp_ReadSegments:
            fr_t[i - mn] = str(mp_ReadSegments[i])
        fr = ''.join(fr_t)
        return mn, fr

    def segment2tuple(self, mp_ReadSegments,mp_ReadSegments_qual):

        mn = min(mp_ReadSegments.keys())
        frr = '{'
        frr_v = '{'
        frr_q = '{'
        for j in sorted(mp_ReadSegments.keys()):
            frr += str(j) + ':' + str(mp_ReadSegments[j]) + ','
            frr_q += str(j) + ':' + str(mp_ReadSegments_qual[j]) + ','
        frr = frr[:-1] + '}'
        frr_q = frr_q[:-1] + '}'
        return (mn, frr, max(mp_ReadSegments.keys()) + 1,frr_q)



    def segment2frag_seg(self, mp_ReadSegments):
        mn = min(mp_ReadSegments.keys())
        mx = max(mp_ReadSegments.keys())
        fr_t = (mx - mn + 1) * ['-']
        for pos in mp_ReadSegments:
            fr_t[pos - mn] = str(mp_ReadSegments[pos][0])
        fr = ''.join(fr_t)
        return mn, fr

    def merge_two_seg(self,l):
        pass

    def convert_to_uniq(self, l_frag):

        frag_l = Counter(l_frag).keys()
        count_l = Counter(l_frag).values()
        l_frag_uniq = [(j, i[0], i[1]) for i, j in zip(frag_l, count_l)]
        l_frag_uniq.sort(key=lambda x: (x[1], x[0]))
        return l_frag_uniq

    def convert_to_uniq_pair(self, l_frag):

        l_frag.sort(key=lambda x: x[0])
        l_frag_pair = []
        for i in xrange(1, len(l_frag)):
            r1 = l_frag[i-1]
            r2 = l_frag[i]
            if r1[0] == r2[0]:
                r1_pos = r1[1][0]
                r2_pos = r2[1][0]

                r1_str = r1[1][1]
                r2_str = r2[1][1]

                r1_end = r1_pos + len(r1_str)
                r2_end = r2_pos + len(r2_str)

                mn = min(r1_pos, r2_pos)
                mx = max(r1_end, r2_end)
                ln = mx - mn
                fr = ln * ['-']

                for j in xrange(r1_pos, r1_end):
                    fr[j - mn] = r1_str[j - r1_pos]

                bo = False
                for j in xrange(r2_pos, r2_end):
                    if r2_str[j - r2_pos] == '-':
                        pass
                    elif fr[j - mn] == '-':
                        fr[j - mn] = r2_str[j - r2_pos]
                    elif fr[j - mn] != r2_str[j - r2_pos]:
                        fr[j - mn] = '-'

                frr = ''.join(fr)
                if ln - frr.count('-') >= 2:
                    l_frag_pair += [(mn, frr)]

        frag_l = Counter(l_frag_pair).keys()
        count_l = Counter(l_frag_pair).values()
        l_frag_uniq = [(j, i[0], i[1]) for i, j in zip(frag_l, count_l)]
        l_frag_uniq.sort(key=lambda x: (x[1], x[0]))
        return l_frag_uniq

    def convert_to_uniq_pair_from_seg(self, l_seg):

        l_seg.sort(key=lambda x: x[0])
        l_frag_pair = []

        for i in xrange(1, len(l_seg)):
            t1 = l_seg[i - 1]
            t2 = l_seg[i]
            if t1[0] == t2[0]:
                print 'here'
                exit()

                r1 = self.segment2frag(ast.literal_eval(t1[1]))#not reachable code to change ast library
                r2 = self.segment2frag(ast.literal_eval(t2[1]))

                r1_pos = r1[0]
                r2_pos = r2[0]

                r1_str = r1[1]
                r2_str = r2[1]

                r1_end = r1_pos + len(r1_str)
                r2_end = r2_pos + len(r2_str)

                mn = min(r1_pos, r2_pos)
                mx = max (r1_end, r2_end)
                ln = mx - mn
                fr = ln * ['-']

                for j in xrange(r1_pos, r1_end):
                    fr[j - mn] = r1_str[j - r1_pos]

                bo = False
                for j in xrange(r2_pos, r2_end):
                    if r2_str[j - r2_pos] == '-':
                        pass
                    elif fr [j - mn] == '-':
                        fr [j - mn] = r2_str[j - r2_pos]
                    elif fr [j - mn] != r2_str[j - r2_pos]:
                        fr [j - mn] = '-'

                while len(fr) > 0 and fr[0] == '-':
                    mn += 1
                    del fr[0]

                while len( fr) > 0 and fr[ -1] == '-':
                    del fr[-1]

                frr = ''.join(fr)
                if len(frr) - frr.count('-') >= 2:
                    l_frag_pair += [(mn, frr)]

        A = Counter(l_frag_pair)
        l_frag_uniq = [(value , key[0],key[1]) for key, value in A.iteritems()]
        l_frag_uniq.sort(key=lambda x: (x[1], x[0], x[2]))
        return l_frag_uniq

    def convert_to_uniq_pairSEG_from_seg(self, l_seg, l_seg_single):
        import numpy as np

        l_seg.sort(key=lambda x: x[0])
        l_frag_pair = []

        for i in xrange(1, len(l_seg)):
            t1 = l_seg[i - 1]
            t2 = l_seg[i]
            if t1[0] == t2[0]:
                # r1 = self.segment2frag(ast.literal_eval(t1[1]))
                # r2 = self.segment2frag(ast.literal_eval(t2[1]))
                # mp1 = ast.literal_eval(t1[1]) #here
                # mp2 = ast.literal_eval(t2[1])
                mp1 = t1[1] #here
                mp2 = t2[1]

                mp1_q = t1[2]
                mp2_q = t2[2]

                mp_all = dict(mp1)
                mp_all_q = dict(mp1_q)
                for j in mp2:
                    if j in mp_all:
                        if mp_all[j] == mp2[j]:
                            mp_all_q[j] = mp_all_q[j] + mp2_q[j]
                        else:
                            if mp_all_q[j] < mp2_q[j]:
                                mp_all[j] = mp2[j]
                                mp_all_q[j] = mp2_q[j] - mp_all_q[j]
                            else :
                                mp_all_q[j] = mp_all_q[j] - mp2_q[j]
                    else:
                        mp_all[j] = mp2[j]
                        mp_all_q[j] = mp2_q[j]
                if len(mp_all) >= 2:
                    mn = min(mp_all.keys())
                    frr = '{'
                    frrq = '{'
                    for j in sorted(mp_all.keys()):
                        frr += str(j) + ':' + str(mp_all[j]) + ','
                        frrq += str(j) + ':' + str(mp_all_q[j]) + ','

                    frr = frr[:-1] + '}'
                    frrq = frrq[:-1] + '}'
                    l_frag_pair += [(mn, frr, max(mp_all.keys()) + 1, frrq )]
                    # l_frag_pair += [(mn, mp_all)]

        #A = Counter(l_seg_single)

        A = Counter(l_frag_pair+l_seg_single)
        # B = Counter(l_frag_pair)
        # print len(A),len(B)
        l_frag_uniq = [(cpy, frag[0], ast.literal_eval(frag[1]), frag[2],ast.literal_eval(frag[3])) for frag, cpy in A.iteritems()] #here
        # l_frag_uniq = [(value, key[0], key[1], max (key[1]) + 1 ) for key, value in A.iteritems()] #here
        l_frag_uniq.sort(key=lambda x: (x[1], x[3], -x[0]))
        return l_frag_uniq


    def refine_haplotype_length(self,l,scfName,folderIndex):

        for b in l:
            for i in  b:
                i.my_initial_blocks = i.name








        inx_b = 0

        for b in l:
            ploidy = 0
            inx_h = 0
            for i in b:
                if i.flag != -1 :
                    ploidy += 1

            for i in b:
                if i.flag != -1 :
                    i.name = scfName+"_ploidy_"+str(ploidy)+"_Bl_" + str(inx_b) + '_Hap_'+str(inx_h)+'_foInx_'+folderIndex
                    inx_h += 1
            inx_b += 1


    def write_subhaps_to_file(self, l, fo, scfName, scfLen):

        fo.write('>>>\t'+scfName + '\t' + str(scfLen) + '\n')
        for b in l:
            for i in b:
                if i.flag != -1 :
                    mytype_pos = ''.join([str(t)+'-' for t in i.type if t!= -1]) # -1 is type separator
                    mytype = ''.join([str(i.haplotype[t][0]) for t in i.type if t != -1])
                    hap_seg = str(i.haplotype).replace(" ", "")
                    hap_sup = str(i.supporting_fragments).replace(" ","")
                    fo.write(str(i.flag) + '\t' + i.name + '\t' + str(i.start_pos) + '\t'
                             +   self.segment2frag_seg(i.haplotype)[1]
                             + '\t' + str(i.corrected_error) + '\t'+ i.my_initial_blocks + '\tQ:'+ hap_seg + '\tS:' + hap_sup
                             + '\t' + mytype + '\t'+ mytype_pos + '\t' + str(i.block_no).replace(" ",'') + '\n')

    def combined_cigar(self, frag_cigar):

        frag_cigar = frag_cigar.replace('X','M')
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
            ans += str(no) + ch
        return ans

    def extract_dna_of_a_fragment(self, s, hap_name, scf, vcf_info):
        l_allele__index_2_realPos = vcf_info[0]
        l_index_2_alleles = vcf_info[1]
        l_index_2_cigar = vcf_info[2]
        name = hap_name
        frag_start = s.start_pos
        frag = self.segment2frag_seg(s.haplotype)[1]
        scf = '-' + scf

        ref_allele = 0
        if frag_start == 0:
            frag_nc = scf[1:l_allele__index_2_realPos[0]]
            frag_nc_start_pos = 1
            frag_cigar = str(l_allele__index_2_realPos[0] - 1) + 'M'
        else:
            frag_nc_start_pos = l_allele__index_2_realPos[frag_start - 1] + len(l_index_2_alleles[frag_start - 1][ref_allele])
            frag_nc = scf[frag_nc_start_pos:l_allele__index_2_realPos[frag_start]]
            frag_cigar = str(l_allele__index_2_realPos[frag_start] - frag_nc_start_pos) + 'M'

        end_pos_on_noneallelic_region = start_pos_on_noneallelic_region = -1
        for i in range(len(frag)):
            if frag[i] == '0':
                allele = l_index_2_alleles[frag_start + i][int(frag[i])]
                frag_nc += allele
                frag_cigar += str(len(allele)) + 'M'

            elif frag[i] == '-':
                alleleSize = len(l_index_2_alleles[frag_start + i][ref_allele])
                frag_nc += alleleSize * 'N'
                frag_cigar += str(alleleSize) + 'M'

            else:
                frag_cigar += l_index_2_cigar[frag_start + i][int(frag[i]) - 1]
                frag_nc += l_index_2_alleles[frag_start + i][int(frag[i])]

            start_pos_on_noneallelic_region = l_allele__index_2_realPos[frag_start + i] \
                                              + len(l_index_2_alleles[frag_start + i][ref_allele])

            if frag_start + i + 1 < len(l_index_2_alleles):
                end_pos_on_noneallelic_region = l_allele__index_2_realPos[frag_start + i + 1]
            elif frag_start + i + 1 == len(l_index_2_alleles):
                end_pos_on_noneallelic_region = len(scf)
            else:
                print 'error: fragment passed the boundary of scaffold'

            frag_nc += scf[start_pos_on_noneallelic_region:end_pos_on_noneallelic_region]
            frag_cigar += str(end_pos_on_noneallelic_region - start_pos_on_noneallelic_region) + 'M'

        return frag_nc_start_pos, frag_nc, self.combined_cigar(frag_cigar)

    def hap2fasta_and_sam_no_SNP(self, ref, scf_name, fo, ploidy):
        for i in range(ploidy):
            hap_name = scf_name + '_hap' + str(i)
            t_dna = (0, ref, str(len(ref)) + 'M')
            fo.write(hap_name + '\t16\t' + scf_name + '\t' + str(t_dna[0] + 1) + '\t100\t' + t_dna[2] + '\t*\t0\t0\t'
                     + t_dna[1] + '\t' + len(t_dna[1])*'9' + '\tNM:i:9\tRG:Z:454\n')

    def hap2fasta_and_sam(self, vcf_info, ref, l_subhap, scf_name, fo, fo_fasta, scf_title):
        l_allele__index_2_realPos, l_index_2_alleles, l_index_2_cigar = vcf_info

        inx_hap = 0
        inx_blc = 0
        for b in l_subhap:
            for i in b:
                if i.flag != -1 :
                    hap_name = i.name
                    t_dna = self.extract_dna_of_a_fragment(i, hap_name, ref, vcf_info)
                    fo.write(hap_name + '\t16\t' + scf_name + '\t' + str(t_dna[0]) + '\t100\t' + t_dna[2] + '\t*\t0\t0\t'
                             + t_dna[1] + '\t' + len(t_dna[1])*'9' + '\tNM:i:9\tRG:Z:454\n')
                    fo_fasta.write('>' + hap_name + '\n' + t_dna[1] + '\n')
                    inx_hap += 1
            inx_blc += 1
            inx_hap = 0
        fo.flush()
        fo_fasta.flush()
        return

    def compute_score_hap_454(self, hap, seq454):
        no_sim , no_dissim = 0, 0
        bo_overlap = False
        for i in set(hap.haplotype.keys()) & set(seq454[1].keys()):
            bo_overlap = True
            if seq454[1][i] == int(hap.haplotype[i]):
                no_sim += 1
            else:
                no_dissim += 1
        return (no_sim, no_dissim),bo_overlap

    def evaluate_haplotype_blocks(self, l_454, scf, mpHap, fo, hapFile):
        import methods as mth

        ploidy = self.ploidy
        fi = open(hapFile)

        l_subhap = []
        for i in mpHap[scf]:
            fi.seek(i)
            hap_scf = fi.readline().split()[1]
            line = fi.readline()
            while line != '' and line[:3] != '>>>' and scf == hap_scf:
                a = line.split()
                #
                # def __init__(self, haplotype, supporting_fragments, flag, corrected_error, name, type, block_no,
                #              frag_list_to_merge_later):
                convert = mth._convert()
                seg = convert.frag2seq(a[3], int(a[2]))
                hap = mth._subhap_seg(seg, [], int(a[0]), -1, a[1], (), -1, [])
                l_subhap.append(hap)
                line = fi.readline()

        l_error = []
        for i in range(len(l_subhap)):
            l_error.append([])

        cnt = 1

        for i_block in xrange(len(l_subhap)/ploidy):
            mn_main = sys.maxint
            mx = -1
            for i_sp in range(ploidy):
                if l_subhap[i_block*ploidy + i_sp].start_pos < mn_main:
                    mn_main = l_subhap[i_block*ploidy + i_sp].start_pos
                if l_subhap[i_block*ploidy + i_sp].end_pos > mx:
                    mx = l_subhap[i_block*ploidy + i_sp].end_pos

            for i_454 in range(len(l_454)):
                frg = mth._fragment(l_454[i_454][2], l_454[i_454][1], l_454[i_454][0])
                mn = mn_main

                l_agreement = []
                for i_sp in range(ploidy):
                    tmp_tuple = self.compute_score_subhap_454(l_subhap[i_block*ploidy + i_sp], frg)
                    l_agreement.append(tmp_tuple)

                bo_conflict = True
                no_conflict = 0
                for tu in l_agreement:
                    if tu[0] == 0 and tu[1] == 0:
                        bo_conflict = False
                    # if tu[1] > 0:
                        # no_conflict += 1

                if bo_conflict:  # and no_conflict == ploidy:
                    # sd similarity dissimilarity score
                    mx_scr = l_agreement[0][0] - l_agreement[0][1]**2

                    inx_mx_scr = 0
                    for inx_l in range(len(l_agreement)):
                        scr = l_agreement[inx_l][0] - l_agreement[inx_l][1]**2
                        if scr > mx_scr:
                            mx_scr = scr
                            inx_mx_scr = inx_l

                    if frg.start_pos < mn:
                        mn = frg.start_pos
                    if frg.end_pos > mx:
                        mx = frg.end_pos

                    fo.write(str(cnt) + '\n')
                    cnt += 1
                    tmp_hap = '-'*(frg.start_pos - mn) + frg.fragment + '-'*(mx - frg.end_pos)
                    fo.write('454\n' + '{0:5d}'.format(frg.start_pos) + '\t' + tmp_hap + '\n')
                    fo.write('illumina\n')
                    for i_sp in range(ploidy):
                        tmp_hap = '-'*(l_subhap[i_block*ploidy + i_sp].start_pos - mn)\
                                  + l_subhap[i_block * ploidy + i_sp].haplotype\
                                  + '-'*(mx - l_subhap[i_block*ploidy + i_sp].end_pos)
                        fo.write('{0:5d}'.format(l_subhap[i_block*ploidy + i_sp].start_pos) + '\t' + tmp_hap + '\t'
                                 + str(l_agreement[i_sp]))
                        if i_sp == inx_mx_scr:
                            fo.write('     <--- \n')
                            # if l_agreement[i_sp] in mp_all_result:
                                #mp_all_result[l_agreement[i_sp]] += 1
                            # else:
                                # mp_all_result[l_agreement[i_sp]] = 1
                        else:
                            fo.write('\n')

    def evaluate_haplotype_not_block(self, l_454, scf, mpHap, fo, hapFile):


        if scf not in mpHap:
            return
        import methods as mth

        fi = open(hapFile)

        l_subhap = []
        cnt_hap_for_name = 0
        for i in mpHap[scf]:
            fi.seek(i)
            hap_scf = fi.readline().split()[1]
            line = fi.readline()
            while line != '' and line[:3] != '>>>' and scf == hap_scf:
                a = line.split()
                convert = mth._convert()
                seg = convert.frag2seq(a[3], int(a[2]))
                hap = mth._subhap_seg(seg, [], int(a[0]), -1, a[1], (), -1, [])
                cnt_hap_for_name +=1
                l_subhap.append(hap)
                line = fi.readline()
        l_error = []
        for i in range(len(l_subhap)):
            l_error.append([])

        for i_454, seg_454 in enumerate(l_454):
            best_scr = -1000000
            bo_overlapWithOne =False
            for i_hap, hap  in enumerate(l_subhap):
                scr,bo_overlap = self.compute_score_hap_454(hap, seg_454)
                bo_overlapWithOne = True if bo_overlap else bo_overlapWithOne
                if bo_overlap:
                    scr_cur = scr[0] - scr[1]**2
                    if scr_cur > best_scr:
                        best_scr = scr_cur
                        simDis_best = scr
                        names = (seg_454[0], hap.name)
                    fo.write(">>>>\t"+ seg_454[0]+ '\t' + hap.name + '\t' + str(scr[0]) + '\t' + str(scr[1])+ '\t' + str(scr_cur) + '\n')

            if bo_overlapWithOne:
                fo.write("++++\t"+names[0]+'\t'+names[1]+'\t'+ str(simDis_best[0])+'\t'+ str(simDis_best[1])+'\t' + str(best_scr) +'\n')

    def merge_hap_and_nano(self, l_454, scf, mpHap, fo, hapFile):

        if scf not in mpHap:
            return
        import methods as mth

        fi = open(hapFile)

        l_subhap = []
        for i in mpHap[scf]:
            fi.seek(i)
            hap_scf = fi.readline().split()[1]
            line = fi.readline()
            while line != '' and line[:3] != '>>>' and scf == hap_scf:
                a = line.split()
                l_subhap += [[]]
                for j in range(int(a[1].split('_')[2])):
                    convert = mth._convert()
                    seg = convert.frag2seq(a[3], int(a[2]))
                    hap = mth._subhap_seg(seg, [], int(a[0]), -1, a[1], (), -1, [])
                    hap.my_initial_blocks = [ int(b.split('_')[0]) for b in a[5].split('__') ]
                    l_subhap[-1].append(hap)
                    line = fi.readline()
                    a = line.split()


        l_connection = []
        for i_454, seg_454 in enumerate(l_454):
            m_lower = 5
            m_rate = 0.2
            best_scr = -1000000

            for i_b, b in enumerate(l_subhap):
                cnt_overlap = 0
                for i_hap, hap in enumerate (b):
                    (m,mm), bo_overlap = self.compute_score_hap_454(hap, seg_454)
                    cnt_overlap += 1 if bo_overlap else 0
                    if bo_overlap and m >=m_lower and (float(mm)/(m+mm))<m_rate:
                        scr_cur = m - mm ** 2
                        if scr_cur > best_scr:
                            best_scr = scr_cur
                            simDis_best = (m,mm)
                            names = (seg_454[0], hap.name)
                        l_connection += [(scr_cur,m,mm,i_454,i_b,i_hap,seg_454[0], hap.name,hap.my_initial_blocks)]
                        #fo.write(">>>>\t" + seg_454[0] + '\t' + hap.name + '\t' + str(scr[0]) + '\t' + str(scr[1]) + '\t' + str(scr_cur) + '\n')


        l_connection.sort(key=lambda x:x[3])

        #for ii in l_connection:


        n = len(l_connection)
        print n
        for i in range(n-1):
            for j in range(i+1,n):
                if l_connection[i][3]==l_connection[j][3]:
                    if len(set(l_connection[i][-1]) & set(l_connection[j][-1])) ==0:
                        fo.write(str(l_connection[i]) + '\t' + str(l_connection[j])+'\n')

        import methods
        m = methods.methods()
        m.find_3_4angular_justconnections_for_nanopore_not_completed(l_subhap,l_connection,7,6)

        exit()


            # if bo_overlapWithOne:
            #     fo.write("++++\t" + names[0] + '\t' + names[1] + '\t' + str(simDis_best[0]) + '\t' + str(
            #         simDis_best[1]) + '\t' + str(best_scr) + '\n')
            #
#
    # def haplotyping(self, l):
    #
    #     fastafile = pysam.FastaFile(self.ref)
    #     samfile = pysam.AlignmentFile(self.bam, 'rb')
    #     mp_vcf = self.load_inx_file(self.vcf_inx)
    #
    #     for scfinfo in l:
    #         print scfinfo
    #         scf = scfinfo[0]
    #         scf_start = scfinfo[1]
    #         scf_end = scfinfo[2]
    #         scf_title = scfinfo[3]
    #
    #         refSeq = fastafile.fetch(scf)
    #         if scf not in mp_vcf.keys():
    #             if self.is_eval:
    #                 continue
    #
    #             t_time0 = datetime.datetime.now()
    #             self.hap2fasta_and_sam_no_SNP(refSeq, scf, self.outfiles.fo_hap_pair_sam, self.ploidy)
    #             t_time1 = datetime. datetime.now()
    #             print 'No SNP'
    #             print '>>>\t' + scf + '\treadBAM\t0:00:00.0\tRanbow_single\t0:00:00.0\tsingle2file\t0:00:00.0\tpair_data\t0:00:00.0\tranbow_pair\t0:00:00.0\tpair_2file\t0:00:00.0\n'
    #             continue
    #
    #         noR2hits = noPairR2hits = HIP_all = noR1hit = ReadWith1Hit = noR1plusHits = 0
    #
    #         self.vcf.seek(mp_vcf[scf])
    #         vcf_info = self.indexingVCF(scf, self.vcf)
    #         scf_hap_len = len(vcf_info[0])
    #
    #         l_frag_single = []
    #         l_frag_pair = []
    #         l_seg_pair = []
    #         t1 = datetime.datetime.now()
    #
    #         cnt1 = 0
    #         for read in samfile.fetch(scf, scf_start - 1, scf_end + 1):
    #             aSAM = str(read).split()
    #             if aSAM[5] != '*':
    #                 lfs = '-'
    #                 lfs = 's' if not read.is_paired else 'r' if read.is_reverse else 'l'
    #                 r = self.readd(pos=read.pos + 1, cigar=aSAM[5], read=read.query_sequence, name=read.query_name,
    #                                leftRightSingle=lfs)
    #                 mp_ReadSegments = self.make_segments_for_read(r, vcf_info)
    #
    #                 if len(mp_ReadSegments) >= 1 and read.is_paired:
    #                     l_seg_pair += [(read.query_name,str(mp_ReadSegments))]
    #
    #                 if len(mp_ReadSegments) >= 2:
    #                     l_frag_single += [self.segment2frag(mp_ReadSegments)]
    #                     # print l_frag_single
    #                     # exit()
    #
    #                 cnt1 += 1
    #
    #         l_frag_single_uniq = self.convert_to_uniq(l_frag_single)
    #         del l_frag_single
    #
    #         if self.is_eval:
    #             if len(l_frag_single_uniq) > 0:
    #                 self.evaluate_haplotype_blocks(l_frag_single_uniq, scf, self.mpHapSing,
    #                                                self.outfiles.fo_eval_singel, self.hapSing)
    #                 # self.evaluate_haplotype_blocks(l_frag_single_uniq, scf, self.mpHapPair, self.outfiles.fo_eval_pair, self.hapPair)
    #                 # self.evaluate_haplotype_not_block(l_frag_single_uniq, scf, self.mpHapPair, self.outfiles.fo_eval_pair_NoBlock, self.hapPair)
    #                 # self.evaluate_haplotype_not_block(l_frag_single_uniq, scf, self.mpHapPair,
    #                 #                                   self.outfiles.fo_eval_pair, self.hapPair)
    #             continue
    #
    #         per_single_readBAM = datetime.datetime.now() - t1
    #         print 'Run ranbow - single mode'
    #         t1 = datetime.datetime.now()
    #
    #         import methods
    #         m = methods.methods()
    #         m.algName = self.algType
    #         m.algorithm_mode = 'single'
    #         m.uniq_frag_mode = 'frag'
    #         m.ploidy = self.ploidy
    #         m.max_type_len = self.max_type_len
    #         m.single_frag = l_frag_single_uniq
    #         m.pair_frag = []
    #         m.scf_hap_len = scf_hap_len
    #
    #         m.run_method_()
    #
    #         per_single_run = datetime.datetime.now() - t1
    #         t1 = datetime.datetime.now()
    #
    #         self.write_subhaps_to_file(m.l_haps_single, self.outfiles.fo_hap_single, scf, scf_hap_len)
    #         self.hap2fasta_and_sam(vcf_info, refSeq, m.l_haps_single, scf, self.outfiles.fo_hap_single_sam,
    #                                self.outfiles.fo_hap_single_fasta, scf_title)
    #         del l_frag_single_uniq
    #         per_single_run2file = datetime.datetime.now() - t1
    #         print 'Read data for Paired-end mode'
    #         t1 = datetime.datetime.now()
    #
    #         m.algorithm_mode = 'pair'
    #         m.uniq_frag_mode = 'seg'
    #         if m.uniq_frag_mode == 'seg':
    #             m.l_frag_pair_uniq = self.convert_to_uniq_pairSEG_from_seg(l_seg_pair)
    #         else:
    #             m.l_frag_pair_uniq = self.convert_to_uniq_pair_from_seg(l_seg_pair)
    #         del l_seg_pair
    #
    #         per_dataPE = datetime.datetime.now() - t1
    #         print 'Run Ranbow - paired mode'
    #         t1 = datetime.datetime.now()
    #
    #         m.run_method_()
    #
    #         per_PE = datetime.datetime.now() - t1
    #         t1 = datetime.datetime.now()
    #
    #         self.write_subhaps_to_file(m.l_haps_pair, self.outfiles.fo_hap_pair, scf, scf_hap_len)
    #         self.hap2fasta_and_sam(vcf_info, refSeq, m.l_haps_pair, scf, self.outfiles.fo_hap_pair_sam,
    #                                self.outfiles.fo_hap_pair_fasta, scf_title)
    #
    #         per_PE2file = datetime.datetime.now() - t1
    #
    #         print scf + '\treadBAM\t' + str(per_single_readBAM) + '\tRanbow_single\t' + str(per_single_run)\
    #               + '\tsingle2file\t' + str(per_single_run2file) + '\tpair_data\t' + str(per_dataPE)\
    #               + '\tranbow_pair\t' + str(per_PE) + '\tpair_2file\t' + str(per_PE2file)
    #

    def haplotyping_change_for_ranbow_new(self, l):
        #print  '___________________', 'in go ', '____________', self.max_type_len
        fastafile = pysam.FastaFile(self.ref)
        samfile = pysam.AlignmentFile(self.bam, 'rb')
        mp_vcf = self.load_inx_file(self.vcf_inx)
        for scfinfo in l:
            print scfinfo
            scf = scfinfo[0]
            scf_start = scfinfo[1]
            scf_end = scfinfo[2]
            scf_title = scfinfo[3]

            refSeq1 = fastafile.fetch(scf)
            refSeq = refSeq1.lower().replace("n","N")

            
            if scf not in mp_vcf.keys():
                if self.is_eval:
                    continue

                t_time0 = datetime.datetime.now()
                self.hap2fasta_and_sam_no_SNP(refSeq, scf, self.outfiles.fo_hap_single_sam, self.ploidy)
                t_time1 = datetime. datetime.now()
                print 'No SNP'
                print '>>>\t' + scf + '\treadBAM\t0:00:00.0\tRanbow_single\t0:00:00.0\tsingle2file\t0:00:00.0\tpair_data\t0:00:00.0\tranbow_pair\t0:00:00.0\tpair_2file\t0:00:00.0\n'
                continue

            noR2hits = noPairR2hits = HIP_all = noR1hit = ReadWith1Hit = noR1plusHits = 0

            self.vcf.seek(mp_vcf[scf])
            vcf_info = self.indexingVCF(scf, self.vcf)
            scf_hap_len = len(vcf_info[0])

            l_frag_single = []
            l_frag_pair = []
            l_seg_pair = []
            l_seg = []
            l_seg_single = []
            l_seg_454 = []
            t1 = datetime.datetime.now()

            cnt1 = 0



            if scf in samfile.references:
                for read in samfile.fetch(scf, scf_start - 1, scf_end + 1):
                    if (read.template_length > 35000 or read.template_length < -35000):
                        #print read.template_length
                        continue



                    if read.cigarstring != None and  read.query_sequence != None:
                        lfs = '-'
                        lfs = 's' if not read.is_paired else 'r' if read.is_reverse else 'l'
                        r = self.readd(pos=read.pos + 1, cigar=read.cigarstring, read=read.query_sequence, name=read.query_name,
                                       leftRightSingle=lfs,qual=read.query_qualities)

                        mp_ReadSegments, mp_ReadSegments_qual= self.make_segments_for_read(r, vcf_info)

                        if read.is_paired and len(mp_ReadSegments) >= 1 and read.next_reference_name == read.reference_name :
                            l_seg += [(read.query_name, mp_ReadSegments,mp_ReadSegments_qual)]





                        if read.is_paired and len(mp_ReadSegments) >= 2 and read.next_reference_name != read.reference_name :
                            l_seg_single += [self.segment2tuple(mp_ReadSegments,mp_ReadSegments_qual)]

                        if not read.is_paired  and len(mp_ReadSegments) >= 2 and not self.is_eval:
                            l_seg_single += [self.segment2tuple(mp_ReadSegments,mp_ReadSegments_qual)]

                        # if len(mp_ReadSegments) >= 1 and read.is_paired:
                        #     l_seg_pair += [(read.query_name,mp_ReadSegments)]
                        #
                        # if len(mp_ReadSegments) >= 2:
                        #     l_frag_single += [self.segment2frag(mp_ReadSegments)]


                        cnt1 += 1

                # l_frag_single_uniq = self.convert_to_uniq(l_frag_single)
                # del l_frag_single
            l_seg_pair_merged = self.convert_to_uniq_pairSEG_from_seg(l_seg, l_seg_single)





            del l_seg
            del l_seg_single





            per_single_readBAM = datetime.datetime.now() - t1
            print 'Run ranbow - single mode'
            t1 = datetime.datetime.now()

            import methods
            m = methods.methods()
            m.algName = self.algType
            m.algorithm_mode = 'single'
            m.uniq_frag_mode = 'frag'
            m.ploidy = self.ploidy
            m.max_type_len = self.max_type_len
            #m.single_frag = l_frag_single_uniq
            m.single_frag = l_seg_pair_merged

            m.pair_frag = []
            m.scf_hap_len = scf_hap_len
            # import cProfile
            # cProfile.runctx('m.run_method_()', None, locals())
            m.run_method_()


            per_single_run = datetime.datetime.now() - t1
            t1 = datetime.datetime.now()



            self.refine_haplotype_length(m.l_haps_single, scfName=scf, folderIndex=self.outfiles.folder_index)



            self.write_subhaps_to_file(m.l_haps_single, self.outfiles.fo_hap_single, scf, scf_hap_len)
            self.hap2fasta_and_sam(vcf_info, refSeq, m.l_haps_single, scf, self.outfiles.fo_hap_single_sam,
                                   self.outfiles.fo_hap_single_fasta, scf_title)
            per_single_run2file = datetime.datetime.now() - t1
                #del l_frag_single_uniq

            if False:
                print 'Read data for Paired-end mode'
                t1 = datetime.datetime.now()

                m.algorithm_mode = 'pair'
                m.uniq_frag_mode = 'seg'
                if m.uniq_frag_mode == 'seg':
                    m.l_frag_pair_uniq = self.convert_to_uniq_pairSEG_from_seg(l_seg_pair)
                else:
                    m.l_frag_pair_uniq = self.convert_to_uniq_pair_from_seg(l_seg_pair)
                del l_seg_pair

                per_dataPE = datetime.datetime.now() - t1
                print 'Run Ranbow - paired mode'
                t1 = datetime.datetime.now()

                m.run_method_()

                per_PE = datetime.datetime.now() - t1
                t1 = datetime.datetime.now()

                self.write_subhaps_to_file(m.l_haps_pair, self.outfiles.fo_hap_pair, scf, scf_hap_len)
                self.hap2fasta_and_sam(vcf_info, refSeq, m.l_haps_pair, scf, self.outfiles.fo_hap_pair_sam,
                                       self.outfiles.fo_hap_pair_fasta, scf_title)

                per_PE2file = datetime.datetime.now() - t1

                print scf + '\treadBAM\t' + str(per_single_readBAM) + '\tRanbow_single\t' + str(per_single_run)\
                      + '\tsingle2file\t' + str(per_single_run2file) + '\tpair_data\t' + str(per_dataPE)\
                      + '\tranbow_pair\t' + str(per_PE) + '\tpair_2file\t' + str(per_PE2file)
            else:

                print scf + '\treadBAM\t' + str(per_single_readBAM) + '\tRanbow_single\t' + str(per_single_run)\
                      + '\tsingle2file\t' + str(per_single_run2file)

    def evaluate_haplotypes(self, l):

        fastafile = pysam.FastaFile(self.ref)
        samfile = pysam.AlignmentFile(self.bam, 'rb')
        mp_vcf = self.load_inx_file(self.vcf_inx)


        for scfinfo in l:
            scf = scfinfo[0]
            if scf not in self.mpHapSing:
                continue
            print scfinfo
            scf_start = scfinfo[1]
            scf_end = scfinfo[2]
            scf_title = scfinfo[3]

            refSeq1 = fastafile.fetch(scf)
            refSeq = refSeq1.lower().replace("n", "N")

            if scf not in mp_vcf.keys():
                if self.is_eval:
                    continue

                t_time0 = datetime.datetime.now()
                self.hap2fasta_and_sam_no_SNP(refSeq, scf, self.outfiles.fo_hap_single_sam, self.ploidy)
                t_time1 = datetime.datetime.now()
                print 'No SNP'
                print '>>>\t' + scf + '\treadBAM\t0:00:00.0\tRanbow_single\t0:00:00.0\tsingle2file\t0:00:00.0\tpair_data\t0:00:00.0\tranbow_pair\t0:00:00.0\tpair_2file\t0:00:00.0\n'
                continue

            noR2hits = noPairR2hits = HIP_all = noR1hit = ReadWith1Hit = noR1plusHits = 0

            self.vcf.seek(mp_vcf[scf])
            vcf_info = self.indexingVCF(scf, self.vcf)
            scf_hap_len = len(vcf_info[0])

            l_frag_single = []
            l_frag_pair = []
            l_seg_pair = []
            l_seg = []
            l_seg_single = []
            l_seg_454 = []
            t1 = datetime.datetime.now()

            cnt1 = 0
            new_ref_list_no_pipe = []
            for scaffolds in samfile.references:
                new_ref_list_no_pipe += [scaffolds.replace ("|","")]

            if scf in new_ref_list_no_pipe:
                for read in samfile.fetch(scf, scf_start - 1, scf_end + 1):

                    if read.cigarstring != None and read.query_sequence !=None:
                        lfs = '-'
                        lfs = 's' if not read.is_paired else 'r' if read.is_reverse else 'l'
                        r = self.readd(pos=read.pos + 1, cigar=read.cigarstring, read=read.query_sequence,
                                       name=read.query_name,
                                       leftRightSingle=lfs, qual=read.query_qualities)

                        mp_ReadSegments, mp_ReadSegments_qual = self.make_segments_for_read(r, vcf_info)

                        if read.is_paired and len(
                                mp_ReadSegments) >= 1 and read.next_reference_name == read.reference_name:
                            l_seg += [(read.query_name, mp_ReadSegments, mp_ReadSegments_qual)]

                        if read.is_paired and len(
                                mp_ReadSegments) >= 2 and read.next_reference_name != read.reference_name:
                            l_seg_single += [self.segment2tuple(mp_ReadSegments, mp_ReadSegments_qual)]

                        if not read.is_paired and len(mp_ReadSegments) >= 2 and not self.is_eval:
                            l_seg_single += [self.segment2tuple(mp_ReadSegments, mp_ReadSegments_qual)]

                        # if len(mp_ReadSegments) >= 1 and read.is_paired:
                        #     l_seg_pair += [(read.query_name,mp_ReadSegments)]
                        #
                        # if len(mp_ReadSegments) >= 2:
                        #     l_frag_single += [self.segment2frag(mp_ReadSegments)]

                        if self.is_eval and len(mp_ReadSegments) >= 1:
                            l_seg_454 += [(read.query_name, mp_ReadSegments, min(mp_ReadSegments.keys()),
                                           max(mp_ReadSegments.keys()) + 1)]

                        cnt1 += 1

                        # l_frag_single_uniq = self.convert_to_uniq(l_frag_single)
                        # del l_frag_single
            l_seg_pair_merged = self.convert_to_uniq_pairSEG_from_seg(l_seg, l_seg_single)

            del l_seg
            del l_seg_single

            if self.is_eval:
                if len(l_seg_454) > 0:
                    self.evaluate_haplotype_not_block(l_seg_454, scf, self.mpHapSing,
                                                      self.outfiles.fo_eval_singel, self.hapSing)
                    # self.evaluate_haplotype_blocks(l_frag_single_uniq, scf, self.mpHapPair, self.outfiles.fo_eval_pair, self.hapPair)
                    # self.evaluate_haplotype_not_block(l_frag_single_uniq, scf, self.mpHapPair, self.outfiles.fo_eval_pair_NoBlock, self.hapPair)
                    # self.evaluate_haplotype_not_block(l_seg_454, scf, self.mpHapPair,
                    #                                   self.outfiles.fo_eval_pair, self.hapPair)
                continue

    def addSMRT(self, l):

        fastafile = pysam.FastaFile(self.ref)
        samfile = pysam.AlignmentFile(self.bam, 'rb')
        mp_vcf = self.load_inx_file(self.vcf_inx)


        for scfinfo in l:




            scf = scfinfo[0]
            if scf not in self.mpHapSing:
                continue
            print scfinfo




            scf_start = scfinfo[1]
            scf_end = scfinfo[2]
            scf_title = scfinfo[3]

            refSeq1 = fastafile.fetch(scf)
            refSeq = refSeq1.lower().replace("n", "N")

            if scf not in mp_vcf.keys():
                continue


            noR2hits = noPairR2hits = HIP_all = noR1hit = ReadWith1Hit = noR1plusHits = 0

            self.vcf.seek(mp_vcf[scf])
            vcf_info = self.indexingVCF(scf, self.vcf)
            scf_hap_len = len(vcf_info[0])

            l_frag_single = []
            l_frag_pair = []
            l_seg_pair = []
            l_seg = []
            l_seg_single = []
            l_seg_454 = []
            t1 = datetime.datetime.now()

            cnt1 = 0

            # new_ref_list_no_pipe = []
            # for scf in samfile.references:
            #     new_ref_list_no_pipe += [scf.replace ("|","")]

            if scf in samfile.references:
                for read in samfile.fetch(scf, scf_start - 1, scf_end + 1):

                    if read.cigarstring != None:
                        lfs = '-'
                        lfs = 's' if not read.is_paired else 'r' if read.is_reverse else 'l'
                        r = self.readd(pos=read.pos + 1, cigar=read.cigarstring, read=read.query_sequence,
                                       name=read.query_name,
                                       leftRightSingle=lfs, qual=read.query_qualities)

                        mp_ReadSegments, mp_ReadSegments_qual = self.make_segments_for_read(r, vcf_info)

                        if read.is_paired and len(
                                mp_ReadSegments) >= 1 and read.next_reference_name == read.reference_name:
                            l_seg += [(read.query_name, mp_ReadSegments, mp_ReadSegments_qual)]

                        if read.is_paired and len(
                                mp_ReadSegments) >= 2 and read.next_reference_name != read.reference_name:
                            l_seg_single += [self.segment2tuple(mp_ReadSegments, mp_ReadSegments_qual)]



                        # if len(mp_ReadSegments) >= 1 and read.is_paired:
                        #     l_seg_pair += [(read.query_name,mp_ReadSegments)]
                        #
                        # if len(mp_ReadSegments) >= 2:
                        #     l_frag_single += [self.segment2frag(mp_ReadSegments)]

                        if  len(mp_ReadSegments) >= 1:
                            l_seg_454 += [(read.query_name, mp_ReadSegments, min(mp_ReadSegments.keys()),
                                           max(mp_ReadSegments.keys()) + 1)]

                        cnt1 += 1

                        # l_frag_single_uniq = self.convert_to_uniq(l_frag_single)
                        # del l_frag_single
            l_seg_pair_merged = self.convert_to_uniq_pairSEG_from_seg(l_seg, l_seg_single)

            del l_seg
            del l_seg_single

            if len(l_seg_454) > 0:
                self.merge_hap_and_nano(l_seg_454, scf, self.mpHapSing,
                                                  self.outfiles.fo_addSMRT_singel, self.hapSing)
                # self.evaluate_haplotype_blocks(l_frag_single_uniq, scf, self.mpHapPair, self.outfiles.fo_eval_pair, self.hapPair)
                # self.evaluate_haplotype_not_block(l_frag_single_uniq, scf, self.mpHapPair, self.outfiles.fo_eval_pair_NoBlock, self.hapPair)
                # self.evaluate_haplotype_not_block(l_seg_454, scf, self.mpHapPair,
                #                                   self.outfiles.fo_eval_pair, self.hapPair)
            continue

    def go(self):
        #print '___________________', 'in go ', '____________', self.is_eval,self.is_addSMRT

        cpu_id = self.processorNumber
        l_l_selected_scf = self.distribute_selected_scaffolds()
        if self.is_eval:
            self.evaluate_haplotypes(l_l_selected_scf[cpu_id])
        elif self.is_addSMRT:
            self.addSMRT(l_l_selected_scf[cpu_id])
        else:
            #print '___________________', 'in go ', '____________', self.is_eval, self.is_addSMRT
            #print self.max_type_len
            self.haplotyping_change_for_ranbow_new(l_l_selected_scf[cpu_id])
