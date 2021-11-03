###---### mahale behine kardan
import ast
from select import select
from threading import _enumerate

import numpy as np
#from Cython.Shadow import typeof
#from jcvi.projects.napus import ploidy
#from jsonschema.validators import extend
#from numpy.f2py.auxfuncs import throw_error
#from scipy.signal.lti_conversion import cont2discrete

import haplotyper
#import networkx as nx
import sys
import itertools
import operator
#from methods import _subhap_seg


class _subhap_seg:
    def __init__(self, haplotype, supporting_fragments, flag, corrected_error, name, type, block_no,
                 frag_list_to_merge_later):
        self.name = name
        self.block_no = []
        self.type = tuple(type)
        self.haplotype = haplotype.copy()
        self.start_pos = min(haplotype.keys())
        self.end_pos = max(haplotype.keys()) + 1
        self.length = self.end_pos - self.start_pos
        self.flag = flag
        self.corrected_error = corrected_error
        self.supporting_fragments = list(supporting_fragments)
        self.frag_list_to_merge_later = list(frag_list_to_merge_later)
        self.mask = ''
        self.connect_to = {}
        self.connect_to_backward = {}
        self.dad = -1
        self.me = -1
        self.my_initial_blocks = ""


    # def __add__(self, other):
    #     ans = self.copy()
    #     ans.name = ans.name + '*' + other.name
    #     ans.block_no = min(ans.block_no, other.block_no)
    #     ans.type = tuple(sorted(ans.type + other.type))
    #
    #
    #     corrected_error = 0
    #     for pos in other.haplotype:
    #         new_allele, new_copies = other.haplotype[pos]
    #         if pos in ans.haplotype:
    #             old_allele, old_copies = ans.haplotype[pos]
    #             if new_allele == old_allele:
    #                 ans.haplotype[pos] = (old_allele, old_copies + new_copies)
    #             else:
    #                 if old_copies >= new_copies:
    #                     ans.haplotype[pos] = (old_allele, old_copies )
    #                     corrected_error += new_copies
    #                 else:
    #                     ans.haplotype[pos] = (new_allele,  new_copies)
    #                     corrected_error += old_copies
    #         else:
    #             ans.haplotype[pos] = (new_allele, new_copies)
    #
    #
    #     ans.start_pos = min(ans.haplotype.keys())
    #     ans.end_pos = max(ans.haplotype.keys()) + 1
    #     ans.length =  ans.end_pos - ans.start_pos
    #     ans.supporting_fragments += other.supporting_fragments
    #     ans.frag_list_to_merge_later += other.frag_list_to_merge_later

        # return ans

    def empty_me(self):
        self.haplotype = -1
        self.start_pos = -1
        self.end_pos = -1
        self.length = -1
        self.flag = -1
        self.corrected_error = -1
        self.supporting_fragments = -1
        self.frag_list_to_merge_later = -1
        self.mask = -1
        self.connect_to = -1
        self.connect_to_backward = -1
        self.type = -1

    def __iadd__(self, other):
        self.name = self.name + '__' + other.name
        self.block_no += other.block_no
        # print self.type , other.type
        # self.print_subhaps_full()
        # other.print_subhaps_full()

        self.type = self.type + tuple([-1]) + other.type


        corrected_error = 0
        for pos in other.haplotype:
            new_allele, new_copies = other.haplotype[pos]
            if pos in self.haplotype:
                old_allele, old_copies = self.haplotype[pos]
                if new_allele == old_allele:
                    self.haplotype[pos] = (old_allele, old_copies + new_copies)
                else:
                    if old_copies >= new_copies:
                        self.haplotype[pos] = (old_allele, old_copies )
                        corrected_error += new_copies
                    else:
                        self.haplotype[pos] = (new_allele,  new_copies)
                        corrected_error += old_copies
            else:
                self.haplotype[pos] = (new_allele, new_copies)


        self.start_pos = min(self.haplotype.keys())
        self.end_pos = max(self.haplotype.keys()) + 1
        self.length =  self.end_pos - self.start_pos
        self.supporting_fragments += other.supporting_fragments
        self.frag_list_to_merge_later += other.frag_list_to_merge_later

        for i in other.connect_to:
            if i in self.connect_to:
                self.connect_to[i] += other.connect_to[i]
            else:
                self.connect_to[i] = other.connect_to[i]


        for i in other.connect_to_backward:
            blc, hap = i
            # if len(other.block_no) > 1:
            #     print 'Test: the other block should not have more than one block_no because all nodes are gathering to the front'
            #     print other.print_subhaps_full()
            #     exit()
            if blc < other.block_no[0] and blc not in self.block_no:
                self.connect_to[i]  = self.connect_to.get(i,0) + other.connect_to_backward[i]


        other.dad = self.me

        other.empty_me()
        return self




    def __lt__(self, other):
        return self.start_pos < other.start_pos

    def copy_from(self, sp):
        self.name = sp.name
        self.block_no = sp.block_no
        self.type = tuple(sp.type)
        self.haplotype = sp.haplotype.copy()
        self.start_pos = sp.start_pos
        self.length = sp.length
        self.end_pos = sp.end_pos
        self.flag = sp.flag
        self.end_pos = sp.end_pos
        self.corrected_error = sp.corrected_error
        self.supporting_fragments = list(sp.supporting_fragments)
        self.frag_list_to_merge_later = list(sp.frag_list_to_merge_later)

    def reverse(self, scffoldLength):
        l_quality = self.quality[2:].split('-')  # to delete 'Q:'
        new_quality = 'Q:' + '-'.join(reversed(l_quality))
        return _subhap(self.haplotype[::-1], scffoldLength - 1 - (self.start_pos + self.length - 1), new_quality,
                       supporting_fragments, 'subhap')

    def find_err_of_merging_subhap_of_a_type_on_supporting_frags(self, s_h, l_subhaps):  # it means that for sure there is no doubt of merging them
        mn = sys.maxint
        mx = -1
        for r in self.supporting_fragments:
            if l_subhaps[r].start_pos < mn:
                mn = l_subhaps[r].start_pos
            if l_subhaps[r].end_pos > mx:
                mx = l_subhaps[r].end_pos

        for r in s_h.supporting_fragments:
            if l_subhaps[r].start_pos < mn:
                mn = l_subhaps[r].start_pos
            if l_subhaps[r].end_pos > mx:
                mx = l_subhaps[r].end_pos

        l_consensus = []
        for ii in range(mn, mx):
            l_consensus.append({})
        for r in self.supporting_fragments:
            for ii in range(l_subhaps[r].length):
                real_pos = ii + l_subhaps[r].start_pos - mn
                if l_subhaps[r].haplotype[ii] in l_consensus[real_pos]:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] += l_subhaps[r].quality[ii]
                else:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] = l_subhaps[r].quality[ii]

        for r in s_h.supporting_fragments:
            for ii in range(l_subhaps[r].length):
                real_pos = ii + l_subhaps[r].start_pos - mn
                if l_subhaps[r].haplotype[ii] in l_consensus[real_pos]:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] += l_subhaps[r].quality[ii]
                else:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] = l_subhaps[r].quality[ii]

        c_frag = ''
        c_frag_quality = []

        error = 0
        for mp_allele_freq in l_consensus:
            frequent_allele = max(mp_allele_freq.iteritems(), key=operator.itemgetter(1))[0]
            for allele in mp_allele_freq:
                if allele != frequent_allele:
                    error += mp_allele_freq[allele]

        return error

    def merge_with(self, subhap_to_del, l_subhaps):

        mn = sys.maxint
        mx = -1
        for r in self.supporting_fragments:
            if l_subhaps[r].start_pos < mn:
                mn = l_subhaps[r].start_pos
            if l_subhaps[r].end_pos > mx:
                mx = l_subhaps[r].end_pos

        for r in subhap_to_del.supporting_fragments:
            if l_subhaps[r].start_pos < mn:
                mn = l_subhaps[r].start_pos
            if l_subhaps[r].end_pos > mx:
                mx = l_subhaps[r].end_pos

        l_consensus = []
        for ii in range(mn, mx):
            l_consensus.append({})
        for r in self.supporting_fragments:
            for ii in range(l_subhaps[r].length):
                real_pos = ii + l_subhaps[r].start_pos - mn
                if l_subhaps[r].haplotype[ii] in l_consensus[real_pos]:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] += l_subhaps[r].quality[ii]
                else:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] = l_subhaps[r].quality[ii]

        for r in subhap_to_del.supporting_fragments:
            for ii in range(l_subhaps[r].length):
                real_pos = ii + l_subhaps[r].start_pos - mn
                if l_subhaps[r].haplotype[ii] in l_consensus[real_pos]:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] += l_subhaps[r].quality[ii]
                else:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] = l_subhaps[r].quality[ii]

        c_frag = ''
        c_frag_quality = []

        error = 0
        for mp_allele_freq in l_consensus:
            frequent_allele = max(mp_allele_freq.iteritems(), key=operator.itemgetter(1))[0]
            for allele in mp_allele_freq:
                if allele != frequent_allele:
                    error += mp_allele_freq[allele]
            c_frag += frequent_allele
            c_frag_quality.append(mp_allele_freq[frequent_allele])

        self = _subhap(c_frag, mn, c_frag_quality, self.supporting_fragments + subhap_to_del.supporting_fragments, 1,
                       error)
        subhap_to_del.flag = -1

    def merge_with_frag_seg(self, fragment, inx_frag):
        # copies, start_pos, tmp_frag, end_pos = fragment
        # frag = ast.literal_eval(tmp_frag)#here
        copies, start_pos, frag, end_pos, _ = fragment


        for i in frag:
            if i in self.haplotype:
                if self.haplotype[i][0] == frag[i]:
                    self.haplotype[i] = (frag[i], self.haplotype[i][1] + copies)
            else:
                self.haplotype[i] = (frag[i], copies)


        self.start_pos = min(self.haplotype.keys())
        self.end_pos = max(self.haplotype.keys()) + 1
        self.length =  self.end_pos - self.start_pos
        self.supporting_fragments.append(inx_frag)

    def update_hap_by_stored_fragments(self):
        for f in self.frag_list_to_merge_later:
            frag = f[3]
            start_pos = int(f[2])
            copies = int(f[1])
            frag_len = len(frag)
            end_pos = frag_len + start_pos

            mn = min(self.start_pos, start_pos)
            mx = max(self.end_pos, end_pos)

            hap = []
            hap_q = []

            for i in range(mn, mx):
                hap.append('-')
                hap_q.append(0)

            for i in range(self.start_pos, self.end_pos):
                hap[i - mn] = self.haplotype[i - self.start_pos]
                hap_q[i - mn] += self.quality[i - self.start_pos]

            for i in range(start_pos, end_pos):
                if frag[i - start_pos] == '-':
                    continue
                if hap[i - mn] == '-':
                    hap[i - mn] = frag[i - start_pos]
                elif hap[i - mn] == frag[i - start_pos]:
                    hap_q[i - mn] += copies
                else:
                    if copies > hap_q[i - mn]:
                        hap[i - mn] = frag[i - start_pos]
                        hap_q[i - mn] = copies - hap_q[i - mn]
                    else:
                        hap_q[i - mn] -= copies

            self.haplotype = ''.join(hap)
            self.length = len(hap)
            self.start_pos = mn
            self.end_pos = mn + len(hap)
            self.quality = hap_q
        self.frag_list_to_merge_later = []

    def print_subhaps_compact(self):
        print str(self.start_pos) + '\t' + str(self.end_pos) + '\t' + str(self.flag) + '\t' + self.name + '\t' + str(len(self.supporting_fragments)) + '\t' + str(self.type) + '\t' + str(self.block_no)
        return

    def print_subhaps_full(self):
        print '---Begin--------'
        print 'name:  ' + self.name
        print ' block#    :   ' + str(self.block_no)
        print ' type    :   ' + str( self.type)
        print ' haplotype    :   ' + str(self.haplotype)
        print ' start    :   ' + str(self.start_pos)
        print ' end    :   ' + str(self.end_pos)
        print ' support    :   ' + str(self.supporting_fragments)
        print ' correct error    :   ' + str(self.corrected_error)
        print ' frag    :   ' + str(self.frag_list_to_merge_later)
        print ' flag    :   ' + str(self.flag)
        print ' flag    :   ' + str(self.connect_to)
        print ' flag    :   ' + str(self.connect_to_backward)
        print '---End-----------'
        return

    def print_hap_with_start(self, start):
        print (self.start_pos - start) * '-' + self.haplotype + '\t' + str(self.start_pos) + '\t' + self.name

    def print_hap_with_start_just_hap(self, start):
        print (self.start_pos - start) * '-' + self.haplotype

    def copy(self):
        return _subhap_seg(self.haplotype, self.supporting_fragments, self.flag,
                       self.corrected_error, self.name, self.type, self.block_no, self.frag_list_to_merge_later)



    def make_mask(self, scf_len):
        l = []
        for i in range(scf_len):
            l.append('-')
        for i in range(self.start_pos, self.end_pos):
            if self.haplotype[i - self.start_pos] != '-':
                l[i] = '*'
        self.mask = ''.join(l)
        return


#
# class _subhap:
#     def __init__(self, haplotype, start_pos, quality, supporting_fragments, flag, corrected_error, name, type, block_no,
#                  frag_list_to_merge_later):
#         self.name = name
#         self.block_no = block_no
#         self.type = tuple(type)
#         self.haplotype = haplotype
#         self.start_pos = start_pos
#         self.length = len(haplotype)
#         self.end_pos = self.start_pos + self.length
#         self.flag = flag
#         self.corrected_error = corrected_error
#         self.quality = list(quality)
#         self.supporting_fragments = list(supporting_fragments)
#         self.frag_list_to_merge_later = list(frag_list_to_merge_later)
#         self.mask = ''
#
#     def __add__(self, other):
#         ans = self.copy()
#         ans.name = ans.name + '*' + other.name
#         ans.block_no = min(ans.block_no, other.block_no)
#         ans.type = tuple(sorted(ans.type + other.type))
#         mn = min(ans.start_pos, other.start_pos)
#         mx = max(ans.end_pos, other.end_pos)
#         hap = []
#         qual = []
#         for i in range(mn, mx):
#             hap.append('-')
#             qual.append(0)
#
#         for i in range(ans.start_pos, ans.end_pos):
#             hap[i - mn] = ans.haplotype[i - ans.start_pos]
#             qual[i - mn] = ans.quality[i - ans.start_pos]
#         corrected_error = 0
#         for i in range(other.start_pos, other.end_pos):
#             if hap[i - mn] == '-':
#                 hap[i - mn] = other.haplotype[i - other.start_pos]
#                 qual[i - mn] = other.quality[i - other.start_pos]
#             elif hap[i - mn] == other.haplotype[i - other.start_pos]:
#                 qual[i - mn] += other.quality[i - other.start_pos]
#             else:
#                 if qual[i - mn] > other.quality[i - other.start_pos]:
#                     corrected_error = qual[i - mn] - other.quality[i - other.start_pos]
#                     qual[i - mn] -= other.quality[i - other.start_pos]
#                 else:
#                     corrected_error = qual[i - mn]
#                     qual[i - mn] = other.quality[i - other.start_pos] - qual[i - mn]
#                     hap[i - mn] = other.haplotype[i - other.start_pos]
#
#         ans.start_pos = mn
#         ans.end_pos = mx
#         ans.haplotype = ''.join(hap)
#         ans.quality = qual
#         ans.length = len(hap)
#         ans.flag = ans.flag
#         ans.corrected_error += (other.corrected_error + corrected_error)
#         ans.supporting_fragments += other.supporting_fragments
#         ans.frag_list_to_merge_later += other.frag_list_to_merge_later
#         return ans
#
#     def __lt__(self, other):
#         return self.start_pos < other.start_pos
#
#     def copy_from(self, sp):
#         self.name = sp.name
#         self.block_no = sp.block_no
#         self.type = tuple(sp.type)
#         self.haplotype = sp.haplotype
#         self.start_pos = sp.start_pos
#         self.length = sp.length
#         self.end_pos = sp.end_pos
#         self.quality = list(sp.quality)
#         self.flag = sp.flag
#         self.end_pos = sp.end_pos
#         self.corrected_error = sp.corrected_error
#         self.supporting_fragments = list(sp.supporting_fragments)
#         self.frag_list_to_merge_later = list(sp.frag_list_to_merge_later)
#
#     def reverse(self, scffoldLength):
#         l_quality = self.quality[2:].split('-')  # to delete 'Q:'
#         new_quality = 'Q:' + '-'.join(reversed(l_quality))
#         return _subhap(self.haplotype[::-1], scffoldLength - 1 - (self.start_pos + self.length - 1), new_quality,
#                        supporting_fragments, 'subhap')
#
#     def find_err_of_merging_subhap_of_a_type_on_supporting_frags(self, s_h, l_subhaps):  # it means that for sure there is no doubt of merging them
#
#         mn = sys.maxint
#         mx = -1
#         for r in self.supporting_fragments:
#             if l_subhaps[r].start_pos < mn:
#                 mn = l_subhaps[r].start_pos
#             if l_subhaps[r].end_pos > mx:
#                 mx = l_subhaps[r].end_pos
#
#         for r in s_h.supporting_fragments:
#             if l_subhaps[r].start_pos < mn:
#                 mn = l_subhaps[r].start_pos
#             if l_subhaps[r].end_pos > mx:
#                 mx = l_subhaps[r].end_pos
#
#         l_consensus = []
#         for ii in range(mn, mx):
#             l_consensus.append({})
#         for r in self.supporting_fragments:
#             for ii in range(l_subhaps[r].length):
#                 real_pos = ii + l_subhaps[r].start_pos - mn
#                 if l_subhaps[r].haplotype[ii] in l_consensus[real_pos]:
#                     l_consensus[real_pos][l_subhaps[r].haplotype[ii]] += l_subhaps[r].quality[ii]
#                 else:
#                     l_consensus[real_pos][l_subhaps[r].haplotype[ii]] = l_subhaps[r].quality[ii]
#
#         for r in s_h.supporting_fragments:
#             for ii in range(l_subhaps[r].length):
#                 real_pos = ii + l_subhaps[r].start_pos - mn
#                 if l_subhaps[r].haplotype[ii] in l_consensus[real_pos]:
#                     l_consensus[real_pos][l_subhaps[r].haplotype[ii]] += l_subhaps[r].quality[ii]
#                 else:
#                     l_consensus[real_pos][l_subhaps[r].haplotype[ii]] = l_subhaps[r].quality[ii]
#
#         c_frag = ''
#         c_frag_quality = []
#
#         error = 0
#         for mp_allele_freq in l_consensus:
#             frequent_allele = max(mp_allele_freq.iteritems(), key=operator.itemgetter(1))[0]
#             for allele in mp_allele_freq:
#                 if allele != frequent_allele:
#                     error += mp_allele_freq[allele]
#
#         return error
#
#     def merge_with(self, subhap_to_del, l_subhaps):
#
#         mn = sys.maxint
#         mx = -1
#         for r in self.supporting_fragments:
#             if l_subhaps[r].start_pos < mn:
#                 mn = l_subhaps[r].start_pos
#             if l_subhaps[r].end_pos > mx:
#                 mx = l_subhaps[r].end_pos
#
#         for r in subhap_to_del.supporting_fragments:
#             if l_subhaps[r].start_pos < mn:
#                 mn = l_subhaps[r].start_pos
#             if l_subhaps[r].end_pos > mx:
#                 mx = l_subhaps[r].end_pos
#
#         l_consensus = []
#         for ii in range(mn, mx):
#             l_consensus.append({})
#         for r in self.supporting_fragments:
#             for ii in range(l_subhaps[r].length):
#                 real_pos = ii + l_subhaps[r].start_pos - mn
#                 if l_subhaps[r].haplotype[ii] in l_consensus[real_pos]:
#                     l_consensus[real_pos][l_subhaps[r].haplotype[ii]] += l_subhaps[r].quality[ii]
#                 else:
#                     l_consensus[real_pos][l_subhaps[r].haplotype[ii]] = l_subhaps[r].quality[ii]
#
#         for r in subhap_to_del.supporting_fragments:
#             for ii in range(l_subhaps[r].length):
#                 real_pos = ii + l_subhaps[r].start_pos - mn
#                 if l_subhaps[r].haplotype[ii] in l_consensus[real_pos]:
#                     l_consensus[real_pos][l_subhaps[r].haplotype[ii]] += l_subhaps[r].quality[ii]
#                 else:
#                     l_consensus[real_pos][l_subhaps[r].haplotype[ii]] = l_subhaps[r].quality[ii]
#
#         c_frag = ''
#         c_frag_quality = []
#
#         error = 0
#         for mp_allele_freq in l_consensus:
#             frequent_allele = max(mp_allele_freq.iteritems(), key=operator.itemgetter(1))[0]
#             for allele in mp_allele_freq:
#                 if allele != frequent_allele:
#                     error += mp_allele_freq[allele]
#             c_frag += frequent_allele
#             c_frag_quality.append(mp_allele_freq[frequent_allele])
#
#         self = _subhap(c_frag, mn, c_frag_quality, self.supporting_fragments + subhap_to_del.supporting_fragments, 1,
#                        error)
#         subhap_to_del.flag = -1
#
#     def merge_with_frag(self, fragment, inx_frag):
#         frag = fragment[2]
#         start_pos = int(fragment[1])
#         copies = int(fragment[0])
#         frag_len = len(frag)
#         end_pos = frag_len + start_pos
#         mn = min(self.start_pos, start_pos)
#         mx = max(self.end_pos, end_pos)
#
#         hap = []
#         hap_q = []
#         for i in range(mn, mx):
#             hap.append('-')
#             hap_q.append(0)
#
#         for i in range(self.start_pos, self.end_pos):
#             hap[i - mn] = self.haplotype[i - self.start_pos]
#             hap_q[i - mn] += self.quality[i - self.start_pos]
#
#         for i in range(start_pos, end_pos):
#             if hap[i - mn] == '-':
#                 hap[i - mn] = frag[i - start_pos]
#             hap_q[i - mn] += copies
#
#         self.haplotype = ''.join(hap)
#         self.length = len(hap)
#         self.start_pos = mn
#         self.end_pos = mn + len(hap)
#         self.quality = hap_q
#         self.supporting_fragments.append(inx_frag)
#
#     def update_hap_by_stored_fragments(self):
#         for f in self.frag_list_to_merge_later:
#             frag = f[3]
#             start_pos = int(f[2])
#             copies = int(f[1])
#             frag_len = len(frag)
#             end_pos = frag_len + start_pos
#
#             mn = min(self.start_pos, start_pos)
#             mx = max(self.end_pos, end_pos)
#
#             hap = []
#             hap_q = []
#
#             for i in range(mn, mx):
#                 hap.append('-')
#                 hap_q.append(0)
#
#             for i in range(self.start_pos, self.end_pos):
#                 hap[i - mn] = self.haplotype[i - self.start_pos]
#                 hap_q[i - mn] += self.quality[i - self.start_pos]
#
#             for i in range(start_pos, end_pos):
#                 if frag[i - start_pos] == '-':
#                     continue
#                 if hap[i - mn] == '-':
#                     hap[i - mn] = frag[i - start_pos]
#                 elif hap[i - mn] == frag[i - start_pos]:
#                     hap_q[i - mn] += copies
#                 else:
#                     if copies > hap_q[i - mn]:
#                         hap[i - mn] = frag[i - start_pos]
#                         hap_q[i - mn] = copies - hap_q[i - mn]
#                     else:
#                         hap_q[i - mn] -= copies
#
#             self.haplotype = ''.join(hap)
#             self.length = len(hap)
#             self.start_pos = mn
#             self.end_pos = mn + len(hap)
#             self.quality = hap_q
#         self.frag_list_to_merge_later = []
#
#     def print_subhaps_compact(self):
#         print str(self.start_pos) + '\t' + self.haplotype + '\t' + str(self.length) + '\t' + self.name
#         return
#
#     def print_subhaps_full(self):
#         print '---Begin--------'
#         print 'name:  ' + self.name
#         print ' block#    :   ' + str(self.block_no)
#         print ' type    :   ' + str(str(self.type))
#         print ' haplotype    :   ' + str(self.haplotype)
#         print ' start    :   ' + str(self.start_pos)
#         print ' end    :   ' + str(self.end_pos)
#         print ' qual    :   ' + str(str(self.quality))
#         print ' support    :   ' + str(self.supporting_fragments)
#         print ' sum_qual    :   ' + str(np.sum(self.quality))
#         print ' correct error    :   ' + str(self.corrected_error)
#         print ' frag    :   ' + str(self.frag_list_to_merge_later)
#         print '---End-----------'
#         return
#
#     def print_hap_with_start(self, start):
#         print (self.start_pos - start) * '-' + self.haplotype + '\t' + str(self.start_pos) + '\t' + self.name
#
#     def print_hap_with_start_just_hap(self, start):
#         print (self.start_pos - start) * '-' + self.haplotype
#
#     def copy(self):
#         return _subhap(self.haplotype, self.start_pos, self.quality, self.supporting_fragments, self.flag,
#                        self.corrected_error, self.name, self.type, self.block_no, self.frag_list_to_merge_later)
#
#     def make_mask(self, scf_len):
#         l = []
#         for i in range(scf_len):
#             l.append('-')
#         for i in range(self.start_pos, self.end_pos):
#             if self.haplotype[i - self.start_pos] != '-':
#                 l[i] = '*'
#         self.mask = ''.join(l)
#         return
#
#
class _type:
    def __init__(self, length, sequence, positions, quality, supporting_fragments):
        self.sequence = sequence
        self.length = length
        self.quality = quality
        self.supporting_fragments = []
        for i in supporting_fragments:
            self.supporting_fragments.append(i)

        self.positions = []
        for i in positions:
            self.positions.append(i)

    def add_supporting_fragment(self, fragment_number, quality):
        self.quality += quality
        self.supporting_fragments.append(fragment_number)

    def make_subhap_form_supporting_fragments_for_a_type(self, l_subhaps):

        mn = sys.maxint
        mx = -1

        for r in self.subhap.supporting_fragments:
            if l_subhaps[r].start_pos < mn:
                mn = l_subhaps[r].start_pos
            if l_subhaps[r].end_pos > mx:
                mx = l_subhaps[r].end_pos

        l_consensus = []
        for ii in range(mn, mx):
            l_consensus.append({})
        for r in self.subhap.supporting_fragments:

            for ii in range(l_subhaps[r].length):
                real_pos = ii + l_subhaps[r].start_pos - mn
                if l_subhaps[r].haplotype[ii] in l_consensus[real_pos]:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] += l_subhaps[r].quality[ii]
                else:
                    l_consensus[real_pos][l_subhaps[r].haplotype[ii]] = l_subhaps[r].quality[ii]

        c_frag = ''
        c_frag_quality = []

        error = 0
        for mp_allele_freq in l_consensus:
            frequent_allele = max(mp_allele_freq.iteritems(), key=operator.itemgetter(1))[0]
            for allele in mp_allele_freq:
                if allele != frequent_allele:
                    error += mp_allele_freq[allele]
            c_frag += frequent_allele
            c_frag_quality.append(mp_allele_freq[frequent_allele])

        self.subhap = _subhap(c_frag, mn, c_frag_quality, self.subhap.supporting_fragments, 1, error, [])

    def merge_with_type(self, type_to_del, l_subhaps):
        self.subhap.merge_with(type_to_del.subhap, l_subhaps)

    def print_type(self):
        print '*****'
        print self.sequence
        print str(self.quality)
        print str(self.supporting_fragments)
        for j in self.supporting_fragments:
            print str(l_frag[j])


class _fragment:
    def __init__(self, fragment, start_pos, repeat):
        self.fragment = fragment
        self.start_pos = start_pos
        self.length = len(fragment)
        self.end_pos = len(fragment) + start_pos
        self.repeat = repeat

    def __lt__(self, other):
        return self.start_pos < other.start_pos

    def reverse(self, scffoldLength):
        return _fragment(self.fragment[::-1], scffoldLength - 1 - (self.start_pos + self.length - 1), self.repeat)

    def pprint(self):
        print self.fragment,
        print self.start_pos,
        print self.end_pos,
        print


class _convert:
    def __init__(self):
        pass
    def frag2seq(self, frag, start_pos):
        mp = {}
        for i,snp in enumerate(frag):
            if snp != '-':
               mp[i+start_pos] = snp
        return mp



class methods:
    def __init__(self):
        self.algName = ''
        self.uniq_frag_mode = ''
        self.typeExtension = False
        self.bo_connect_adjacent = False
        self.ploidy = -1
        self.max_type_len = -1
        self.WinLen = -1
        self.single_frags = []
        self.pair_frags = []
        self.scf_hap_len = -1

    def load_frags(self, l_frag, l_frag_used, fi):  # phase_by_sorted_type

        line = fi.readline()
        fragNo = 0
        while line != '' and line.rstrip()[-3:] != '<<<':
            a = line.split()
            l_frag.append((int(a[0]), int(a[1]), a[2]))
            l_frag_used.append(False)
            fragNo += 1
            line = fi.readline()
        return line

    def obtain_types_nolimit(self, l_frag, typeLenMin, typeLenMax, ploidy):  # phase_by_sorted_type


        mp = []
        for i in range(typeLenMin, typeLenMax):
            mp.append({})

        for fragment in l_frag:
            # (copies, start_pos, frag_mp, end_pos) = fragment
            #
            # frag = ast.literal_eval(frag_mp)#here
            (copies, start_pos, frag, end_pos) = fragment
            l = frag.keys()
            l.sort()

            # (copies, start_pos, frag_mp) = fragment
            # frag_len = len(frag_)
            # l = []
            # for i in range(frag_len):
            #     if frag[i] != '-':
            #         l.append(i)



            for t_len in range(typeLenMin, typeLenMax):
                if len(l) < t_len:
                    continue

                mp_inx = t_len - 2

                for p in itertools.combinations(l, t_len):  # f as frag
                    type = ''
                    for p_i in p:
                        type += str(frag[p_i])
                    tup = tuple(p)

                    if tup in mp[mp_inx]:
                        if type in mp[mp_inx][tup]:
                            mp[mp_inx][tup][type] += copies
                        else:
                            mp[mp_inx][tup][type] = copies
                    else:
                        mp[mp_inx][tup] = {type: copies}

        return self.extract_informative_types(mp, ploidy)



    def obtain_types(self, l_frag, typeLenMin, typeLenMax, ploidy,win_type_limit_size):  # phase_by_sorted_type
        
        # win_type_limit_size = 8 submitted 

        mp = []
        for i in range(typeLenMin, typeLenMax):
            mp.append({})

        for fragment in l_frag:
            # (copies, start_pos, frag_mp, end_pos) = fragment
            #
            # frag = ast.literal_eval(frag_mp)#here
            (copies, start_pos, frag, end_pos,frag_q) = fragment
            l = frag.keys()
            l.sort()

            # (copies, start_pos, frag_mp) = fragment
            # frag_len = len(frag_)
            # l = []
            # for i in range(frag_len):
            #     if frag[i] != '-':
            #         l.append(i)



            for t_len in range(typeLenMin, typeLenMax):
                if len(l) < t_len:
                    continue
                if len(l) > win_type_limit_size:
                    all_types = []
                    for lower_bound in range(len(l) - win_type_limit_size + 1):
                        for p_1 in itertools.combinations(l[lower_bound + 1 :lower_bound + win_type_limit_size], t_len - 1):
                            all_types += [tuple([l[lower_bound]]+list(p_1))]
                    for p in itertools.combinations(l[len(l) - win_type_limit_size + 1:], t_len ):
                        all_types += [p]
                else:
                    all_types = itertools.combinations(l, t_len)

                mp_inx = t_len - 2
                for p in all_types:  # f as frag
                    type = ''
                    typeq = []
                    for p_i in p:
                        type += str(frag[p_i])
                        typeq += [frag_q[p_i]]
                    tup = tuple(p)
                    if tup in mp[mp_inx]:
                        if type in mp[mp_inx][tup]:

                            #mp[mp_inx][tup][type] += copies
                            mp[mp_inx][tup][type] =tuple([xx + yy for xx , yy in zip(typeq,mp[mp_inx][tup][type])])
                        else:
                            mp[mp_inx][tup][type] = tuple(typeq)
                    else:
                        mp[mp_inx][tup] = {type: tuple(typeq)}



        return self.extract_informative_types(mp, ploidy)



    def extract_informative_types_hight_JUSTploidy(self, mp, ploidy):

        l_all_types_inx = []
        l_all_types_info = []
        cnt = []

        for i in range(ploidy + 1):
            l_all_types_inx.append([])
            l_all_types_info.append([])
            cnt.append(0)

        for m in mp:
            for t in m:
                height_type = len(m[t])
                if height_type >= ploidy:
                    sixth_freq = sorted(m[t].values(), reverse=True)[ploidy - 1]
                    l_all_types_inx[ploidy].append((len(m[t]), sixth_freq, cnt[ploidy]))
                    l_all_types_info[ploidy].append((t, m[t]))
                    cnt[ploidy] += 1

        l_all_types = []
        for j in range(ploidy, ploidy - 1, -1):
            for i in sorted(l_all_types_inx[j], key=lambda x: (x[1], x[0]), reverse=True):
                l_all_types.append(((i[0], i[1]), l_all_types_info[j][i[2]]))

        return l_all_types


    def extract_informative_types(self, mp, ploidy):
        l_all_types = []
        for m in mp:
            for t in m:
                real_h = len(m[t])
                if real_h >= 2:
                    used_h = ploidy if real_h > ploidy else real_h
                    score = sorted([min(xx) for xx in m[t].values()], reverse=True)[used_h - 1]
                    l_all_types += [((real_h, score, used_h),(t,m[t]))]

        l_all_types.sort(key=lambda x: (x[0][2], x[0][1]), reverse=True)
        return l_all_types



    def make_arr_start_from_frag(self, ans, l_frag):  # phase_by_sorted_type
        # for i in range(len(l_frag)):
        #     _ , _, seg, _ = l_frag[i]
        #     seg = ast.literal_eval(seg)
        #     for j in seg.keys():
        #         if ans[j] == -1:
        #             ans[j] = i
        #     # start_pos = l_frag[i][1]
        #     # frag_len = len(l_frag[i][2])
        #     # for j in range(frag_len):
        #     #     if ans[start_pos + j] == -1:
        #     #         ans[start_pos + j] = i
        # prev = 0
        # for i in range(len(ans)):
        #     if ans[i] == -1:
        #         ans[i] = prev
        #     else:
        #         prev = ans[i]



        for i in range(len(l_frag) - 1, -1, -1):
            st , _, _, en, _ = l_frag[i]
            for j in range(st,en):
                ans[j] = i
            # start_pos = l_frag[i][1]
            # frag_len = len(l_frag[i][2])
            # for j in range(frag_len):
            #     if ans[start_pos + j] == -1:
            #         ans[start_pos + j] = i
        prev = len(ans) - 1
        for i in range(len(ans) - 1, -1 , -1):
            if ans[i] == -1:
                ans[i] = prev
            else:
                prev = ans[i]

    def make_arr_end_to_frag(self, ans, l_frag):  # phase_by_sorted_type

        # l = []
        # for i in range(len(l_frag)):
        #     l.append((l_frag[i][1] + len(l_frag[i][2]), l_frag[i][2], i))
        l = [(i[-1],inx) for inx,i in enumerate(l_frag)]



        for i in sorted(l, key=lambda x:x[0]):
            st , _, _, en,_ = l_frag[i[1]]
            for j in range(st,en):
                ans[j] = i[1]

        prev = 0
        for i,v in enumerate(ans):
            if v == -1:
                ans[i] = prev
            else:
                prev = v

        #
        # n = len(l_frag)
        # for i in xrange(n-1, -1, -1 ):
        #     _, be, seg, en = l_frag[i]
        #     for j in xrange(be, en):
        #         if ans[j] == -1:
        #             ans [j] = i



        # for i in sorted(l, key=lambda x: (x[2]), reverse=True):
        #     for j in range(i[0] - len(i[1]), i[0]):
        #         if ans[j] == -1:
        #             ans[j] = i[2]
        #


        #
        # prev = len(ans) - 1
        # for i in range(len(ans) - 1, -1, -1):
        #     if ans[i] == -1:
        #         ans[i] = prev
        #     else:
        #         prev = ans[i]
        #

    def make_hap_from_type(self, l_frag, t, ploidy, cnt_type, cnt_hap):

        frag_merged = {}
        for f_inx in t.supporting_fragments:
            # copies, bn, frag_str , en = l_frag[f_inx]
            # frag = ast.literal_eval(frag_str)#here
            copies, bn, frag , en, _ = l_frag[f_inx]

            for pos in frag:
                allele = frag[pos]
                if pos in frag_merged:
                    frag_merged[pos][allele] = frag_merged[pos].get(allele,0)
                    frag_merged[pos][allele] += copies
                else:
                    frag_merged[pos] = {allele:copies}

        hap = {}
        for pos in frag_merged:
            mx = 0
            corrected_error = 0
            for allele in frag_merged[pos]:
                if mx < frag_merged[pos][allele]:
                    corrected_error += mx
                    mx = frag_merged[pos][allele]
                    mx_a = allele
                else:
                    corrected_error += frag_merged[pos][allele]

            hap[pos] = (mx_a, mx)


        return _subhap_seg(hap,t.supporting_fragments, ploidy, corrected_error,
                       str(cnt_type) + '_' + str(cnt_hap), t.positions, cnt_type, [])

    def make_hap_from_type_new(self, l_frag, t, ploidy, cnt_type, cnt_hap):

        type = {}
        extended_type = {}
        for allele, pos in zip(t.sequence,sorted(t.positions)):
            type [pos] = int(allele)


        for f_inx in t.supporting_fragments:
            copies, bn, frag_str , en = l_frag[f_inx]

            frag = ast.literal_eval(frag_str)#not reachable code for changing ast function

            is_frag_coveres_type = True
            for pos in type.keys():
                if pos not in frag:
                    is_frag_coveres_type = False
                    break

            if is_frag_coveres_type:
                for pos in frag:
                    allele = frag[pos]
                    if pos in extended_type:
                        extended_type[pos][allele] = extended_type[pos].get(allele,0)
                        extended_type[pos][allele] += copies
                    else:
                        extended_type[pos] = {allele:copies}

        hap = {}

        for pos in type:
            copies = extended_type[pos][type[pos]]
            extended_type[pos] = {allele:copies}

        for pos in extended_type:
            mx = 0
            corrected_error = 0
            for allele in extended_type[pos]:
                if mx < extended_type[pos][allele]:
                    corrected_error += mx
                    mx = extended_type[pos][allele]
                    mx_a = allele
                else:
                    corrected_error += extended_type[pos][allele]

            hap[pos] = (mx_a, mx)


        return _subhap_seg(hap,t.supporting_fragments, ploidy, corrected_error,
                       str(cnt_type) + '-' + str(cnt_hap), t.positions, cnt_type, [])



    def phase_a_type(self, l_frag, l_frag_used, type_info, first_frag, last_frag, skyline, cnt_type, extention_allowed,
                     ploidy_genome):  # phase_scaffold



        l_frag_used.fill(False)


        ploidy = type_info[0][2]
        no_types = type_info[0][0]
        type_pos_real = type_info[1][0]

        mp_type_support = {}
        if no_types > ploidy:
            l_tup_sorted = sorted(type_info[1][1].items(), key=operator.itemgetter(1), reverse=True)
            for i in range(ploidy):
                mp_type_support[l_tup_sorted[i][0]] = l_tup_sorted[i][1]
        else:
            mp_type_support = type_info[1][1].copy()

        inx_first_frag = first_frag[type_pos_real[0]]
        inx_last_frag = last_frag[type_pos_real[-1]]

        typeLen = len(type_pos_real)

        l_types = []

        for f_inx in range(inx_first_frag, min(inx_last_frag + 1,len(l_frag))):
        #for f_inx in range(len(l_frag)):
            if l_frag_used[f_inx]:
                continue

            # f_copies = l_frag[f_inx][0]
            # f_start_pos = l_frag[f_inx][1]
            # f_frag = l_frag[f_inx][2]
            # f_frag_len = len(l_frag[f_inx][2])
            # f_end_pos = f_start_pos + f_frag_len

            # f_copies, f_start_pos, f_frag_str, f_end_pos = l_frag[f_inx]
            # f_frag_len = f_end_pos - f_start_pos
            #
            # f_frag = ast.literal_eval(f_frag_str)#here
            f_copies, f_start_pos, f_frag, f_end_pos,_ = l_frag[f_inx]
            f_frag_len = f_end_pos - f_start_pos





            # type = ''
            # bo_farg_contains_type = True
            # for i in type_pos_real:
            #     if f_start_pos <= i < f_end_pos and f_frag[i - f_start_pos] != '-':
            #         type += f_frag[i - f_start_pos]
            #     else:
            #         bo_farg_contains_type = False
            # if not bo_farg_contains_type:
            #     continue


            type = ''
            bo_farg_contains_type = True
            for i in type_pos_real:
                if i in f_frag:
                    type += str(f_frag[i])
                else:
                    bo_farg_contains_type = False
                    break
            if not bo_farg_contains_type:
                continue



            if type in mp_type_support:
                bo_added = False
                for t in xrange(len(l_types)):
                    if l_types[t].sequence == type:
                        l_types[t].add_supporting_fragment(f_inx, f_copies)
                        l_frag_used[f_inx] = True
                        bo_added = True
                    if bo_added:
                        break
                if not bo_added:
                    l_types.append(_type(typeLen, type, list(type_pos_real), f_copies, [f_inx]))
                    l_frag_used[f_inx] = True

        # mp_tmp= {}
        # for tmp in l_types:
        #     mp_tmp[tmp.sequence] = len(tmp.supporting_fragments)
        #
        # for tmp in mp_type_support:
        #     if mp_type_support[tmp] != mp_tmp[tmp]:
        #         print mp_type_support,mp_tmp
        #         exit()
        #

        l_haps = []
        for i in range(len(l_types)):
            l_haps.append(self.make_hap_from_type(l_frag, l_types[i], len(l_types), cnt_type, i))


        if len(l_haps) < 2:
            return l_haps

        changing = True
        while changing and extention_allowed:
            changing = False

            mn = sys.maxint
            mx = -1
            for i in range(len(l_haps)):
                if mn > l_haps[-1].start_pos:
                    mn = l_haps[-1].start_pos
                if mx < l_haps[-1].end_pos:
                    mx = l_haps[-1].end_pos - 1

            for  i_frag in xrange(first_frag[mn], last_frag[mx]):

                if l_frag_used[i_frag]:
                    continue


                cnt_shared = 0
                whichHap = -1
                bo_overlapAll = True
                for i_hap in range(len(l_haps)):
                    bo_overlap, bo_conflict = self.overlap_conflict_Region_seg(l_haps[i_hap], l_frag[i_frag])
                    if not bo_overlap:
                        bo_overlapAll = False
                        break
                    if not bo_conflict:
                        cnt_shared += 1
                        whichHap = i_hap
                    if cnt_shared > 1:
                        break

                if cnt_shared == 1 and bo_overlapAll:
                    l_haps[whichHap].merge_with_frag_seg(l_frag[i_frag], i_frag)
                    l_frag_used[i_frag] = True
                    changing = True

        #if len(l_haps) < ploidy:
        #    print "possible error phased to less than ploidy --> Expected to phase the region to {}, but it was phased to {}".format(ploidy, len(l_haps))
#            exit()

        sky_tmp = np.zeros(skyline.__len__(), dtype=np.uint8)
        for i in l_haps:
            for j in i.haplotype:
                sky_tmp[j] += 1

        for i in range(skyline.__len__()):
            if sky_tmp[i] != 0:
                skyline[i] += [(ploidy, sky_tmp[i])]


        return l_haps

    def no_conflict_with_skyline(self, type, skyline, strategy):
        assert strategy in ['tail block', 'fully phased block']
        for j in type[1][0]:
            if strategy == 'fully phased block':
                if skyline[j].__len__() != 0:
                    for i in skyline[j]:
                        if i[0] <= i[1]:
                            return False
            elif strategy == 'tail block':
                if skyline[j].__len__() != 0:
                    return False

        return True


    def blocks_to_haps(self, l_b):
        l_obtaind_hap = []
        for i in l_b:
            for j in i:
                if j.flag != -1:
                    l_obtaind_hap += [j]
        return l_obtaind_hap

    def is_valid_to_merge_triangle(self,H,H1,H2):
        bo_merge = True
        bo_merge = len(set(H.block_no) & set(H1.block_no)) == 0 and bo_merge
        bo_merge = len(set(H1.block_no) & set(H2.block_no)) == 0 and bo_merge
        bo_merge = len(set(H2.block_no) & set(H.block_no)) == 0 and bo_merge
        bo_merge = (H.flag != -1 or H1.flag != -1 or H2.flag != -1) and bo_merge
        return bo_merge




    def find_3_4angular_connections_by_merge_gready(self, l_b, threshold , ploidy):

        noBlock = len(l_b)

        merged_triangle = 0
        for b1, b2 in itertools.combinations(range(len(l_b)), 2):
            len_b1 = len(l_b[b1])
            len_b2 = len(l_b[b2])
            for h1, h2 in itertools.product(range(len_b1), range(len_b2)):
                if l_b[b1][h1].flag == -1 or l_b[b2][h2].flag == -1:
                    continue
                scr = len(set(l_b[b1][h1].supporting_fragments) & set(l_b[b2][h2].supporting_fragments))
                if scr != 0:
                    l_b[b1][h1].connect_to[(b2, h2)] = scr
                    l_b[b2][h2].connect_to_backward[(b1, h1)] = scr


        for b in xrange(noBlock - 2):  # 2 for checking triangle
            l_triangle = []
            for h in xrange(ploidy):
                H = l_b[b][h]

                if H.flag == -1:
                    continue

                if len(H.connect_to) <= 1:
                    continue

                for hap1, hap2 in itertools.combinations(H.connect_to.keys(), 2):
                    b1, h1, b2, h2 = hap1 + hap2 if hap1[0] < hap2[0] else hap2 + hap1
                    H1, H2 = l_b[b1][h1], l_b[b2][h2]

                    if b1 == b2 or H1.flag == -1 or H2.flag == -1:
                        continue

                    if (b2, h2) in H1.connect_to:
                        scrs = (H.connect_to[(b1, h1)], H.connect_to[(b2, h2)], H1.connect_to[(b2, h2)])
                        l_triangle += [(b, h, b1, h1, b2, h2, min(scrs))]

            # print len(l_triangle)



            while len(l_triangle) > 0:
                l_triangle.sort( key=lambda x: x[-1])
                b, h, b1, h1, b2, h2, scr = l_triangle.pop()
                if scr > threshold:
                    H, H1, H2 = l_b[b][h], l_b[b1][h1], l_b[b2][h2]

                    if self.is_valid_to_merge_triangle(H,H1,H2):
                        # print b2, h2, b1, h1, b,h
                        if H1.flag != -1 and H2.flag != -1 :
                            H1H2_connections = set(H1.connect_to.keys()) & set(H2.connect_to.keys())
                            for b3, h3 in H1H2_connections:
                                H3 = l_b[b3][h3]
                                if H3.flag != -1 :
                                    scrs = (H1.connect_to[(b2, h2)], H1.connect_to[(b3, h3)], H2.connect_to[(b3, h3)])
                                    l_triangle += [(b, h, b1, h1, b3, h3, min(scrs))]
                                    l_triangle += [(b, h, b2, h2, b3, h3, min(scrs))]

                        if H1.flag != -1:
                            H += H1
                        else:
                            pass
                            #print 'A'
                        if H2.flag != -1:
                            H += H2
                        else:
                            pass
                            #print 'B'

                        merged_triangle += 1
                        if len(H.connect_to) > 1:

                            for hap1, hap2 in itertools.combinations(H.connect_to.keys(), 2):
                                b1, h1, b2, h2 = hap1 + hap2 if hap1[0] < hap2[0] else hap2 + hap1
                                H1, H2 = l_b[b1][h1], l_b[b2][h2]

                                if b1 == b2 or H1.flag == -1 or H2.flag == -1:
                                    continue



                                if (b2, h2) in H1.connect_to:
                                    scrs = (H.connect_to[(b1, h1)], H.connect_to[(b2, h2)], H1.connect_to[(b2, h2)])
                                    l_triangle += [(b, h, b1, h1, b2, h2, min(scrs))]

        return merged_triangle



    def find_3_4angular_justconnections_for_nanopore_not_completed(self, l_b,connections, threshold , ploidy):

        noBlock = len(l_b)

        for c1,c2 in connections:
            print c1,c2
        exit()



        for b1, b2 in itertools.combinations(range(len(l_b)), 2):
            len_b1 = len(l_b[b1])
            len_b2 = len(l_b[b2])
            for h1, h2 in itertools.product(range(len_b1), range(len_b2)):
                print l_b[b1][h1]
                print l_b[b2][h2]
                exit()
                l_b[b1][h1].connect_to[(b2, h2)] = scr
                l_b[b2][h2].connect_to_backward[(b1, h1)] = scr


        merged_triangle = 0
        for b1, b2 in itertools.combinations(range(len(l_b)), 2):
            len_b1 = len(l_b[b1])
            len_b2 = len(l_b[b2])
            for h1, h2 in itertools.product(range(len_b1), range(len_b2)):
                if l_b[b1][h1].flag == -1 or l_b[b2][h2].flag == -1:
                    continue
                scr = len(set(l_b[b1][h1].supporting_fragments) & set(l_b[b2][h2].supporting_fragments))
                if scr != 0:
                    l_b[b1][h1].connect_to[(b2, h2)] = scr
                    l_b[b2][h2].connect_to_backward[(b1, h1)] = scr


        for b in xrange(noBlock - 2):  # 2 for checking triangle
            l_triangle = []
            for h in xrange(ploidy):
                H = l_b[b][h]

                if H.flag == -1:
                    continue

                if len(H.connect_to) <= 1:
                    continue

                for hap1, hap2 in itertools.combinations(H.connect_to.keys(), 2):
                    b1, h1, b2, h2 = hap1 + hap2 if hap1[0] < hap2[0] else hap2 + hap1
                    H1, H2 = l_b[b1][h1], l_b[b2][h2]

                    if b1 == b2 or H1.flag == -1 or H2.flag == -1:
                        continue

                    if (b2, h2) in H1.connect_to:
                        scrs = (H.connect_to[(b1, h1)], H.connect_to[(b2, h2)], H1.connect_to[(b2, h2)])
                        l_triangle += [(b, h, b1, h1, b2, h2, min(scrs))]

            # print len(l_triangle)



            while len(l_triangle) > 0:
                l_triangle.sort( key=lambda x: x[-1])
                b, h, b1, h1, b2, h2, scr = l_triangle.pop()
                if scr > threshold:
                    H, H1, H2 = l_b[b][h], l_b[b1][h1], l_b[b2][h2]

                    if self.is_valid_to_merge_triangle(H,H1,H2):
                        # print b2, h2, b1, h1, b,h
                        if H1.flag != -1 and H2.flag != -1 :
                            H1H2_connections = set(H1.connect_to.keys()) & set(H2.connect_to.keys())
                            for b3, h3 in H1H2_connections:
                                H3 = l_b[b3][h3]
                                if H3.flag != -1 :
                                    scrs = (H1.connect_to[(b2, h2)], H1.connect_to[(b3, h3)], H2.connect_to[(b3, h3)])
                                    l_triangle += [(b, h, b1, h1, b3, h3, min(scrs))]
                                    l_triangle += [(b, h, b2, h2, b3, h3, min(scrs))]

                        if H1.flag != -1:
                            H += H1
                        else:
                            print 'A'
                        if H2.flag != -1:
                            H += H2
                        else:
                            print 'B'

                        merged_triangle += 1
                        if len(H.connect_to) > 1:

                            for hap1, hap2 in itertools.combinations(H.connect_to.keys(), 2):
                                b1, h1, b2, h2 = hap1 + hap2 if hap1[0] < hap2[0] else hap2 + hap1
                                H1, H2 = l_b[b1][h1], l_b[b2][h2]

                                if b1 == b2 or H1.flag == -1 or H2.flag == -1:
                                    continue



                                if (b2, h2) in H1.connect_to:
                                    scrs = (H.connect_to[(b1, h1)], H.connect_to[(b2, h2)], H1.connect_to[(b2, h2)])
                                    l_triangle += [(b, h, b1, h1, b2, h2, min(scrs))]

        return merged_triangle



    def is_valid_to_merge_2angle(self,H,H1,):
        bo_merge = True
        bo_merge = len(set(H.block_no) & set(H1.block_no)) == 0 and bo_merge
        bo_merge = (H.flag != -1 and H1.flag != -1 ) and bo_merge
        return bo_merge

    def find_2_angular_connections_by_merge_gready(self, l_b, threshold ,  window_size): #_______________too slow

        noBlock = len(l_b)

        l_candidate = []

        merged_2angle = 0
        for i in xrange(noBlock - window_size ):
            for b1, b2 in itertools.combinations(range(i,i+window_size), 2):
                ploidy_list1 = range(len(l_b[b1]))
                ploidy_list2 = range(len(l_b[b2]))
                for h1, h2 in itertools.product(ploidy_list1, ploidy_list2):

                    if l_b[b1][h1].flag == -1 or l_b[b2][h2].flag == -1:
                        continue
                    scr = len(set(l_b[b1][h1].supporting_fragments) & set(l_b[b2][h2].supporting_fragments))
                    if scr != 0:
                        l_b[b1][h1].connect_to[(b2, h2)] = scr
                        l_b[b2][h2].connect_to_backward[(b1, h1)] = scr
                        l_candidate += [(scr, b1,h1,b2,h2)]

            for scr, b1,h1,b2,h2 in sorted(l_candidate, key=lambda x:x[0], reverse=True):

                H = l_b[b1][h1]
                H1 = l_b[b2][h2]
                if scr > threshold and self.is_valid_to_merge_2angle(H,H1):
                    H += H1
                    merged_2angle +=1

        return merged_2angle


    def convert_haplotypes_into_blocks(self,l_haps, ploidy):



        l_haps.sort( key = lambda x:x.type[0])
        l_b = []
        for i,hap in enumerate(l_haps):
            my_block_no = i / ploidy
            my_inx = i % ploidy
            hap.block_no = [my_block_no]
            hap.name = str(my_block_no) + '_' + str(my_inx)
            hap.me = (my_block_no,my_inx)
            if my_inx == 0:
                l_b.append([hap])
            else:
                l_b[-1] += [hap]
        return l_b

    def update_skyline(self, ploidy,skyline,l_b):
        sky_tmp = np.zeros(skyline.__len__(), dtype=np.uint8)
        for b in l_b:
            for h in b:
                if h.flag == -1:
                    continue
                for j in h.haplotype:
                    sky_tmp[j] += 1

        for i in range(skyline.__len__()):
            if sky_tmp[i] != 0:
                skyline[i] = [(ploidy, sky_tmp[i])]

    def unused_frags(self, l_b,l_frag, l_frag_used):
        print np.sum(l_frag_used)
        print (len(l_frag_used) - np.sum(l_frag_used) ) / float(len(l_frag_used))

        return
        exit()
        for b in l_b:
            for h in b:
                h.supporting_fragment = []

        for i_frag, frag in enumerate(l_frag):
            if not l_frag_used[i_frag]:
                f = ast.literal_eval(frag[2])
                for b in l_b:
                    for h in b:
                        if h.flag == -1:
                            continue
                        bo_con = bo_over = False
                        cnt_overlap = 0
                        for pos in h.haplotype:
                            if pos in f:
                                bo_over = True

                                if f[pos] == h.haplotype[pos][0]:
                                    cnt_overlap += 1
                                else:
                                    bo_con = True
                        if not bo_con and bo_over and cnt_overlap != len(f):
                            h.supporting_fragments += [i_frag]

        for b in l_b:
            for h in b:
                if h.flag == -1:
                    continue
                print len(h.supporting_fragments)
        noBlock = len(l_b)

        merged_triangle = 0

        ploidy_list = range(ploidy)
        for b1, b2 in itertools.combinations(range(len(l_b)), 2):
            for h1, h2 in itertools.product(ploidy_list, ploidy_list):

                if l_b[b1][h1].flag == -1 or l_b[b2][h2].flag == -1:
                    continue
                scr = len(set(l_b[b1][h1].supporting_fragments) & set(l_b[b2][h2].supporting_fragments))
                if scr != 0:
                    l_b[b1][h1].connect_to[(b2, h2)] = scr
                    l_b[b2][h2].connect_to_backward[(b1, h1)] = scr


        exit()

    def sort_block_containing_less_ploidy(self,l_b, ploidy, no_blc_with_P):
        l_b_new = []
        cnt_blc = no_blc_with_P
        for b in l_b:
            this_block = []
            start_pos = sys.maxint
            if len(b) != ploidy:
                cur_block = cnt_blc
                cnt_blc += 1
                my_inx = 0
                for h in b:
                    if h.flag != -1:
                        start_pos = min (h.start_pos , start_pos)
                        h.block_no = [cur_block]
                        h.name = str(cur_block) + '_' + str(my_inx)
                        h.me = (cur_block,my_inx)
                        this_block += [h]
                        my_inx += 1
            else:
                for h in b:
                    if h.flag != -1:
                        start_pos = min (h.start_pos , start_pos)
                        this_block += [h]

            if len(this_block) != 0:
                l_b_new += [(start_pos,this_block)]



        l_b_new.sort(key = lambda x: x[0])
        l_return = []
        for i in l_b_new:
            l_return += [i[1]]
        return l_return


        # print merged_triangle

    def fit2ploidy(self,l_b,ploidy, scfLen):
        haps_overlap = []
        for i in range(scfLen):
            haps_overlap += [[]]
        for i_b, b in enumerate(l_b):
            for i_h, h in enumerate(b):
                if h.flag == -1:
                    continue
                for n in h.haplotype:
                    haps_overlap[n] += [(i_b, i_h)]

        for overlaped_blcs in haps_overlap:
            height = len(overlaped_blcs)
            if height > ploidy:
                mp_blc = {}
                for i in overlaped_blcs:
                    mp_blc[i[0]] = mp_blc.get(i[0], 0) + 1

                b_base = -1
                for k, v in mp_blc.items():
                    if v == ploidy:
                        b_base = k

                if b_base == -1:
                    continue

                l_base = []
                l_merge = []
                for i in overlaped_blcs:
                    if i[0] == b_base:
                        l_base += [i]
                    else:
                        l_merge += [i]

                for merge in l_merge:
                    mp_to_merge = {}

                    b2, h2 = merge
                    H1 = l_b[b2][h2]
                    mp2 = H1.haplotype

                    for base in l_base:
                        b1, h1 = base
                        H = l_b[b1][h1]
                        if b2 in H.block_no or H.flag == -1 or H1.flag == -1:
                            continue

                        mp1 = H.haplotype
                        m = mm = ml = mml = 0
                        for n in mp1:
                            if n in mp2:
                                if mp1[n][0] == mp2[n][0]:
                                    m += min(mp1[n][1], mp2[n][1])
                                    ml += 1
                                else:
                                    mm += min(mp1[n][1], mp2[n][1])
                                    mml += 1
                        mp_to_merge[mml] = mp_to_merge.get(mml, []) + [(base, ml, m, mml, mm)]
                    if  len((mp_to_merge)) == 0:
                         continue
                    sel = min(mp_to_merge.keys())
                    sel_b, sel_h = sorted(mp_to_merge[sel], key=lambda x: (-1 * x[4], x[1], x[2])[0])[0][0]
                    H = l_b[sel_b][sel_h]
                    H += H1
                    for pos, blc in mp2.items():
                        if merge not in haps_overlap[pos]:
                            print 'errrroor'
                        if (sel_b, sel_h) in haps_overlap[pos]:
                            haps_overlap[pos].remove(merge)
                        else:
                            haps_overlap[pos].remove(merge)
                            haps_overlap[pos] += [(sel_b, sel_h)]


    def phase_type_region(self, l_sorted_types, l_frag, l_frag_used, scf_hap_len, extention_allowed,
                          ploidy):  # phase_by_sorted_type

        scfLen = scf_hap_len
        first_frag_with_last_pos = [-1] * scfLen
        self.make_arr_start_from_frag(first_frag_with_last_pos, l_frag)
        last_frag_with_first_pos = [-1] * scfLen
        self.make_arr_end_to_frag(last_frag_with_first_pos, l_frag)

        l_haps_on_blocks = []

        skyline = [[] for i in range(scfLen)]
        cnt_type = 0
        collect_blocks_with_P = True
        down_to_two_phase = True
        down_to_two_connector = False
        no_blc_with_P = -1
        l_b = []
        for i in l_sorted_types:
            if i[0][2] == ploidy:
                if self.no_conflict_with_skyline(i, skyline, 'fully phased block'):
                    l_new_haps = self.phase_a_type(l_frag, l_frag_used, i, first_frag_with_last_pos,
                                                   last_frag_with_first_pos, skyline, cnt_type, extention_allowed, ploidy)
                    l_haps_on_blocks += l_new_haps
                continue

            elif collect_blocks_with_P:
                collect_blocks_with_P = False
                # l_blocks_sorted = self.sort_blocks(l_haps_on_blocks)
                l_b = self.convert_haplotypes_into_blocks(l_haps_on_blocks, ploidy)


                # self.bridge_blocks_max_weighted_matching(l_b, l_frag, l_frag_used,
                #                                                           first_frag_with_last_pos,
                #                                                           last_frag_with_first_pos,
                #                                                           ploidy)

                # l_obtained_haps = self.merge_consecutive_blocks_by_overlaps(l_blocks_sorted, ploidy,
                #
                #                                                          l_hap2hap_adj_blocks)
                import datetime
                t1 = datetime.datetime.now()
                self.find_3_4angular_connections_by_merge_gready(l_b, 7 ,ploidy) #parameter
                #print 'One triangle time:', datetime.datetime.now() - t1

                self.find_3_4angular_connections_by_merge_gready(l_b, 2 ,ploidy) #parameter
                self.find_2_angular_connections_by_merge_gready(l_b, 15 , 7) #parameter
                self.find_2_angular_connections_by_merge_gready(l_b, 5 , 7) #parameter
                self.find_2_angular_connections_by_merge_gready(l_b, 2 , 7) #parameter
                self.find_2_angular_connections_by_merge_gready(l_b, 1 , 7) #parameter


                self.update_skyline(ploidy,skyline, l_b)
                no_blc_with_P = len(l_b)
                self.fit2ploidy(l_b, ploidy, scfLen)

            if down_to_two_phase:
                cur_ploidy = i[0][2]
                if self.no_conflict_with_skyline(i, skyline, 'tail block'):

                    extention_allowed = False
                    l_new_haps = self.phase_a_type(l_frag, l_frag_used, i, first_frag_with_last_pos,
                                                   last_frag_with_first_pos, skyline, cnt_type, extention_allowed,
                                                   cur_ploidy)
                    if len(l_new_haps) >= 2:
                        cnt_blocks = len(l_b)
                        l_b += [l_new_haps]
                        for no_bl, bl in enumerate(l_b[-1]):
                            bl.name = str(cnt_blocks) + '_' + str(no_bl)
            else:
                break
        if down_to_two_connector :
            l_b = self.sort_block_containing_less_ploidy(l_b,ploidy, no_blc_with_P)
            self.find_2_angular_connections_by_merge_gready(l_b, 5, 7)  # parameter
            self.find_2_angular_connections_by_merge_gready(l_b, 2, 7)  # parameter


        #self.fit2ploidy(l_b, ploidy,scfLen)





        return l_b#, l_hap2hap, l_blocks_sorted

    def sort_blocks(self, l_ans):  # phase_by_sorted_type
        l = []
        l_block_to_inx = []
        prev_block = -1
        cnt = 0
        while cnt < len(l_ans):
            bc_no = l_ans[cnt].block_no[0]
            l.append((l_ans[cnt].type[0], bc_no, l_ans[cnt].flag))
            l_block_to_inx.append(cnt)
            cnt += l_ans[cnt].flag
        l = sorted(l, key=lambda x: x[0])
        l_blocks = []
        cnt = 0
        for i in l:
            l_blocks.append([])
            for j in range(i[2]):
                l_blocks[-1].append(l_ans[l_block_to_inx[i[1]] + j])
        return l_blocks

    def overlap(self, h1, h2):  # bridge_blocks_alg2
        bo = False
        no_sim = 0
        no_dis = 0
        qual_sim = 0
        qual_dis = 0

        for pos in h1.haplotype.viewkeys() & h2.haplotype.viewkeys():
            bo = True
            allele1 = h1.haplotype[pos]
            allele2 = h2.haplotype[pos]
            if allele1[0] == allele2[0]:
                no_sim += 1
                qual_sim += (allele1[1] + allele2[1])
            else:
                no_dis += 1
                qual_dis += (allele1[1] + allele2[1])

        # for i in range(max(h1.start_pos, h2.start_pos), min(h1.end_pos, h2.end_pos)):
        #     nc1 = h1.haplotype[i - h1.start_pos]
        #     nc1_qual = h1.quality[i - h1.start_pos]
        #
        #     nc2 = h2.haplotype[i - h2.start_pos]
        #     nc2_qual = h2.quality[i - h2.start_pos]
        #
        #     if nc1 == '-' or nc2 == '-':
        #         continue
        #     elif nc1 == nc2:
        #         no_sim += 1
        #         qual_sim += (nc1_qual + nc2_qual)
        #         bo = True
        #     else:
        #         no_dis += 1
        #         qual_dis += (nc1_qual + nc2_qual)
        #         bo = True


        return (bo, no_sim, no_dis, qual_sim, qual_dis)


    def overlap_score(self, h1, h2):  # bridge_blocks_alg2
        bo = False
        no_sim = 0
        no_dis = 0
        qual_sim = 0
        qual_dis = 0

        for pos in h1.haplotype.viewkeys() & h2.haplotype.viewkeys():
            bo = True
            allele1 = h1.haplotype[pos]
            allele2 = h2.haplotype[pos]
            if allele1[0] == allele2[0]:
                no_sim += 1
                qual_sim += (allele1[1] + allele2[1])
            else:
                no_dis += 1
                qual_dis += (allele1[1] + allele2[1])

        # for i in range(max(h1.start_pos, h2.start_pos), min(h1.end_pos, h2.end_pos)):
        #     nc1 = h1.haplotype[i - h1.start_pos]
        #     nc1_qual = h1.quality[i - h1.start_pos]
        #
        #     nc2 = h2.haplotype[i - h2.start_pos]
        #     nc2_qual = h2.quality[i - h2.start_pos]
        #
        #     if nc1 == '-' or nc2 == '-':
        #         continue
        #     elif nc1 == nc2:
        #         no_sim += 1
        #         qual_sim += (nc1_qual + nc2_qual)
        #         bo = True
        #     else:
        #         no_dis += 1
        #         qual_dis += (nc1_qual + nc2_qual)
        #         bo = True
        return 0 if not bo else -1 if no_dis > 0 else qual_sim




    def bridge_blocks(self, l_b, l_frags, l_frag_used, scf_info, start_from_frag, end_till_frag):

        for i in range(0, len(l_b) - 1):
            any_change = True
            cnt_tmp = 0
            l_2nd_block = []

            while any_change:
                l_overlap = []
                l_sim = []
                l_dis = []
                l_sim_qual = []
                l_dis_qual = []
                any_change = False
                cnt_j = 0

                for j in l_b[i]:
                    cnt_j += 1
                    l_overlap.append([])
                    l_sim.append([])
                    l_dis.append([])
                    l_sim_qual.append([])
                    l_dis_qual.append([])
                    cnt_k = 0

                    for k in l_b[i + 1]:
                        cnt_k += 1
                        mn = min(j.start_pos, k.start_pos)
                        tmp = overlap(j, k)
                        l_sim[-1].append(tmp[1])
                        l_dis[-1].append(tmp[2])
                        l_sim_qual[-1].append(tmp[3])
                        l_dis_qual[-1].append(tmp[4])
                        if tmp[0]:
                            l_overlap[-1].append(1)
                        else:
                            l_overlap[-1].append(0)

                no_row = len(l_b[i])
                no_col = len(l_b[i + 1])

                if no_row == 0 or no_col == 0:
                    continue

                l_cand = []
                for j in range(no_col):
                    no_sim = 0
                    no_dis = 0
                    no_sim_nc = 0
                    no_dis_nc = 0
                    no_sim_qual = 0
                    no_dis_qual = 0
                    bo = True
                    for k in range(no_row):
                        if l_overlap[k][j] == 0:
                            bo = False
                        else:
                            no_sim += 0 if l_sim[k][j] == 0 else 1
                            no_dis += 0 if l_dis[k][j] == 0 else 1
                            no_sim_nc += l_sim[k][j]
                            no_dis_nc += l_dis[k][j]
                            no_sim_qual += l_sim_qual[k][j]
                            no_dis_qual += l_dis_qual[k][j]

                    if bo:
                        l_cand.append(('col', j, no_sim, -1*abs(abs(no_row - no_dis) - 1), no_sim_nc, -1 * no_dis_nc,
                                       no_sim_qual, -1 * no_dis_qual))

                for k in range(no_row):
                    no_sim = 0
                    no_dis = 0
                    no_sim_nc = 0
                    no_dis_nc = 0
                    no_sim_qual = 0
                    no_dis_qual = 0
                    bo = True
                    for j in range(no_col):
                        if l_overlap[k][j] == 0:
                            bo = False
                        else:
                            no_sim += 0 if l_sim[k][j] == 0 else 1
                            no_dis += 0 if l_dis[k][j] == 0 else 1
                            no_sim_nc += l_sim[k][j]
                            no_dis_nc += l_dis[k][j]
                            no_sim_qual += l_sim_qual[k][j]
                            no_dis_qual += l_dis_qual[k][j]

                    if bo:
                        l_cand.append(('row', k, no_sim, -1*abs(abs(no_row - no_dis) - 1), no_sim_nc, -1 * no_dis_nc,
                                       no_sim_qual, -1 * no_dis_qual))

                if len(l_cand) != 0:
                    sel = sorted(l_cand, key=lambda x: (x[3], x[4], x[6] + x[7]), reverse=True)[0]
                else:
                    continue

                if sel[0] == 'row':
                    any_change = True
                    j = sel[1]
                    l_sel_tmp = []
                    for k in range(no_col):
                        if l_dis[j][k] == 0:
                            l_sel_tmp.append(k)
                    if len(l_sel_tmp) == 1:
                        sel_col = l_sel_tmp[0]
                    elif len(l_sel_tmp) > 1:
                        tmp_mx = -1
                        for k in l_sel_tmp:
                            if l_sim_qual[j][k] > tmp_mx:
                                tmp_mx = l_sim_qual[j][k]
                                sel_col = k
                    elif len(l_sel_tmp) == 0:
                        tmp_score_mx = -1 * sys.maxint
                        sel_col = -1
                        for k in range(no_col):
                            scr = l_sim_qual[j][k] - l_dis_qual[j][k]
                            if scr > tmp_score_mx:
                                tmp_score_mx = scr
                                sel_col = k
                    mergedHaps = l_b[i + 1][sel_col] + l_b[i][j]
                    del l_b[i][j]
                    del l_b[i + 1][sel_col]
                    l_2nd_block.append(mergedHaps)

                if sel[0] == 'col':
                    any_change = True
                    k = sel[1]
                    l_sel_tmp = []
                    for j in range(no_row):
                        if l_dis[j][k] == 0:
                            l_sel_tmp.append(j)
                    if len(l_sel_tmp) == 1:
                        sel_row = l_sel_tmp[0]
                    elif len(l_sel_tmp) > 1:
                        tmp_mx = -1
                        for j in l_sel_tmp:
                            if l_sim_qual[j][k] > tmp_mx:
                                tmp_mx = l_sim_qual[j][k]
                                sel_row = j
                    elif len(l_sel_tmp) == 0:
                        tmp_score_mx = -1 * sys.maxint
                        sel_row = -1
                        for j in range(no_row):
                            scr = l_sim_qual[j][k] - l_dis_qual[j][k]
                            if scr > tmp_score_mx:
                                tmp_score_mx = scr
                                sel_row = j

                    mergedHaps = l_b[i][sel_row] + l_b[i + 1][k]
                    del l_b[i][sel_row]
                    del l_b[i + 1][k]
                    l_2nd_block.append(mergedHaps)

            if len(l_b[i]) == len(l_b[i + 1]) == 1:
                mergedHaps = l_b[i][0] + l_b[i + 1][0]
                del l_b[i][0]
                del l_b[i + 1][0]
                l_2nd_block.append(mergedHaps)

            for j in l_2nd_block:
                l_b[i + 1].append(j)




    def bridge_blocks_max_weighted_matching(self, l_b, l_frags, l_frag_used, start_from_frag, end_till_frag,
                           ploidy):  # phase_by_sorted_type


        l_hap2hap = [] ; l_hap2hap_by_blocks = []
        # bo_overlap, no_sim, no_dis, qual_sim, qual_dis
        l_overlap_two_block = np.zeros((ploidy,ploidy),dtype=np.int32)




        for i_block in range(0, len(l_b) - 1):
            l_candidate = []
            for cnt_j, j in enumerate(l_b[i_block]):
                if j.flag == -1:
                    continue
                l_candidate += [[]]
                for cnt_k, k in enumerate(l_b[i_block + 1]):
                    if k.flag == -1:
                        continue
                    l_overlap_two_block[cnt_j][cnt_k] += [self.overlap_score(j, k)]
                    if l_overlap_two_block[cnt_j][cnt_k] >= 0 :
                        l_candidate[-1] += [((cnt_j, cnt_k), l_overlap_two_block[cnt_j][cnt_k])]

            bo_stay_in_loop = True
            l_connect = []
            while bo_stay_in_loop:
                bo_stay_in_loop = False
                mx_val_cand = -1
                selected_cand = -1


                for i_cand, cands  in enumerate(l_candidate):
                    if cands.__len__() == 1 and cands[0][1] > mx_val_cand:
                        selected_cand, mx_val_cand = cands[0] , cands[0][1]


                if selected_cand != -1:
                    bo_stay_in_loop = True
                    l_connect += [selected_cand[0]]


                    ### remove selected cand from both sets
                    l_new_candidate = []
                    for i_cand, cands in enumerate(l_candidate):
                        l_new_candidate += [[]]
                        rm_i , rm_j = selected_cand[0]
                        for cand in cands:
                            if cand[0][0] == rm_i or cand[0][1] == rm_j:
                                continue
                            else:
                                l_new_candidate[-1] += [cand]
                    l_candidate = l_new_candidate
            for i, j in l_connect:
                l_b [i_block + 1][j] += l_b [i_block][i]







    def bridge_blocks_alg2(self, l_b, l_frags, l_frag_used, start_from_frag, end_till_frag,
                           ploidy):  # phase_by_sorted_type

        l_hap2hap = [] ; l_hap2hap_by_blocks = []

        for i in range(0, len(l_b) - 1):
            l_hap2hap_just_cur_block = []
            any_change = True
            cnt_tmp = 0
            l_2nd_block = []



            used_hap = np.zeros((2, ploidy), dtype=bool)


            while any_change:
                any_change = False
                #l_overlap = []; l_sim = [] ; l_dis = []; l_sim_qual = []; l_dis_qual = []

                l_overlap_two_block = []
                for cnt_j, j in enumerate(l_b[i]):
                    if used_hap[0][cnt_j]:
                        k = len(l_b[i + 1])
                        l_overlap.append(k * [2])
                        l_sim.append(k * [-1])
                        l_dis.append(k * [-1])
                        l_sim_qual.append(k * [-1])
                        l_dis_qual.append(k * [-1])
                        continue

                    # l_overlap.append([]); l_sim.append([]); l_dis.append([]); l_sim_qual.append([]); l_dis_qual.append([])
                    l_overlap_two_block.append([])

                    for cnt_k, k in enumerate(l_b[i + 1]):

                        if used_hap[1][cnt_k]:
                            l_sim[-1].append(-1)
                            l_dis[-1].append(-1)
                            l_sim_qual[-1].append(-1)
                            l_dis_qual[-1].append(-1)
                            l_overlap[-1].append(2)
                            continue
                        #bo_overlap, no_sim, no_dis, qual_sim, qual_dis
                        l_overlap_two_block[-1] += [self.overlap(j, k)]
                        # l_sim[-1].append(tmp[1])
                        # l_dis[-1].append(tmp[2])
                        # l_sim_qual[-1].append(tmp[3])
                        # l_dis_qual[-1].append(tmp[4])
                        # l_overlap[-1].append(1)

                no_row = len(l_b[i])
                no_col = len(l_b[i + 1])

                if used_hap[0].count(False) == 0 or used_hap[1].count(False) == 0:
                    continue

                l_cand = []
                for j in range(no_col):

                    if used_hap[1][j]:
                        continue

                    no_sim = no_dis = no_sim_nc = no_dis_nc = no_sim_qual = no_dis_qual = 0
                    bo = False
                    for k in range(no_row):
                        if used_hap[0][k]:
                            continue
                        bo = True

                        if l_overlap[k][j] != 1:
                            bo = False
                            break
                        else:
                            no_sim += 0 if l_sim[k][j] == 0 else 1
                            no_dis += 0 if l_dis[k][j] == 0 else 1
                            no_sim_nc += l_sim[k][j]
                            no_dis_nc += l_dis[k][j]
                            no_sim_qual += l_sim_qual[k][j]
                            no_dis_qual += l_dis_qual[k][j]

                    if bo:
                        l_cand.append(('col', j, no_sim, -1*abs(abs(no_row - no_dis) - 1), no_sim_nc, -1 * no_dis_nc,
                                       no_sim_qual, -1 * no_dis_qual))

                for k in range(no_row):
                    if used_hap[0][k]:
                        continue

                    no_sim = 0
                    no_dis = 0
                    no_sim_nc = 0
                    no_dis_nc = 0
                    no_sim_qual = 0
                    no_dis_qual = 0
                    bo = False
                    for j in range(no_col):
                        if used_hap[1][j]:
                            continue

                        bo = True

                        if l_overlap[k][j] != 1:
                            bo = False
                            break
                        else:
                            no_sim += 0 if l_sim[k][j] == 0 else 1
                            no_dis += 0 if l_dis[k][j] == 0 else 1
                            no_sim_nc += l_sim[k][j]
                            no_dis_nc += l_dis[k][j]
                            no_sim_qual += l_sim_qual[k][j]
                            no_dis_qual += l_dis_qual[k][j]

                    if bo:
                        l_cand.append(('row', k, no_sim, -1*abs(abs(no_row - no_dis) - 1), no_sim_nc, -1 * no_dis_nc,
                                       no_sim_qual, -1 * no_dis_qual))

                if len(l_cand) != 0:
                    sel = sorted(l_cand, key=lambda x: (x[3], x[4], x[6] + x[7]), reverse=True)[0]
                else:
                    continue

                if sel[0] == 'row':
                    any_change = True
                    j = sel[1]
                    l_sel_tmp = []
                    for k in range(no_col):
                        if not used_hap[1][k] and l_dis[j][k] == 0:
                            l_sel_tmp.append(k)
                    if len(l_sel_tmp) == 1:
                        sel_col = l_sel_tmp[0]
                    elif len(l_sel_tmp) > 1:
                        tmp_mx = -1
                        for k in l_sel_tmp:
                            if l_sim_qual[j][k] > tmp_mx:
                                tmp_mx = l_sim_qual[j][k]
                                sel_col = k
                    elif len(l_sel_tmp) == 0:
                        tmp_score_mx = -1 * sys.maxint
                        sel_col = -1
                        for k in range(no_col):
                            if used_hap[1][k]:
                                continue
                            scr = l_sim_qual[j][k] - l_dis_qual[j][k]
                            if scr > tmp_score_mx:
                                tmp_score_mx = scr
                                sel_col = k
                    mergedHaps = l_b[i + 1][sel_col] + l_b[i][j]
                    l_hap2hap.append(((i, j), (i + 1, sel_col)))
                    l_hap2hap_just_cur_block.append(((i, j), (i + 1, sel_col)))
                    used_hap[0][j] = True
                    used_hap[1][sel_col] = True

                    l_2nd_block.append(mergedHaps)

                if sel[0] == 'col':
                    any_change = True
                    k = sel[1]
                    l_sel_tmp = []
                    for j in range(no_row):
                        if not used_hap[0][j] and l_dis[j][k] == 0:
                            l_sel_tmp.append(j)

                    if len(l_sel_tmp) == 1:
                        sel_row = l_sel_tmp[0]
                    elif len(l_sel_tmp) > 1:
                        tmp_mx = -1
                        for j in l_sel_tmp:
                            if l_sim_qual[j][k] > tmp_mx:
                                tmp_mx = l_sim_qual[j][k]
                                sel_row = j
                    elif len(l_sel_tmp) == 0:
                        tmp_score_mx = -1 * sys.maxint
                        sel_row = -1
                        for j in range(no_row):
                            if used_hap[0][j]:
                                continue
                            scr = l_sim_qual[j][k] - l_dis_qual[j][k]
                            if scr > tmp_score_mx:
                                tmp_score_mx = scr
                                sel_row = j

                    mergedHaps = l_b[i][sel_row] + l_b[i + 1][k]
                    l_hap2hap.append(((i, sel_row), (i + 1, k)))
                    l_hap2hap_just_cur_block.append(((i, sel_row), (i + 1, k)))

                    used_hap[0][sel_row] = True
                    used_hap[1][k] = True

                    l_2nd_block.append(mergedHaps)

            if used_hap[0].count(False) == used_hap[1].count(False) == 1:
                inx1 = used_hap[0].index(False)
                inx2 = used_hap[1].index(False)
                mergedHaps = l_b[i][0] + l_b[i + 1][0]
                l_hap2hap.append(((i, inx1), (i + 1, inx2)))
                l_hap2hap_just_cur_block.append(((i, inx1), (i + 1, inx2)))
                used_hap[0][inx1] = True
                used_hap[1][inx2] = True
                l_2nd_block.append(mergedHaps)

            l_hap2hap_by_blocks.append(l_hap2hap_just_cur_block)
        return l_hap2hap, l_hap2hap_by_blocks

    def score_model(self, w1, w2):  # infer_edges_from_blocks
        model_sum = True
        mdeol_min = False
        if model_sum:
            return w1 + w2

    def overlap_conflict_weighted_with_frag(self, subHap, fragment):  # infer_edges_from_blocks
        start_pos = int(fragment[1])
        frag = fragment[2]
        copies = int(fragment[0])
        frag_len = len(frag)
        t = ast.literal_eval(frag)#not reachable for change ast lib
        w = 0
        bo_overlap = False
        bo_conflict = False
        for i in range(frag_len):
            if i + start_pos in range(subHap.start_pos, subHap.end_pos):
                if subHap.haplotype[i + start_pos - subHap.start_pos] != '-' and frag[i] != '-':
                    bo_overlap = True
                    w += copies  # simplest model
                if subHap.haplotype[i + start_pos - subHap.start_pos] != frag[i]\
                        and subHap.haplotype[i + start_pos - subHap.start_pos] != '-'\
                        and frag[i] != '-':
                    bo_conflict = True
        return bo_overlap, bo_conflict, w

    def overlap_conflict_weighted(self, subHap, fragment):  # infer_edges_from_blocks
        copies, start_pos, frag_tmp, end_pos = fragment
        t = {}

        frag = ast.literal_eval(frag_tmp)#not reachable for changing the ast library

        bo_overlap = False
        bo_conflict = False
        w = 0
        for pos in frag:
            if pos in subHap.haplotype:
                allele, copies = subHap.haplotype[pos]
                bo_overlap = True
                bo_conflict = True if allele != frag[pos] else False
                w += copies  ###---### in vazn majmooe vazn oonhaist ke overlap daran o overlap nadaran



        # for i in seg:
        #     t[i] = str(seg[i])
        #
        # w = 0
        # bo_overlap = False
        # bo_conflict = False
        # for i in sorted(t.keys()):
        #     if subHap.start_pos <= i < subHap.end_pos:
        #         if subHap.haplotype[i - subHap.start_pos] != '-':  # and frag[i] != '-':
        #             bo_overlap = True
        #             w += copies  # simplest model
        #         if subHap.haplotype[i - subHap.start_pos] != t[i]\
        #                 and subHap.haplotype[i - subHap.start_pos] != '-':  # and frag[i] != '-':
        #             bo_conflict = True
        return bo_overlap, bo_conflict, w

    def seg_to_frag(self, seg):
        t = ast.literal_eval(seg)#not reachable for changing ast library
        mn = min(t.keys())
        mx = max(t.keys()) + 1
        ln = mx - mn
        fr = ln * ['-']

        for j in t:
            fr[j - mn] = str(t[j])
        return ''.join(fr)

    def infer_edges_from_blocks(self, l_b, l_frag, l_frag_used, mp_g):  # make_block_graph_alg2

        cnt = 0
        for f_inx in xrange(len(l_frag)): ###---### mitoonam az index estefade konam faghat read hai ke overlap daran biaram
            copies = l_frag[f_inx][0]
            start_pos = l_frag[f_inx][1]
            # if self.uniq_frag_mode == 'seg':
            #	frag = self.seg_to_frag(l_frag[f_inx][2])
            # else:
            #	frag			=	l_frag[f_inx][2]
            #	#print frag
            # frag_len		=	len(frag)
            # end_pos		= 	start_pos + frag_len

            l_match_to_blocks = []
            l_match_to_blocks_weight = []

            for b in range(0, len(l_b)):
                cnt_shared = 0
                whichHap = -1
                bo_overlapAll = True
                for i in range(0, len(l_b[b])):

                    bo_overlap, bo_conflict, w = self.overlap_conflict_weighted(l_b[b][i], l_frag[f_inx])
                    if not bo_overlap:
                        bo_overlapAll = False
                        break
                    if not bo_conflict:
                        cnt_shared += 1
                        whichHap = i
                    if cnt_shared > 1:
                        break

                if cnt_shared == 1 and bo_overlapAll:
                    l_match_to_blocks.append((b, whichHap))
                    l_match_to_blocks_weight.append(w)

            no_nodes = len(l_match_to_blocks)
            if no_nodes > 1:
                for k in itertools.combinations(range(no_nodes), 2):
                    n1 = l_match_to_blocks[k[0]]
                    n2 = l_match_to_blocks[k[1]]
                    i_n1 = k[0]
                    i_n2 = k[1]

                    w1 = l_match_to_blocks_weight[k[0]]
                    w2 = l_match_to_blocks_weight[k[1]]

                    if n1 in mp_g:
                        if n2 in mp_g[n1]:
                            mp_g[n1][n2] += self.score_model(w1, w2)
                        else:
                            mp_g[n1][n2] = self.score_model(w1, w2)
                    else:
                        mp_g[n1] = {}
                        mp_g[n1][n2] = self.score_model(w1, w2)

                    n1, n2 = n2, n1

                    if n1 in mp_g:
                        if n2 in mp_g[n1]:
                            mp_g[n1][n2] += self.score_model(w1, w2)
                        else:
                            mp_g[n1][n2] = self.score_model(w1, w2)
                    else:
                        mp_g[n1] = {}
                        mp_g[n1][n2] = self.score_model(w1, w2)

    def convert_mp_mp_to_graph(self, mp_g, g, no_blocks, l_adjacent, ploidy):  # make_block_graph_alg2
        for i in range(no_blocks):
            for j in range(ploidy):
                g.add_node((i, j))
        for i in l_adjacent:
            g.add_edge(i[0], i[1], weight=1000)

        for i in mp_g:
            for j in mp_g[i]:
                g.add_edge(i, j, weight=mp_g[i][j])

    def convert_l_to_graph(self, g, no_blocks, l_adjacent, ploidy):  # make_block_graph_alg2
        for i in range(no_blocks):
            for j in range(ploidy):
                g.add_node((i, j))
        for i in l_adjacent:
            g.add_edge(i[0], i[1], weight=1000)

    def produce_conflict(self, g, n1, n2, no_blocks):

        faster_method = True
        if faster_method:
            sp = nx.descendants(g, n1)
            sw = [0] * no_blocks
            sw[n1[0]] += 1

            for i in sp:
                sw[i[0]] += 1

            for i in sw:
                if i > 1:
                    return True
            return False
        else:
            sp, _ = nx.single_source_dijkstra(g, n1)
            mp = {}
            for i in sp:
                if i[0] in mp:
                    if mp[i[0]] != i[1]:
                        return True
                else:
                    mp[i[0]] = i[1]
            return False

    def kruskal(self, g, gk, no_blocks):
        for i in g:
            gk.add_node(i)

        for a, b, dct in sorted(g.edges(data=True), key=lambda (a, b, dct): dct['weight'], reverse=True):
            gk.add_edge(a, b, weight=1)
            if self.produce_conflict(gk, a, b, no_blocks):
                gk.remove_edge(a, b)
        return

    def ordered_edge(self, u, v):
        if u[0] > v[0]:
            return v, u
        else:
            return u, v

    def nodes_of_a_color(self, nodes):
        l = []
        for i in nodes:
            l += [i[0]]
        return set(l)

    def print_matchig_parts(self, mp, no_blocks):
        print mp
        for i in mp:
            print str(i) + ': ',
            l = []
            for j in mp[i]:
                l += [j[0]]
            for j in sorted(l):
                print j,
            print

        return

    def graph_coloring(self, g, g1, ploidy, no_b):
        mx_weight = -1

        for (u, v, d) in g.edges(data=True):
            if u[0] == v[0]:
                print 'error12'
                exit()

            mx_weight = d['weight'] if d['weight'] > mx_weight else mx_weight

            for j in range(ploidy):
                if j != v[1]:
                    v1 = (v[0], j)
                    u1, v1 = self.ordered_edge(u, v1)
                    if (u1, v1) in g1.edges() or (v1, u1) in g1.edges():
                        g1[u1][v1]['weight'] = g1[u1][v1]['weight'] + d['weight']
                    else:
                        g1.add_edge(u1, v1, weight=d['weight'])

            for j in range(ploidy):
                if j != u[1]:
                    u1 = (u[0], j)
                    u1, v1 = self.ordered_edge(u1, v)

                    if (u1, v1) in g1.edges() or (v1, u1) in g1.edges():
                        g1[u1][v1]['weight'] = g1[u1][v1]['weight'] + d['weight']
                    else:
                        g1.add_edge(u1, v1, weight=d['weight'])

        for i in range(no_b):
            for j in range(ploidy):
                for k in range(j + 1, ploidy):
                    g1.add_edge((i, j), (i, k), weight=3 * mx_weight)

        print 'strategy_largest_first'

        d = nx.coloring.greedy_color(g1, strategy=nx.coloring.strategy_largest_first)
        colors = max(d.values()) + 1

        mp = {}
        for i in d:
            mp[d[i]] = mp.get(d[i], []) + [i]

        if len(mp) == ploidy:
            l_cc = []
            for i in mp:
                l_cc.append(mp[i])
            return l_cc

        ch = True

        while ch:
            ch = False
            cnt = mp.keys()[0]
            for cnt in mp.keys():
                s = self.nodes_of_a_color(mp[cnt])

                which_col = []
                no_col = 0
                for i in mp.keys():
                    if i != cnt:
                        s2 = self.nodes_of_a_color(mp[i])
                        if len(s & s2) == 0:
                            no_col += 1
                            which_col += [i]
                if no_col == 1:
                    if which_col[0] > cnt:
                        mp[cnt] += mp[which_col[0]]
                        del mp[which_col[0]]
                    else:
                        mp[which_col[0]] += mp[cnt]
                        del mp[cnt]

                    ch = True
                    break

        l_cc = []
        for i in mp:
            l_cc.append(mp[i])
        return l_cc

        exit()

        mp_best_seed = {}
        for i in mp.keys():
            for j in mp[i]:
                exit()

        for i in range(len(colors)):
            for j in range(i + 1, len(colors)):
                pass

        for i in range(max(d.values()) + 1):
            for j in d:
                if d[j] == i:
                    print j,
            print

        print 'greedy_color'
        print 'strategy_random_sequential'
        d = nx.coloring.greedy_color(g1, strategy=nx.coloring.strategy_random_sequential)
        print max(d.values())
        print 'strategy_smallest_last'
        d = nx.coloring.greedy_color(g1, strategy=nx.coloring.strategy_smallest_last)
        print max(d.values())
        print 'strategy_independent_set'
        d = nx.coloring.greedy_color(g1, strategy=nx.coloring.strategy_independent_set)
        print max(d.values())
        print 'strategy_connected_sequential'
        d = nx.coloring.greedy_color(g1, strategy=nx.coloring.strategy_connected_sequential)
        print max(d.values())
        print 'strategy_connected_sequential_dfs'
        d = nx.coloring.greedy_color(g1, strategy=nx.coloring.strategy_connected_sequential_dfs)
        print max(d.values())
        print 'strategy_connected_sequential_bfs'
        d = nx.coloring.greedy_color(g1, strategy=nx.coloring.strategy_connected_sequential_bfs)
        print max(d.values())
        print 'strategy_saturation_largest_first'
        d = nx.coloring.greedy_color(g1, strategy=nx.coloring.strategy_saturation_largest_first)
        print max(d.values())
        exit()

    def merge_consecutive_blocks_by_overlaps(self, l_b, ploidy, l_hap2hap_adj_blocks):  # phase_by_sorted_type
#        from visualClass import vis

        l_consecutive_blocks = []  # when the adjacent blocks have 6 matches
        for i in l_hap2hap_adj_blocks:
            if len(i) == ploidy:
                for j in i:
                    l_consecutive_blocks.append(j)

        no_blocks = len(l_b)

        g = nx.Graph()
        self.convert_l_to_graph(g, no_blocks, l_consecutive_blocks, ploidy)

        if True:
            from visualClass import vis
            v = vis()
            v.plot_multipartide_g(g)

        l_cc = []
        l_cc = nx.connected_components(g)

        l_ans = []
        for i in l_cc:
            l_sorted_haps = sorted(i, key=lambda x: x[0])
            n = l_sorted_haps[0]
            merged_hap = l_b[n[0]][n[1]].copy()

            for jj in range(1, len(l_sorted_haps)):
                n = l_sorted_haps[jj]
                merged_hap = merged_hap + l_b[n[0]][n[1]]

            l_ans.append(merged_hap)

        l_for_sort1 = []
        for i in l_ans:
            l_for_sort1.append((i.name, i))

        l_sorted = []
        for i in sorted(l_for_sort1, key=lambda x: (x[0])):
            l_sorted.append(i[1])

        return l_sorted

    def make_block_graph_alg2(self, l_b, l_frag, l_frag_used, l_consecutive_blocks, ploidy):

        mp_g = {}
        self.infer_edges_from_blocks(l_b, l_frag, l_frag_used, mp_g)
        del l_frag
        del l_frag_used

        no_blocks = len(l_b)

        g = nx.Graph()
        self.convert_mp_mp_to_graph(mp_g, g, no_blocks, l_consecutive_blocks, ploidy)

        mp_merged = {}

        l_cc = []
        if self.algName == 'kruskal':
            g_kruscal = nx.Graph()
            self.kruskal(g, g_kruscal, no_blocks)

            for i in nx.connected_components(g_kruscal):
                l_cc.append([])
                for j in i:
                    if j in mp_merged:
                        l_cc[-1] += mp_merged[j]
                    else:
                        l_cc[-1].append(j)
            return l_cc

        elif self.algName == 'coloring':
            g_coloring = nx.Graph()
            l_cc = self.graph_coloring(g, g_coloring, ploidy, no_blocks)
            return l_cc

    def make_haps_on_snp(self, scf_info, l_h):  # filling_gaps
        ln = scf_info.hap_len
        l = []
        for i in range(ln):
            l.append([])

        for i_h, h in enumerate(l_h):
            for i in range(h.start_pos, h.end_pos):
                if h.haplotype[i - h.start_pos] != '-':
                    if i_h not in l[i]:
                        l[i].append(i_h)
        return l

    def overlap_conflict_Region(self, subHap, fragment):  # filling_gaps
        start_pos = int(fragment[1])
        frag = fragment[2]
        copies = int(fragment[0])
        frag_len = len(frag)

        bo_overlap = False
        bo_conflict = False
        for i in range(frag_len):
            if i + start_pos in range(subHap.start_pos, subHap.end_pos):
                if subHap.haplotype[i + start_pos - subHap.start_pos] != '-' and frag[i] != '-':
                    bo_overlap = True
                if subHap.haplotype[i + start_pos - subHap.start_pos] != frag[i]\
                        and subHap.haplotype[i + start_pos - subHap.start_pos] != '-'\
                        and frag[i] != '-':
                    bo_conflict = True

        return bo_overlap, bo_conflict

    def overlap_conflict_Region_seg(self, subHap, fragment):  # read format segment
        # _, start_pos, tmp_frag, end_pos = fragment
        # frag = ast.literal_eval(tmp_frag)#here

        _, start_pos, frag, end_pos,_ = fragment
        bo_overlap = False
        bo_conflict = False
        for pos in frag:
            if pos in subHap.haplotype:
                bo_overlap = True
                if subHap.haplotype[pos][0] != frag[pos]:
                    bo_conflict = True
        return bo_overlap, bo_conflict


    def filling_gaps(self, scf_info, l_blocks_sorted, l_frag, l_frag_used, l_h):  # phase_by_sorted_type

        l_hap_on_snps = self.make_haps_on_snp(scf_info, l_h)

        l_hap_absorb_frag = []
        for i in range(len(l_h)):
            l_hap_absorb_frag.append([])

        for f_inx, f in enumerate(l_frag):
            if l_frag_used[f_inx]:
                continue
            copies = l_frag[f_inx][0]
            start_pos = l_frag[f_inx][1]
            frag = l_frag[f_inx][2]
            frag_len = len(l_frag[f_inx][2])
            end_pos = start_pos + frag_len

            bo_selected = False
            for i in range(start_pos, end_pos):
                if len(l_hap_on_snps[i]) >= 6 and not bo_selected:
                    no_matches = 0
                    cnt_shared = 0
                    whichHap = -1
                    bo_overlapAll = True
                    for sh in l_hap_on_snps[i]:
                        bo_overlap, bo_conflict = self.overlap_conflict_Region(l_h[sh], f)
                        if not bo_overlap:
                            bo_overlapAll = False
                            break
                        if not bo_conflict:
                            cnt_shared += 1
                            whichHap = sh
                        if cnt_shared > 1:
                            break

                    if cnt_shared == 1 and bo_overlapAll:
                        bo_selected = True

                elif bo_selected and len(l_hap_on_snps[i]) >= 6:
                    if whichHap not in l_hap_on_snps[i]:
                        bo_selected = False
                        break

            if bo_selected:
                l_hap_absorb_frag[whichHap].append(f_inx)
                l_frag_used[f_inx] = True
        for inx_i, i in enumerate(l_hap_absorb_frag):
            for j in i:
                l_h[inx_i].merge_with_frag(l_frag[j], j)

        return

    def obtain_haps_from_connected_components(self, l_cc, l_blocks_sorted):  # phase_by_sorted_type
        l_ans = []
        for i in l_cc:
            l_sorted_haps = sorted(i, key=lambda x: x[0])
            n = l_sorted_haps[0]
            merged_hap = l_blocks_sorted[n[0]][n[1]].copy()

            for jj in range(1, len(l_sorted_haps)):
                n = l_sorted_haps[jj]
                merged_hap = merged_hap + l_blocks_sorted[n[0]][n[1]]
            l_ans.append(merged_hap)
        return l_ans

    def use_pairs_to_connect_blocks(self, l_h, l_f, l_f_used, scf_info):  # phase_by_sorted_type

        cnt_used_frag = 0
        lenn = scf_info.hap_len

        l_f_new = []
        cnt_frags = 0
        for f in l_f:
            l_f_new.append((f[0] * (len(f[2]) - f[2].count('-')), f[0], f[1], f[2], cnt_frags))
            cnt_frags += 1

        sorted_l_f_new = sorted(l_f_new, key=lambda x: (x[0]), reverse=True)
        sorted_l_f_new_used = [False] * len(sorted_l_f_new)
        bo_iter = True
        while bo_iter:
            bo_iter = False

            cnt_f = 0
            for f in sorted_l_f_new:
                if sorted_l_f_new_used[cnt_f]:
                    cnt_f += 1
                    continue

                l_roof = []
                for i in range(lenn):
                    l_roof.append([])
                cnt_h = 0
                for h in l_h:
                    for snp in range(h.length):
                        if h.haplotype[snp] != '-':
                            l_roof[snp + h.start_pos].append(cnt_h)
                    cnt_h += 1

                cop = f[1]
                start = f[2]
                frag = f[3]
                inx_frag = f[4]
                l_cand = []
                for snp in range(len(frag)):
                    if frag[snp] != '-':
                        l_cand += l_roof[snp + start]

                if len(l_cand) == 0:
                    continue
                l_cand_uniq = list(set(l_cand))
                l_sel = []
                for i in l_cand_uniq:
                    if mapped_to(frag, start, l_h[i]):
                        l_sel.append(i)
                if len(l_sel) == 0:
                    continue
                if len(l_sel) == 1:
                    continue
                    l_h[l_sel[0]].merge_with_frag((cop, start, frag), inx_frag)

                l_blocks = []
                for i in l_sel:
                    l_tmp = []
                    list_of_blocks(l_h[i].name, l_tmp)
                    l_blocks.append(l_tmp)

                no_conflict = 0
                no_no_conflict = 0
                bo = False
                for i in range(len(l_blocks) - 1):
                    for j in range(i + 1, len(l_blocks)):
                        if any([item in l_blocks[i] for item in l_blocks[j]]):
                            no_conflict += 1
                            bo = True
                        else:
                            no_no_conflict += 1
                            sel = (i, j)

                if no_no_conflict == 1:
                    start_bl = min(l_h[l_sel[sel[0]]].start_pos, l_h[l_sel[sel[1]]].start_pos)

                    new_hap = l_h[l_sel[sel[0]]] + l_h[l_sel[sel[1]]]
                    new_hap.frag_list_to_merge_later.append(f)

                    l_h[l_sel[sel[0]]].copy_from(new_hap)
                    del l_h[l_sel[sel[1]]]
                    sorted_l_f_new_used[cnt_f] = True
                    cnt_used_frag += 1
                    bo_iter = True

                cnt_f += 1

            cnt_iiii = 0
            for i in sorted_l_f_new_used:
                if i:
                    cnt_iiii += 1

        for h in l_h:
            if len(h.frag_list_to_merge_later) > 0:
                start_bl = h.start_pos
                for f in h.frag_list_to_merge_later:
                    if start_bl > f[2]:
                        start_bl = f[2]

                h.update_hap_by_stored_fragments()

        l_f = []
        l_f_used = []
        for i in range(len(sorted_l_f_new)):
            if not sorted_l_f_new_used[i]:
                f = sorted_l_f_new[i]
                l_f.append((f[1], f[2], f[3]))
                l_f_used.append(False)

        return

    def KruscalLike_or_coloring_algorithm(self):

        extention_allowed = True
        bo_connect_adjacent = True
        ploidy = self.ploidy
        typeLenMin = 2
        typeLenMax = self.max_type_len
        WinLen = self.WinLen
        l_frag_single = self.single_frag
        l_frag_pair = self.pair_frag

        l_frag_single_used = np.zeros(len(l_frag_single), dtype= bool)
        l_frag_pair_used = len(l_frag_pair) * [False]

        scf_hap_len = self.scf_hap_len

        if len(l_frag_single) == 0:
            return [], []

        l_all_types = self.obtain_types(l_frag_single, typeLenMin, typeLenMax, ploidy,WinLen)
        #l_all_types = self.obtain_types_nolimit(l_frag_single, typeLenMin, typeLenMax, ploidy)

        if len(l_all_types) == 0:
            return [], []

        l_haps_on_blocks, l_hap2hap, l_blocks_sorted_extenstion = self.phase_type_region(l_all_types, l_frag_single,
                                                                                         l_frag_single_used,
                                                                                         scf_hap_len, extention_allowed,
                                                                                         ploidy)
        l_cc = self.make_block_graph_alg2(l_blocks_sorted_extenstion, l_frag_pair, l_frag_pair_used, l_hap2hap, ploidy)
        l_obtained_haps = self.obtain_haps_from_connected_components(l_cc, l_blocks_sorted_extenstion)
        return l_haps_on_blocks, l_obtained_haps

        gap_clossing = False
        if gap_clossing:
            l_frag_all = l_frag_single + l_frag_pair
            l_frag_all_used = l_frag_single_used + l_frag_pair_used

            self.filling_gaps(scf_info, l_blocks_sorted, l_frag_all, l_frag_all_used, l_obtained_haps)
            pass

    def KruscalLike_or_coloring_algorithm_single(self):

        extention_allowed = True
        bo_connect_adjacent = True
        ploidy = self.ploidy
        typeLenMin = 2
        typeLenMax = self.max_type_len
        WinLen = self.WinLen
        l_frag_single = self.single_frag
        l_frag_single_used = np.zeros(len(l_frag_single),dtype=bool)
        scf_hap_len = self.scf_hap_len

        if len(l_frag_single) == 0:
            self.l_haps_single = []
            self.l_hap2hap = []
            self.l_blocks_sorted_extenstion = []
            return

        l_all_types = self.obtain_types(l_frag_single, typeLenMin, typeLenMax, ploidy, WinLen)

        if len(l_all_types) == 0:
            self.l_haps_single = []
            self.l_hap2hap = []
            self.l_blocks_sorted_extenstion = []
            return

        # self.l_haps_single, self.l_hap2hap, self.l_blocks_sorted_extenstion = self.phase_type_region(l_all_types,
        #                                                                                              l_frag_single,
        #                                                                                              l_frag_single_used,
        #                                                                                              scf_hap_len,
        #                                                                                              extention_allowed,
        #                                                                                              ploidy)
        self.l_haps_single = self.phase_type_region(l_all_types,
                                                                                                     l_frag_single,
                                                                                                     l_frag_single_used,
                                                                                                     scf_hap_len,
                                                                                                     extention_allowed,
                                                                                                     ploidy)

        return

    def KruscalLike_or_coloring_algorithm_pair(self):
        ploidy = self.ploidy
        l_frag_pair = self.l_frag_pair_uniq
        l_frag_pair_used = len(l_frag_pair) * [False]

        scf_hap_len = self.scf_hap_len

        l_cc = self.make_block_graph_alg2(self.l_blocks_sorted_extenstion, l_frag_pair, l_frag_pair_used,
                                          self.l_hap2hap, ploidy)
        self.l_haps_pair = self.obtain_haps_from_connected_components(l_cc, self.l_blocks_sorted_extenstion)
        return

    def KruscalLike_or_coloring_algorithm_gap_close(self):
        return
        gap_clossing = False
        if gap_clossing:
            l_frag_all = l_frag_single + l_frag_pair
            l_frag_all_used = l_frag_single_used + l_frag_pair_used

            self.filling_gaps(scf_info, l_blocks_sorted, l_frag_all, l_frag_all_used, l_obtained_haps)
            pass

    def run_method(self):
        if self.algName in ['kruskal', 'coloring']:
            l1, l2 = self.KruscalLike_or_coloring_algorithm()
            return l1, l2

    def run_method_(self):
        if self.algName in ['kruskal', 'coloring']:
            if self.algorithm_mode == 'single':
                return self.KruscalLike_or_coloring_algorithm_single()
            elif self.algorithm_mode == 'pair':
                return self.KruscalLike_or_coloring_algorithm_pair()
            elif self.algorithm_mode == 'gap close':
                return self.KruscalLike_or_coloring_algorithm_gap_close()
