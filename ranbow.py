#!/usr/bin/env python2.7

__version__ = '2.0'
__author__ = 'M-Hossein Moeinzadeh'

import argparse
import sys


def main(argv):


    options = ['hap', 'eval', 'sim', 'hscaf', 'close_gap', 'util', 'phylo', 'circlarHeatMap', 'repeat_vs_tree','addSMRT','version']
    if len(argv) == 0 or argv[0] not in options:
        for i in open("options.txt"):
            print i.rstrip()

        exit()

    else:
        main_act = argv[0]

    if main_act == 'version':
        print 'Ranbow version 2.0'

    if main_act == 'hap':
        import haplotyper
        hap = haplotyper.haplotyper(argv, main_act)





        hap.is_eval = False
        hap.is_addSMRT = False
        if hap.act in ['hap', 'index']:
            if hap.okey_params:
                #print '___________________', hap.act, '____________', hap.okey_params
                hap.go()
        elif hap.act == 'collect':
            print 'number of folders:', hap.no_processor
            import subprocess
            import os

            print subprocess.check_output(
                ['bash', os.path.dirname(sys.argv[0])+'/hap_collect.sh', hap.outputFolder, hap.prefix, str(hap.no_processor), hap.ref_fai])

            #import IGV
            #igv = IGV.IGV(argv,hap.outputFolder, hap.prefix)
            #igv.makeBam_ATCG()
            #igv.makeSam()


        elif hap.act == 'modVCF':
            import vcf_correction
            v = vcf_correction.vcf_correction(hap)
            v.run()


    if main_act == 'eval':
        import haplotyper
        hap = haplotyper.haplotyper(argv, main_act)
        if hap.act == 'index':
            exit()

        if hap.act == 'run':
            hap.go()

        if hap.act == 'collect':
            print 'number of folders:', hap.no_processor
            print hap.evalOutput
            import subprocess
            import os
            print subprocess.check_output(['bash', os.path.dirname(sys.argv[0])+'/eval_collect.sh', hap.evalOutput])

    if main_act == 'addSMRT':
        import haplotyper
        hap = haplotyper.haplotyper(argv, main_act)
        if hap.act == 'index':
            exit()

        if hap.act == 'run':
            hap.go()

        if hap.act == 'collect':
            print 'number of folders:', hap.no_processor
            print hap.evalOutput
            import subprocess
            print subprocess.check_output(['bash', './eval_collect.sh', hap.evalOutput])



    if main_act == 'phylo':
        import phylo
        p = phylo.phylo(argv)

        if p.act == 'index':
            p.pre_run()

        if p.act == 'run':
            p.run()

        if p.act == 'collect':
            import subprocess
            out = subprocess.check_output(
                ['bash', './phylo_collect.sh', p.phyloFolder, str(p.noCpu - 1), str(p.mnHapLen)])

            def compliment(i):
                ii = i[::-1]
                ii = ii.replace('(', '-')
                ii = ii.replace(')', '(')
                return ii.replace('-', ')')

            a = out.split()
            mp = {}
            for j in xrange(1, len(a), 2):
                c_j = compliment(a[j])
                c_j = c_j if c_j < a[j] else a[j]
                mp[c_j] = mp.get(c_j, 0) + int(a[j - 1])
            for i in sorted(mp):
                print mp[i], i

        if p.act == 'tree':
            p.tree_char(p.phyloFolder + 'collect.txt')

    if main_act == "sim":  # simulate simulation data
        p = project(argv)
        p.run_simulation_myMethod_differentLen()

    if main_act == "hscaf":  # haplo scaffolder
        import haploScaffolder
        hs_p = haploScaffolder.haploScaffolder(argv)
        pass

    if main_act == "close_gap":  # haplo scaffolder
        import gap_close
        hs_p = gap_close.gap_close(argv)
        pass

    if main_act == "circlarHeatMap":  # haplo scaffolder
        import circlarHeatMap
        c = circlarHeatMap.circlarHeatMap()
        # c.hap_percentage()
        # c.SNP_dencity()
        # c.repeat_density()
        # c.genoem_cov()
        # c.gene_density()
        # c.frame_shift_exome()
        c.frame_shift_transcript()
    # c.frame_shift_manual_check()

    if main_act == "repeat_vs_tree":
        import repeatVStree
        c = repeatVStree.repeatVStree()
        c.run()

    if main_act == "util":
        import util
        u = util.util(argv)
        if u.mode == 'allmap_hap_reg':
            u.allmapp_hap_region()


if __name__ == "__main__":
    main(sys.argv[1:])
