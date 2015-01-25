import os
import sys
import urllib2
import hashlib 
from utils import error

cohort = "ACC BLCA BRCA CESC CHOL COAD COADREAD DLBC ESCA FPPP GBM GBMLGG HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PANCANCER PANCAN8 PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM"
cohort = cohort.split(" ")

dryrun = True 

for opt in [x.lower() for x in sys.argv[1:]]:
    print "Option: ", opt
    for d in cohort:
        if opt in ['mrnaseq', 'clinical', 'mutation']:
            root_url = 'http://gdac.broadinstitute.org/runs/stddata__2014_12_06/data/'
            root_url += d + "/20141206/gdac.broadinstitute.org_" + d
            outdir = './stddata__2014_12_06/' + d + '/20141206/'
        if opt in ['gistic2']:
            root_url = "http://gdac.broadinstitute.org/runs/analyses__2014_10_17/data/"
            root_url += d + "/20141017/gdac.broadinstitute.org_" + d 
            outdir = "./analyses__2014_10_17/" + d + "/20141017/"
        
        if opt == 'gistic2':
            url = root_url + "-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0.tar.gz"
        elif opt == 'mrnaseq':
            url = root_url + '.mRNAseq_Preprocess.Level_3.2014120600.0.0.tar.gz'
        elif opt == 'clinical':
            url = root_url + '.Merge_Clinical.Level_1.2014120600.0.0.tar.gz'
        elif opt == 'mutation':
            url = root_url + '.Mutation_Packager_Calls.Level_3.2014120600.0.0.tar.gz'
        else:
            error('Unknown option')

        outf = outdir + os.path.basename(url)
        if os.stat(outf).st_size != 0:
            print d + " already downloaded. Skipping.."
        else:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            try:
                gdac = urllib2.urlopen(url)

                if not dryrun:
                    print "Started ", d
                    with open(outf, 'wb') as f:
                        f.write(gdac.read())
                    # File downloaded check MD%
                    with open(outf, 'r') as f:
                        web_md5 = urllib2.urlopen(url+'.md5').read().split()[0]
                        local_md5 = hashlib.md5(f.read()).hexdigest()
                        if web_md5 != local_md5:
                            error("*****md5sum MISMATCH *****")
                    print "Finished ", d
                    
                else:
                    print "Will download: ", url
                    print "Will save at: ", outf
            except urllib2.HTTPError:
                print "\n****** Failed ***** ", d
                print url, "\n"
                print "*"*20