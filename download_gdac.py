import os
import sys
import urllib2

cohort = "ACC BLCA BRCA CESC CHOL COAD COADREAD DLBC ESCA FPPP GBM GBMLGG HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PANCANCER PANCAN8 PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM"
cohort = cohort.split(" ")
print cohort

for d in cohort:
    url = "http://gdac.broadinstitute.org/runs/analyses__2014_10_17/data/" + d + "/20141017/gdac.broadinstitute.org_" + d + "-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0.tar.gz"

    outdir = "./analyses__2014_10_17/" + d + "/20141017"

    outf = outdir + os.sep + "gdac.broadinstitute.org_" + d + "-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0.tar.gz" 
    if os.path.exists(outf):
        print d + " already downloaded. Skipping.."
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        print "Started ", d
        try:
            gdac = urllib2.urlopen(url)
            with open(outf, 'wb') as f:
                f.write(gdac.read())

            print "Finished ", d
        except urllib2.HTTPError:
            print "****** Failed ***** ", d

