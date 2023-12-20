import sys
import os
import pandas as pd
import jinja2
import re

'''
This scripts completes the html template with the
specific information of each run.

USAGE:  python3 CAPTVRED_create_report.py     \
                samples_definition.tbl        \      
                reports_dir                   \       
                run_ID                        \ 
                samples_summary_tbl           \       
                outfl                                

'''

#define variables and preferences:
myclasstr="mytables table  dataTables_wrapper table-striped table-hover table-borderless"
idir=sys.argv[2]
runid=sys.argv[3]
taxap=sys.argv[7].lower()
covfigs_dir=idir + "/coverage_figures"  #Coverage figures directory

## 1. Charge information:

#samples:
sd=pd.read_csv(sys.argv[1],sep='\t', skip_blank_lines=True, comment='#',
               names=[ 'Sample ID', 'Sequencing ID', 'Sample\nName', 'Reads\nType', 'Sample\nSource', 'Experimental\nMethodology'  ])   
sd_html=sd.to_html(table_id="tblsd", 
                   border = 0, 
                   classes="mytables table  dataTables_wrapper table-striped table-hover table-borderless"
                   ).replace('<thead>','<thead class="thead-dark">')

sd = sd.reset_index()  # make sure indexes pair with number of rows
samplesls=[]

#Sumary table:

sumtbl=pd.read_csv(sys.argv[4],sep=',')   
sumtbl_html=sumtbl.to_html(table_id="summarytbl", 
                   border = 0, 
                   classes="mytables table  dataTables_wrapper table-striped table-hover table-borderless"
                   ).replace('<thead>','<thead class="thead-dark">')
#Taxonomy tables:
for index, row in sd.iterrows():
  idf=row[1]
  aid="sample"+str(index)+"_tablink"
  href="#sample"+str(index)+"_tabpanel"
  ariactrls="sample"+str(index)+"_tabpanel"
  dividsp="sample"+str(index)+"_sp"
  dividsq="sample"+str(index)+"_sq"
  dividrd="sample"+str(index)+"_rd"
  bysp=pd.read_csv( idir + "/" + idf + "." + taxap + "_taxonomysum_byspecie.tbl" , sep="\t" )
  bysp_html=bysp.to_html(  table_id="sptbl_".index, 
                           border = 0, 
                           classes=myclasstr
                        ).replace('<thead>','<thead class="thead-dark">')
                        
  bysq=pd.read_csv( idir + "/" + idf + "." + taxap +"_taxonomysum_byspecie.tbl", sep="\t" )
  bysq_html=bysq.to_html(  table_id="sqtbl_".index, 
                           border = 0, 
                           classes=myclasstr
                        ).replace('<thead>','<thead class="thead-dark">')
  byrd=pd.read_csv( idir + "/" + idf + "." + taxap + "_taxonomysum_byread.tbl", sep="\t" )
  byrd_html=bysp.to_html(  table_id="rdtbl_".index, 
                           border = 0, 
                           classes=myclasstr
                        ).replace('<thead>','<thead class="thead-dark">')
  ## Nota pel futur: Aixoò millor amb un hash! 
  samplesls.append( [idf, index, aid, href, ariactrls, dividsp, bysp_html, dividsq, bysq_html, dividrd, byrd_html] ) 

print(samplesls[1][1], samplesls[1][0], samplesls[1][2], samplesls[1][3])
print(samplesls[1][4], samplesls[1][5])
print(samplesls[1][7], samplesls[1][9])


## Coverage Figures:
images={}
for index, row in sd.iterrows():
   sampid=row[1]
   print(sampid)
   pattern=re.compile('.*(%s).*png$'%sampid)
   images[sampid] = [ '/'.join([idir, "coverage_figures", file]) for file in os.listdir(covfigs_dir) if pattern.match(file)]
   #print (images)


## JINJA:
thyfl   = sys.argv[5]
tpldir  = os.path.dirname(thyfl)
tplname = os.path.basename(thyfl)
templateLoader = jinja2.FileSystemLoader(searchpath=tpldir)
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(tplname)
outhtml = template.render( run_name=runid, 
                           SAMPLESDEF=sd_html, 
                           samples=samplesls,
                           sumtbl=sumtbl_html,
                           nxfrep=idir + "/Nextflow_execution_report.html",
                           mqc_raw=idir + "/" + runid + "_multiqc_raw.html",
                           mqc_filt=idir + "/" + runid + "_multiqc_filt.html",
                           covfigures=images
                           )

with open(sys.argv[6], 'w') as f:
    f.write(outhtml)


####
