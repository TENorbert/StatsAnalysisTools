#!/usr/bin/env python
#############################
#
# fittable_generator.py
#     puts fit parmaters from a log-file into a latex table
#
#  Author : Michael Luk
#  Date   : 19 Jan 2012
#
#  e.g. ./fittable_generator.py -i StdHypo_LogFile.txt -o outputFileName -c channelName
#
#############################

###
#
# banner
#
###########
def banner():
        print '''
        +--------------------------------------------------------------
        |
        | fittable_generator.py
        |
        | Generates table of fit values for nuisance parameters
        | from HyotestInverter log file
        |
        | author: Michael Luk, Jan 2012
        |
        +--------------------------------------------------------------
            '''
banner()
_legend = '[fittable_generator]:'
        

import os

from optparse import OptionParser
add_help_option = "./fittable_generator -ACTION [other options]"
parser = OptionParser(add_help_option)

parser.add_option("-i", "--in-file", dest="in_file", default=None,
                                    help="Input file (.txt) to generate fittable"
                                        )
parser.add_option("-o", "--out-file", dest="out_file", default=None,
                  help="output filename (w/out .tex) to generate fittable" )
parser.add_option("-c", "--channel", dest="channel", default=None,
                  help="output filename (w/out .tex) to generate fittable")
print _legend,'parsing command line options...',
(options, args) = parser.parse_args()
print _legend,'done'


channelname    = str(options.channel) #"combined"
inputfilename  = str(options.in_file)  #"comb_log2.txt"  #log file of ./plots --toys limit > logfile.txt
output_textable= str(options.out_file) + ".tex" #"combined_latextable.tex" # the texfile name to make the table
outputdir      = 'tabletemp'      # temporory storage required for table making - deleted at end
if not os.path.exists('./'+outputdir):
    os.mkdir('./'+outputdir)
    print _legend,'created temporary folder :',outputdir

nuisanceParamNames = open(outputdir+'/nuisanceParamNames.txt','w')
nuisanceAmpersand  = open(outputdir+'/nuisanceAmpersand.txt','w')
nuisanceLine       = open(outputdir+'/nuisanceLine.txt','w')
nuisanceBlock = False

infile        = open(inputfilename)
print _legend, 'opened input file :', inputfilename
deltam        = 25
massinit      = 400

filecounter   = massinit - deltam
linecounter=0
for line in infile:
    if nuisanceBlock:
        nuisanceTempFile.write(line)

        if filecounter==massinit and ("LGM" not in line):
            nameBegin = line.find('"')
            nameEnd   = line[nameBegin+1:].find('"')
            linecounter+=1
            name = line[nameBegin+1:nameEnd+nameBegin+1].replace("_"," ") 
            nuisanceParamNames.write(name+" & \n")
            nuisanceAmpersand.write("&\n")
            nuisanceLine.write("\\\\ \n")
            
        
    if "LGMMID" in line:
        nuisanceBlock= True
        filecounter  = filecounter + deltam
        nuisanceTempFile = open(outputdir+'/log_temp_'+str(filecounter)+'.txt','w')
    if "LGMEND" in line:
        nuisanceBlock= False
        nuisanceTempFile.close()
        nuisanceAmpersand.close()
        nuisanceLine.close()
        nuisanceParamNames.close()
        os.system("awk '{printf("+'"%.4f\\n"'+",$6)}' "+outputdir+'/log_temp_'+str(filecounter)+'.txt  > '+outputdir+'/tempawk_'+str(filecounter)+'.txt')
        

mass = massinit
buildcommand = "paste "

tableheader= '\\documentclass[sort&compress,onecolumn,letterpaper,openany]{revtex4} \n \\usepackage[UKenglish]{babel} \n \\begin{document} \n \\pagenumbering{roman}\n \\begin{table} \n \\caption{Shows fitted values of the parameters for the '+channelname+' channel} \n \\begin{tabular}{'
masstitle  = '\hline\hline \n $t^\prime$ mass & '
while mass <= filecounter:
    tableheader+='c'
    if mass == massinit:
        buildcommand += outputdir+'/nuisanceParamNames.txt '
    if mass!=filecounter:
        masstitle += str(mass)+"GeV & "
        buildcommand += outputdir+'/tempawk_'+str(mass)+'.txt '+outputdir+'/nuisanceAmpersand.txt '
    else:
        tableheader+="c} \n"
        masstitle += str(mass)+"GeV \\\\ \n \hline \n"
        buildcommand += outputdir+'/tempawk_'+str(mass)+'.txt '+outputdir+'/nuisanceLine.txt '        
    mass += deltam

table = open(output_textable,'w')
table.write(tableheader+masstitle)
os.system(buildcommand+" > "+outputdir+'/temptable.txt')
temptable = open(outputdir+'/temptable.txt')

linecounter2=0
for line in temptable:
    linecounter2+=1
    if linecounter>=linecounter2:
        table.write(line)
table.write("\\hline\\hline \n \\end{tabular} \n \\end{table} \n\n \\end{document}")
table.close()


#closing and deleting temporary files
print _legend, 'clean up, closing files and deleting temporary directory'
nuisanceParamNames.close()
nuisanceAmpersand.close()
nuisanceLine.close()
temptable.close()

print _legend, 'finished, output is :',output_textable
os.system('rm -r '+outputdir+"/*")
os.system('rm -r '+outputdir)
