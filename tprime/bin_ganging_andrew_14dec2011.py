#! /usr/bin/env python

#python -i python/tree.py --var=fitMass --input=limit --limit

# =====================================================
#  BIN GANGING ALGORITHM
# =====================================================

# Import everything from ROOT
from ROOT import *
gROOT.Macro("rootlogon.C")

# =============== 
# options
# ===============
from optparse import OptionParser
parser = OptionParser()	  				  				  					  			  					  				  				  			  
parser.add_option('--input', metavar='T', type='string', action='store',
                  default='ele-2D',
                  dest='input',
                  help='input file')
(options,args) = parser.parse_args()
# ==========end: options =============

f = TFile(options.input+".root")

histos = []

histos.append (f.Get("BKG_fitMass_ht35") )
histos.append (f.Get("Top_fitMass_ht35") )
hsm = histos[0].Clone()
hsm.Add(histos[1])

histos.append (f.Get("Tprime350_fitMass_ht35") )
histos.append (f.Get("Tprime400_fitMass_ht35") )
histos.append (f.Get("Tprime450_fitMass_ht35") )
histos.append (f.Get("Tprime500_fitMass_ht35") )
histos.append (f.Get("Tprime550_fitMass_ht35") )
histos.append (f.Get("Tprime600_fitMass_ht35") )

histos.append (f.Get("BKG_JES095_fitMass_ht35") )
histos.append (f.Get("Top_JES095_fitMass_ht35") )
histos.append (f.Get("Tprime350_JES095_fitMass_ht35") )
histos.append (f.Get("Tprime400_JES095_fitMass_ht35") )
histos.append (f.Get("Tprime450_JES095_fitMass_ht35") )
histos.append (f.Get("Tprime500_JES095_fitMass_ht35") )
histos.append (f.Get("Tprime550_JES095_fitMass_ht35") )
histos.append (f.Get("Tprime600_JES095_fitMass_ht35") )


histos.append (f.Get("BKG_JES105_fitMass_ht35") )
histos.append (f.Get("Top_JES105_fitMass_ht35") )
histos.append (f.Get("Tprime350_JES105_fitMass_ht35") )
histos.append (f.Get("Tprime400_JES105_fitMass_ht35") )
histos.append (f.Get("Tprime450_JES105_fitMass_ht35") )
histos.append (f.Get("Tprime500_JES105_fitMass_ht35") )
histos.append (f.Get("Tprime550_JES105_fitMass_ht35") )
histos.append (f.Get("Tprime600_JES105_fitMass_ht35") )

histos.append (f.Get("Data_fitMass_ht35") )

c = TCanvas("c","c",800,600)

gStyle.SetOptStat(0)

NbinsX = hsm.GetNbinsX()
NbinsY = hsm.GetNbinsY()


# flag zero bins, they need to be migrated as well
for binX in range (1,NbinsX+1) :
	for binY in range (1, NbinsY+1) :		
		if ( hsm.GetBinContent(binX, binY) == 0 ) :
			hsm.SetBinError(binX, binY, 1. )


hrel = hsm.Clone()
hrel.SetMaximum(1.0)

while True :

# find a bin with the maximum relative error
	max_rel_error = 0.4
	bin_maxX = -1 
	bin_maxY = -1
 
	for binX in range (1,NbinsX+1) :
		for binY in range (1, NbinsY+1) :
			if (( hsm.GetBinContent(binX, binY) == 0 ) and ( hsm.GetBinError(binX, binY) == 0 )) :
				continue
			elif (( hsm.GetBinContent(binX, binY) == 0 ) and ( hsm.GetBinError(binX, binY) == 1. )) :			 	
				rel_error = 0.41
			else :		
				rel_error = hsm.GetBinError(binX, binY) / hsm.GetBinContent(binX, binY) 
			if ( rel_error > max_rel_error ) :
				max_rel_error = rel_error
				bin_maxX = binX
				bin_maxY = binY

	if (( bin_maxX < 0 ) and ( bin_maxY < 0 )) : 
		break
				
#	print "\n bin_maxX = ", bin_maxX, "; bin_maxY = ", bin_maxY, "; rel_error = ", max_rel_error


# find an adjacent bin with the smallest relative error
	min_rel_error = max_rel_error
	bin_minX = bin_maxX
	bin_minY = bin_maxY

	for offset in range (1,50) :


		for binX in range (bin_maxX-offset, bin_maxX+offset) :
			for binY in range (bin_maxY-offset, bin_maxY+offset) :
										
				if ( hsm.GetBinContent(binX, binY) == 0 ) :
					continue
				if ( abs(binX - bin_maxX) + abs(binY - bin_maxY ) > offset) :
					continue	

				rel_error = hsm.GetBinError(binX, binY) / hsm.GetBinContent(binX, binY) 	
#				print "(",binX,",",binY, ") = ", rel_error
				if ( rel_error < min_rel_error ) :
					min_rel_error = rel_error
					bin_minX = binX
					bin_minY = binY

#		print "\n found an adjacent bin, binX = ", bin_minX, "; binY = ", bin_minY, "; rel_error = ", min_rel_error
			
		if  ( (bin_minX == bin_maxX) and ( bin_minY == bin_maxY) ) :
#			print "\n This is the same bin ! Increasing offset to ", offset+1
			pass
		else :	
			break

	if  ( (bin_minX == bin_maxX) and ( bin_minY == bin_maxY) ) :
#		print "\n ERROR ! This is the same bin ! "
		break
		

	if ( hsm.GetBinContent(bin_maxX,bin_maxY) == 0 ) :
		hsm.SetBinError(bin_maxX,bin_maxY, 0.)
			
	hsm.SetBinContent(bin_minX, bin_minY, hsm.GetBinContent(bin_minX, bin_minY) + hsm.GetBinContent(bin_maxX, bin_maxY) ) 			
	hsm.SetBinError(bin_minX, bin_minY, TMath.Sqrt( TMath.Power(hsm.GetBinError(bin_minX, bin_minY),2) + TMath.Power( hsm.GetBinError(bin_maxX, bin_maxY) , 2)) ) 			
	hsm.SetBinContent(bin_maxX, bin_maxY, 0. ) 					
	hsm.SetBinError(bin_maxX, bin_maxY, 0. )
	
	for hist in range(0,len(histos)) :
		histos[hist].SetBinContent(bin_minX, bin_minY, histos[hist].GetBinContent(bin_minX, bin_minY) + histos[hist].GetBinContent(bin_maxX, bin_maxY) ) 			
		histos[hist].SetBinError(bin_minX, bin_minY, TMath.Sqrt( TMath.Power(histos[hist].GetBinError(bin_minX, bin_minY),2) + TMath.Power( histos[hist].GetBinError(bin_maxX, bin_maxY) , 2)) ) 			
		histos[hist].SetBinContent(bin_maxX, bin_maxY, 0. ) 					
		histos[hist].SetBinError(bin_maxX, bin_maxY, 0. )
	
	for binX in range (1,NbinsX+1) :
		for binY in range (1, NbinsY+1) :
			if (( hsm.GetBinContent(binX, binY) == 0 ) and ( hsm.GetBinError(binX, binY) == 0 )) :
				rel_error = 0.0
			elif (( hsm.GetBinContent(binX, binY) == 0 ) and ( hsm.GetBinError(binX, binY) == 1. )) :			 	
				rel_error = 0.41
			else :		
				rel_error = hsm.GetBinError(binX, binY) / hsm.GetBinContent(binX, binY) 
			hrel.SetBinContent(binX, binY, rel_error )
#			print "\n(",binX,",",binY, ") = ", rel_error
		
	hrel.Draw("COLZ")
	c.Update()



# cross-check	
for i in range (0, len(histos) ) :
	for binX in range (1,NbinsX+1) :
		for binY in range (1, NbinsY+1) :	
			if (( hsm.GetBinContent(binX, binY) == 0 ) and ( histos[i].GetBinContent(binX,binY) != 0 )) :
				print "ERROR! hsm.GetBinContent(",binX,",", binY,") = ", hsm.GetBinContent(binX, binY), " ! = histos[",i,"].GetBinContent(",binX,",",binY,") = ", histos[i].GetBinContent(binX,binY)

	
#	print "\n\n\n\n\n\n"
num_nonzero = 0
for binX in range (1,NbinsX+1) :
	for binY in range (1, NbinsY+1) :	
		if ( hsm.GetBinContent(binX, binY) != 0 ) :
			num_nonzero = num_nonzero + 1
print "num_nonzero = ", num_nonzero


fout = TFile(options.input+"_ganged.root","RECREATE")
for i in range (0, len(histos) ) :
	histos[i].Write()
fout.Close()


#1-dimensional histograms
fout1 = TFile(options.input+"_ganged-1D.root","RECREATE")
hists = []
for i in range (0, len(histos) ) :
	hist = TH1D(histos[i].GetName(), histos[i].GetTitle(),num_nonzero,0,num_nonzero)
	bin = 1
	for binX in range (1,NbinsX+1) :
		for binY in range (1, NbinsY+1) :	
			if ( hsm.GetBinContent(binX, binY) != 0 ) :
				hist.SetBinContent(bin, histos[i].GetBinContent(binX,binY))
				hist.SetBinError(bin, histos[i].GetBinError(binX,binY))
				bin = bin + 1
	hists.append(hist)			
	hists[i].Write()
fout1.Close()	
						
				
				
				
				
				
				
				
				
				
				
				
				
				
