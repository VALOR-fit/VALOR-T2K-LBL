#!/usr/bin/env python

#	 Brazilian plot  macro 
# 	Calculates brazilian dchi2 value at each marginalization point for each pseudo experiment (toy) i.e. dchi2[mpt][toy]
# 	Calculates one and two sigma intervals either for all toys ( one sided plot) 
# 	Or with respect to a median toy ( two-sided plot)
#		Usually we create it in T2K data BF point for normal or inverted mh
# 	But can also be created at any other point for a cross-check 
#	
# 	Maria Antonova <maria.antonova \at ific.uv.es>
# 	February 2021

import ROOT
import sys
import os
import math
from ROOT import gPad
sys.path.append(os.path.abspath(os.environ["VALOR_FMWK"] + "/src/valor/Utilities/"))

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")
valor = ROOT.valor
conv =  valor.conv
util =  valor.util

#dchi2  data fit histogram name 
_dchi2_dhname="chi2_1d_fit"
#creating a fancy golden yellow and green colours  
golden = ROOT.TColor.GetFreeColorIndex()
cg=ROOT.TColor(golden, 1,0.87,0)
bgreen = ROOT.TColor.GetFreeColorIndex()
cgr=ROOT.TColor(bgreen, 0,0.61,0.23)

_dchi2_arr_nh_name="Brazilian_dchi2_dcp_bin_vs_toy_vec_nh" #TObjArray storing chi[toy] vectors for each marginalization point of dcp 
_dchi2_arr_ih_name="Brazilian_dchi2_dcp_bin_vs_toy_vec_ih" #TObjArray storing chi[toy] vectors for each marginalization point of dcp 
#_______________________________________________________________________
def draw_one_sided(in_filenames, out_filename, data_fit_nh, data_fit_ih):
	print "Ready to generate one-sided plot? Off we go!"

	#combine all batches of marg toys in one
	nh_dchi_list, ih_dchi_list = combine(in_filenames)
	
	#Get data dchi distribution 
	fdatan = ROOT.TFile(data_fit_nh)
	hdata_nh = fdatan.Get(_dchi2_dhname)
	fdatai = ROOT.TFile(data_fit_ih)
	hdata_ih = fdatai.Get(_dchi2_dhname)

	#set +- 1/2 sigma histograms 
	dchi2_nh_p1sigma=hdata_nh.Clone("hdchi2_nh_p1sigma") 
	dchi2_ih_p1sigma=hdata_nh.Clone("hdchi2_ih_p1sigma") 

	dchi2_nh_p2sigma=hdata_nh.Clone("hdchi2_nh_p2sigma") 
	dchi2_ih_p2sigma=hdata_nh.Clone("hdchi2_ih_p2sigma") 

	#fill +- 1/2 sigma  histograms 
	print "Filling 1&2 sigma  histogram for each point"
	for ibin in range (0, len(nh_dchi_list)): 
		dchi2_nh_p1sigma.SetBinContent(ibin+1, calc_critical(nh_dchi_list[ibin], conv.kCL1Sigma))
		dchi2_ih_p1sigma.SetBinContent(ibin+1, calc_critical(ih_dchi_list[ibin], conv.kCL1Sigma))
		dchi2_nh_p2sigma.SetBinContent(ibin+1, calc_critical(nh_dchi_list[ibin], conv.kCL2Sigma))
		dchi2_ih_p2sigma.SetBinContent(ibin+1, calc_critical(ih_dchi_list[ibin], conv.kCL2Sigma))

	c_nh = ROOT.TCanvas("Brazilian_nh","Brazilian_nh",800,600)
	c_ih = ROOT.TCanvas("Brazilian_ih", "Brazilian_ih",800,600)
	
	#Plotting settings 
	fout_name = out_filename+".root"
	fout = ROOT.TFile(fout_name, "recreate")
	#NH
	dchi2_nh_p2sigma.GetYaxis().SetRangeUser(0, 25)
	dchi2_nh_p1sigma.GetYaxis().SetRangeUser(0, 25)
	dchi2_nh_p2sigma.SetLabelFont(42)
	dchi2_nh_p2sigma.SetLabelSize(0.05)
	dchi2_nh_p2sigma.GetXaxis().SetTitle("#delta_{CP} (Radians)")
	dchi2_nh_p2sigma.GetYaxis().SetTitle("#Delta #chi^{2}")
	dchi2_nh_p2sigma.GetXaxis().SetTitleFont(42)
	dchi2_nh_p2sigma.GetYaxis().SetTitleFont(42)
	dchi2_nh_p2sigma.GetXaxis().SetTitleSize(0.05)
	dchi2_nh_p2sigma.GetYaxis().SetTitleSize(0.05)
	dchi2_nh_p2sigma.GetYaxis().SetLabelSize(0.05)
	
	#IH
	dchi2_ih_p2sigma.GetYaxis().SetRangeUser(0, 40)
	dchi2_ih_p1sigma.GetYaxis().SetRangeUser(0, 40)
	dchi2_ih_p2sigma.SetLabelFont(42)
	dchi2_ih_p2sigma.SetLabelSize(0.05)
	dchi2_ih_p2sigma.GetXaxis().SetTitle("#delta_{CP} (Radians)")
	dchi2_ih_p2sigma.GetYaxis().SetTitle("#Delta #chi^{2}")
	dchi2_ih_p2sigma.GetXaxis().SetTitleFont(42)
	dchi2_ih_p2sigma.GetYaxis().SetTitleFont(42)
	dchi2_ih_p2sigma.GetXaxis().SetTitleSize(0.05)
	dchi2_ih_p2sigma.GetYaxis().SetTitleSize(0.05)
	dchi2_ih_p2sigma.GetYaxis().SetLabelSize(0.05)


	#Set up a legend 
 	lBrazilianPlot_nh = ROOT.TLegend(0.17,0.64,0.37,0.92)
 	lBrazilianPlot_nh.SetFillColor(0)
 	lBrazilianPlot_nh.SetLineColor(0)
 	lBrazilianPlot_nh.SetTextFont(42)
 	lBrazilianPlot_nh.SetTextSize(0.045)
	lBrazilianPlot_nh.AddEntry(dchi2_nh_p1sigma,"68.27% of  pseudo-experiments","f")
	lBrazilianPlot_nh.AddEntry(dchi2_nh_p2sigma,"95.45% of  pseudo-experiments","f")
	lBrazilianPlot_nh.AddEntry(hdata_nh,"Data","l")

	#Set up a legend 
 	lBrazilianPlot_ih = ROOT.TLegend(0.17,0.64,0.37,0.92)
 	lBrazilianPlot_ih.SetFillColor(0)
 	lBrazilianPlot_ih.SetLineColor(0)
 	lBrazilianPlot_ih.SetTextFont(42)
 	lBrazilianPlot_ih.SetTextSize(0.045)
	lBrazilianPlot_ih.AddEntry(dchi2_ih_p1sigma,"68.27% of  pseudo-experiments","f")
	lBrazilianPlot_ih.AddEntry(dchi2_ih_p2sigma,"95.45% of  pseudo-experiments","f")
	lBrazilianPlot_ih.AddEntry(hdata_ih,"Data","l")
	
	#p2, yellow 
	dchi2_nh_p2sigma.SetFillColor(golden)
	dchi2_ih_p2sigma.SetFillColor(golden)
	dchi2_nh_p2sigma.SetLineColor(golden)
	dchi2_ih_p2sigma.SetLineColor(golden)
	#p1, green 
	dchi2_nh_p1sigma.SetFillColor(bgreen)
	dchi2_ih_p1sigma.SetFillColor(bgreen)
	dchi2_nh_p1sigma.SetLineColor(bgreen)
	dchi2_ih_p1sigma.SetLineColor(bgreen)

 	hdata_nh.SetLineStyle(1)
 	hdata_nh.SetLineColor(1)
 	hdata_nh.SetLineWidth(2)

 	hdata_ih.SetLineStyle(1)
 	hdata_ih.SetLineColor(1)
 	hdata_ih.SetLineWidth(2)

	#drawing and printing NH
	c_nh.cd()
	dchi2_nh_p2sigma.Draw("c")
	dchi2_nh_p1sigma.Draw("c same")
	hdata_nh.Draw("l same")
	lBrazilianPlot_nh.Draw("same")
	ROOT.gPad.RedrawAxis()
	nname = out_filename+"_one_sided_NH.pdf"
	c_nh.Print(nname)

	#drawing and printing NH
	c_ih.cd()
	dchi2_ih_p2sigma.Draw("c")
	dchi2_ih_p1sigma.Draw("c same")
	hdata_ih.Draw("l same")
	lBrazilianPlot_ih.Draw("same")
	ROOT.gPad.RedrawAxis()
	iname = out_filename+"_one_sided_IH.pdf"
	c_ih.Print(iname)

	dchi2_ih_p1sigma.Write()
	dchi2_ih_p2sigma.Write()
	dchi2_nh_p1sigma.Write()
	dchi2_nh_p2sigma.Write()
	c_nh.Write()
	c_ih.Write()

	fout.Write()
	fout.Close() 		
#_______________________________________________________________________
def draw_two_sided(in_filenames, out_filename, data_fit_nh, data_fit_ih):
	print "Ready to generate two-sided plot? Off we go!"
	
	#combine all batches of marg toys in one
	nh_dchi_list, ih_dchi_list = combine(in_filenames)	
	
	#Get data dchi distribution 
	fdatan = ROOT.TFile(data_fit_nh)
	hdata_nh = fdatan.Get(_dchi2_dhname)

	fdatai = ROOT.TFile(data_fit_ih)
	hdata_ih = fdatai.Get(_dchi2_dhname)

	#set histograms for 
	#p1 sigma
	dchi2_nh_p1sigma=hdata_nh.Clone("hdchi2_nh_p1sigma") 
	dchi2_ih_p1sigma=hdata_nh.Clone("hdchi2_ih_p1sigma") 
	#p2 sigma
	dchi2_nh_p2sigma=hdata_nh.Clone("hdchi2_nh_p2sigma") 
	dchi2_ih_p2sigma=hdata_nh.Clone("hdchi2_ih_p2sigma") 
	#m1 sigma
	dchi2_nh_m1sigma=hdata_nh.Clone("hdchi2_nh_p1sigma") 
	dchi2_ih_m1sigma=hdata_nh.Clone("hdchi2_ih_p1sigma") 
	#m2 sigma 
	dchi2_nh_m2sigma=hdata_nh.Clone("hdchi2_nh_p2sigma") 
	dchi2_ih_m2sigma=hdata_nh.Clone("hdchi2_ih_p2sigma") 

	#median toy 
	median_nh=hdata_nh.Clone("median_toy_nh") 
	median_ih=hdata_nh.Clone("median_toy_ih") 
	
	#Filling up 1&2 sigma histos 

	for ibin in range (0, len(nh_dchi_list)): 
		s1_range =[]
		s2_range =[]
		s1_range, median_nh_mpt = calc_critical_binned(nh_dchi_list[ibin], conv.kCL1Sigma)
		s2_range, median_nh_mpt = calc_critical_binned(nh_dchi_list[ibin], conv.kCL2Sigma)
		dchi2_nh_p1sigma.SetBinContent(ibin+1, s1_range[1])
		dchi2_nh_p2sigma.SetBinContent(ibin+1, s2_range[1])
		dchi2_nh_m1sigma.SetBinContent(ibin+1, s1_range[0])
		dchi2_nh_m2sigma.SetBinContent(ibin+1, s2_range[0])
		median_nh.SetBinContent(ibin+1, median_nh_mpt)

		s1_range, median_ih_mpt = calc_critical_binned(ih_dchi_list[ibin], conv.kCL1Sigma)
		s2_range, median_ih_mpt = calc_critical_binned(ih_dchi_list[ibin], conv.kCL2Sigma)
		dchi2_ih_p1sigma.SetBinContent(ibin+1,s1_range[1])
		dchi2_ih_p2sigma.SetBinContent(ibin+1,s2_range[1])
		dchi2_ih_m1sigma.SetBinContent(ibin+1, s1_range[0])
		dchi2_ih_m2sigma.SetBinContent(ibin+1, s2_range[0])
		median_ih.SetBinContent(ibin+1, median_ih_mpt)


	
	c_nh = ROOT.TCanvas("Brazilian_nh","Brazilian_nh",800,600)
	c_ih = ROOT.TCanvas("Brazilian_ih", "Brazilian_ih",800,600)
	
	#Plotting settings 
	#NH
	fout_name = out_filename+".root"
	fout = ROOT.TFile(fout_name, "recreate")
	dchi2_nh_p2sigma.GetYaxis().SetRangeUser(0, 30)
	dchi2_nh_p1sigma.GetYaxis().SetRangeUser(0, 30)
	dchi2_nh_m2sigma.GetYaxis().SetRangeUser(0, 30)
	dchi2_nh_m1sigma.GetYaxis().SetRangeUser(0, 30)
	dchi2_nh_p2sigma.SetLabelFont(42)
	dchi2_nh_p2sigma.SetLabelSize(0.05)
	dchi2_nh_p2sigma.GetXaxis().SetTitle("#delta_{CP} (Radians)")
	dchi2_nh_p2sigma.GetYaxis().SetTitle("#Delta #chi^{2}")
	dchi2_nh_p2sigma.GetXaxis().SetTitleFont(42)
	dchi2_nh_p2sigma.GetYaxis().SetTitleFont(42)
	dchi2_nh_p2sigma.GetXaxis().SetTitleSize(0.05)
	dchi2_nh_p2sigma.GetYaxis().SetTitleSize(0.05)
	dchi2_nh_p2sigma.GetYaxis().SetLabelSize(0.05)
	#IH
	dchi2_ih_p2sigma.GetYaxis().SetRangeUser(0, 40)
	dchi2_ih_p1sigma.GetYaxis().SetRangeUser(0, 40)
	dchi2_ih_m2sigma.GetYaxis().SetRangeUser(0, 40)
	dchi2_ih_m1sigma.GetYaxis().SetRangeUser(0, 40)
	dchi2_ih_p2sigma.SetLabelFont(42)
	dchi2_ih_p2sigma.SetLabelSize(0.05)
	dchi2_ih_p2sigma.GetXaxis().SetTitle("#delta_{CP} (Radians)")
	dchi2_ih_p2sigma.GetYaxis().SetTitle("#Delta #chi^{2}")
	dchi2_ih_p2sigma.GetXaxis().SetTitleFont(42)
	dchi2_ih_p2sigma.GetYaxis().SetTitleFont(42)
	dchi2_ih_p2sigma.GetXaxis().SetTitleSize(0.05)
	dchi2_ih_p2sigma.GetYaxis().SetTitleSize(0.05)
	dchi2_ih_p2sigma.GetYaxis().SetLabelSize(0.05)

	#p2 yellow 
	dchi2_nh_p2sigma.SetFillColor(golden)
	dchi2_ih_p2sigma.SetFillColor(golden) 
	dchi2_nh_p2sigma.SetLineColor(golden)
	dchi2_ih_p2sigma.SetLineColor(golden) 

	#p1, green 
	dchi2_nh_p1sigma.SetFillColor(bgreen)
	dchi2_ih_p1sigma.SetFillColor(bgreen)
	dchi2_nh_p1sigma.SetLineColor(bgreen)
	dchi2_ih_p1sigma.SetLineColor(bgreen)

	#m1 Should be yellow to make a "band effect"
	dchi2_nh_m1sigma.SetFillColor(golden)
	dchi2_ih_m1sigma.SetFillColor(golden)
	dchi2_nh_m1sigma.SetLineColor(golden)
	dchi2_ih_m1sigma.SetLineColor(golden)

	#m2, white as it will be cutting out this region from the band 
	dchi2_nh_m2sigma.SetFillColor(0)
	dchi2_ih_m2sigma.SetFillColor(0)
	dchi2_nh_m2sigma.SetLineColor(golden)
	dchi2_ih_m2sigma.SetLineColor(golden)

	median_nh.SetLineStyle(4)
	median_nh.SetLineColor(1)
	median_nh.SetLineWidth(2)

	median_ih.SetLineStyle(4)
	median_ih.SetLineColor(1)
	median_ih.SetLineWidth(2)

 	hdata_nh.SetLineColor(1)
 	hdata_nh.SetLineWidth(2)
 	hdata_ih.SetLineStyle(1)
 	hdata_ih.SetLineColor(1)
 	hdata_ih.SetLineWidth(2)

 	#Set up a legend 
 	lBrazilianPlot_nh = ROOT.TLegend(0.17,0.64,0.37,0.92)
 	lBrazilianPlot_nh.SetFillColor(0)
 	lBrazilianPlot_nh.SetLineColor(0)
 	lBrazilianPlot_nh.SetTextFont(42)
 	lBrazilianPlot_nh.SetTextSize(0.045)
	lBrazilianPlot_nh.AddEntry(dchi2_nh_p1sigma,"68.27% of  pseudo-experiments","f")
	lBrazilianPlot_nh.AddEntry(dchi2_nh_p2sigma,"95.45% of  pseudo-experiments","f")
	lBrazilianPlot_nh.AddEntry(median_nh,"Median pseudo-experiment","l")
	lBrazilianPlot_nh.AddEntry(hdata_nh,"Data","l")

	#Set up a legend 
 	lBrazilianPlot_ih = ROOT.TLegend(0.17,0.64,0.37,0.92)
 	lBrazilianPlot_ih.SetFillColor(0)
 	lBrazilianPlot_ih.SetLineColor(0)
 	lBrazilianPlot_ih.SetTextFont(42)
 	lBrazilianPlot_ih.SetTextSize(0.045)
	lBrazilianPlot_ih.AddEntry(dchi2_ih_p1sigma,"68.27% of  pseudo-experiments","f")
	lBrazilianPlot_ih.AddEntry(dchi2_ih_p2sigma,"95.45% of  pseudo-experiments","f")
	lBrazilianPlot_ih.AddEntry(median_ih,"Median pseudo-experiment","l")
	lBrazilianPlot_ih.AddEntry(hdata_ih,"Data","l")

	#drawing and printing NH
	c_nh.cd()
	dchi2_nh_p2sigma.Draw("c")
	dchi2_nh_p1sigma.Draw("c same")
	dchi2_nh_m1sigma.Draw("c same")
	dchi2_nh_m2sigma.Draw("c same")

	median_nh.Draw("l same")
	hdata_nh.Draw("l same")
	lBrazilianPlot_nh.Draw("same")
	ROOT.gPad.RedrawAxis()
	nname = out_filename+"_two_sided_NH.pdf"
	c_nh.Print(nname)
	
	#drawing and printing IH
	c_ih.cd()
	dchi2_ih_p2sigma.Draw("c")
	dchi2_ih_p1sigma.Draw("c same")
	dchi2_ih_m1sigma.Draw("c same")
	dchi2_ih_m2sigma.Draw("c same")

	median_ih.Draw("l same")
	hdata_ih.Draw("l same")
	lBrazilianPlot_ih.Draw("same")
	ROOT.gPad.RedrawAxis()
	iname = out_filename+"_two_sided_IH.pdf"
	c_ih.Print(iname)

	#writting to a file 
	dchi2_ih_p1sigma.Write()
	dchi2_ih_p2sigma.Write()
	dchi2_nh_p1sigma.Write()
	dchi2_nh_p2sigma.Write()
	dchi2_ih_m1sigma.Write()
	dchi2_ih_m2sigma.Write()
	dchi2_nh_m1sigma.Write()
	dchi2_nh_m2sigma.Write()
	c_nh.Write()
	c_ih.Write()	

	fout.Write()
	fout.Close() 		

#_______________________________________________________________________
def combine(in_filenames):
	#Lists of lists to store dchi valies per marg point for each toy 	
	_nh_dchi_list = []
	_ih_dchi_list = []

	print "Combining dchi2 values for each marg point form all batch of toys. Can take some time. Be patient"
	file_idx = 0 #need this for brazilian arrays reading, to copy first array and append all next 
	for filename in in_filenames:
		print "Reading from file", filename
		in_file = ROOT.TFile(filename, "read")  
		#addining together info for each batch of 1000 toys for each dcp bin 
		br_nh_arr_tmp = in_file.Get(_dchi2_arr_nh_name)
		br_ih_arr_tmp = in_file.Get(_dchi2_arr_ih_name)
		if file_idx==0: # clone first array 
			_nh_dchi_list = as_list(br_nh_arr_tmp)
			_ih_dchi_list = as_list(br_ih_arr_tmp)
		else: 
			for ibin in range (0, br_nh_arr_tmp.GetEntries()):
				for elem in range (0, br_nh_arr_tmp[ibin].GetNoElements()):	
					_nh_dchi_list[ibin].append(br_nh_arr_tmp[ibin](elem))
			for jbin in range (0, br_ih_arr_tmp.GetEntries()):
				for elem in range (0, br_ih_arr_tmp[ibin].GetNoElements()):	
					_ih_dchi_list[jbin].append(br_ih_arr_tmp[jbin](elem))    
		in_file.Close()
		del in_file
		file_idx+=1 # incrementing counter after reading the first file 
	return [_nh_dchi_list,_ih_dchi_list]
#_______________________________________________________________________
def as_list(in_array): #converts an TObjArray of TVectorTs in to list of lists, as it is way better to work with lists in python 
	dchi_mpt_list = []

	if in_array.GetEntries() ==0:
		print "Empty array! Time to die"
		sys.exit()
			
	for ibin in range (0, in_array.GetEntries()):
		dchi_toy_list_tmp = []
		
		if in_array[ibin].GetNoElements()==0:
			print "Empty vector! Time to die"
			sys.exit()

		for elem in range (0, in_array[ibin].GetNoElements()):
			dchi_toy_list_tmp.append(in_array[ibin](elem))
		dchi_mpt_list.append(dchi_toy_list_tmp)
		del dchi_toy_list_tmp
	return dchi_mpt_list

#_______________________________________________________________________
def calc_critical(dchi2_vec, conf_level_enum):
	conf_level = util.ConfLevelToFrac(conf_level_enum)
	if len(dchi2_vec)==0:
		print"Empty vector! What's the point in sorting it? "
	dchi2_vec.sort()
	ntoys_conf_level = conf_level * len(dchi2_vec)  # The number of toys contained within this confidence level
	index = int(ntoys_conf_level)
	if index >= len(dchi2_vec):
	 print "ERROR: the index is outside of the Delta chisq list range.\nconf_level =", \
	 conf_level, ", index =", index, ", length =", len(dchi2_vec), \
	 "\n\nDelta chisq values =", dchi2_vec
	 sys.exit()
	return dchi2_vec[index]

#_______________________________________________________________________
def calc_critical_binned(dchisq_list, conf_level_enum):

	conf_level = util.ConfLevelToFrac(conf_level_enum)
	dchisq_list.sort()
	ntoys = len(dchisq_list)
		
	# Find median dchisq value and its index
	median_index = None
	median       = None
	if ntoys % 2 != 0:
		# Odd number of toys
		median_index = int(ntoys/2)	# This should round down, as we need for this (because we're counting the list index from 0)
		median = dchisq_list[median_index]
	else:
		# Even number of toys
		median_index = int(ntoys/2)	# Even if the median dchisq value doesn't correspond exactly to this index (as it's the average of two values), this is close enough given the large number of toys used
		median = 0.5 * (dchisq_list[median_index-1] + dchisq_list[median_index])
	
	# Find range of dchisq values around the 
	# median that encompases ntoys * conf_level toys
	range_index = int(math.ceil(0.5 * float(ntoys) * conf_level))
	dchisq_low_index = median_index - range_index
	dchisq_up_index  = median_index + range_index
	if dchisq_low_index < 0:
		dchisq_low_index = 0
	if dchisq_up_index >= ntoys:
		dchisq_up_index = ntoys-1
	
	dchisq_range = [dchisq_list[dchisq_low_index], dchisq_list[dchisq_up_index]]
	return dchisq_range, median
#_______________________________________________
if __name__=="__main__":
 	import argparse

 	parser = argparse.ArgumentParser(description="Generate Brazilian plot")
 	parser.add_argument("-f", "--filenames"			, dest="in_filenames", nargs="+"													, help="Input ROOT filenames containing chi2 arrays.")
 	parser.add_argument("-o", "--output"   			, dest="out_filename"           													, help="Output ROOT filename")
 	parser.add_argument("-d", "--data-fit-nh"   , dest="data_fit_nh"    																	, help="File containing dchi2 data fit for NH")
 	parser.add_argument("-i", "--data-fit-ih"   , dest="data_fit_ih"    																	, help="File containing dchi2 data fit for IH")
	parser.add_argument('--two-sided'						, dest='two_sided'   , default=False, action="store_true"	, help="Set this flag if you want to plot two side band ")
	args = parser.parse_args()
	if args.two_sided:
		draw_two_sided(args.in_filenames, args.out_filename, args.data_fit_nh, args.data_fit_ih)
	else: 
		draw_one_sided(args.in_filenames, args.out_filename, args.data_fit_nh, args.data_fit_ih)