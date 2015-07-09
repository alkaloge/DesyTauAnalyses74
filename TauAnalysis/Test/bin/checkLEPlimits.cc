
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "TauAnalysis/Test/interface/mssm_xs_tools.h"

#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TMath.h>

#include <iostream>
#include <string>
#include <vector>

double getLimit_tanBeta_boundByBelow(mssm_xs_tools& mssmTool, double mA, double hMassMin, double hMassMax)
{
  double limit_tanBeta = 1.e+3;
  for ( double tanBeta = 60.; tanBeta >= 1.; tanBeta -= 0.1 ) {
    double tanBeta_floor = TMath::Floor(tanBeta);
    double hMass_floor = mssmTool.Give_Mass_h(mA, tanBeta_floor);
    double tanBeta_ceil = TMath::Ceil(tanBeta);
    if ( !(tanBeta_ceil - tanBeta_floor) > 0.5 ) tanBeta_ceil += 1.;
    double hMass_ceil = mssmTool.Give_Mass_h(mA, tanBeta_ceil);
    double hMass = (tanBeta - tanBeta_floor)*hMass_ceil + (tanBeta_ceil - tanBeta)*hMass_floor;
    //std::cout << "tanBeta = " << tanBeta << ": tanBeta_floor = " << tanBeta_floor << ", tanBeta_ceil = " << tanBeta_ceil << std::endl;
    //std::cout << " hMass_floor = " << hMass_floor << ", hMass_ceil = " << hMass_ceil << " --> hMass = " << hMass << std::endl; 
    if ( hMass > hMassMin && hMass < hMassMax ) {
      limit_tanBeta = tanBeta;
    } else {
      std::cout << "mA = " << mA << ": tanBeta = " << limit_tanBeta << " (mh = " << hMass << ")" << std::endl;
      break;
    }
  }
  return limit_tanBeta;
}

double getLimit_tanBeta_boundByAbove(mssm_xs_tools& mssmTool, double mA, double hMassMin, double hMassMax)
{
  double limit_tanBeta = 0.;
  for ( double tanBeta = 0.5; tanBeta <= 0.9; tanBeta += 0.01 ) {
    double tanBeta_floor = 0.1*TMath::Floor(10.*tanBeta);
    double hMass_floor = mssmTool.Give_Mass_h(mA, tanBeta_floor);
    double tanBeta_ceil = 0.1*TMath::Ceil(10.*tanBeta);
    if ( !(tanBeta_ceil - tanBeta_floor) > 0.05 ) tanBeta_ceil += 0.1;
    double hMass_ceil = mssmTool.Give_Mass_h(mA, tanBeta_ceil);
    double hMass = 10.*(tanBeta - tanBeta_floor)*hMass_ceil + 10.*(tanBeta_ceil - tanBeta)*hMass_floor;
    std::cout << "tanBeta = " << tanBeta << ": tanBeta_floor = " << tanBeta_floor << ", tanBeta_ceil = " << tanBeta_ceil << std::endl;
    std::cout << " hMass_floor = " << hMass_floor << ", hMass_ceil = " << hMass_ceil << " --> hMass = " << hMass << std::endl; 
    if ( hMass > hMassMin && hMass < hMassMax ) {
      limit_tanBeta = tanBeta;
    } else {
      std::cout << "mA = " << mA << ": tanBeta = " << limit_tanBeta << " (mh = " << hMass << ")" << std::endl;
      break;
    }
  }
  return limit_tanBeta;
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		TGraph* graph5, const std::string& legendEntry5,
		TGraph* graph6, const std::string& legendEntry6,
		int colors[], int markerStyles[], 
		double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
		std::vector<std::string>& labelTextLines, double labelTextSize,
		double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", 100, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw("axis");

  graph1->SetLineColor(colors[0]);
  graph1->SetLineWidth(2);
  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->SetMarkerSize(1);
  graph1->Draw("L");

  if ( graph2 ) {
    graph2->SetLineColor(colors[1]);
    graph2->SetLineWidth(2);
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetMarkerSize(1);
    graph2->Draw("L");
  }
  
  if ( graph3 ) {
    graph3->SetLineColor(colors[2]);
    graph3->SetLineWidth(2);
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetMarkerSize(1);
    graph3->Draw("L");
  }

  if ( graph4 ) {
    graph4->SetLineColor(colors[3]);
    graph4->SetLineWidth(2);
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetMarkerSize(1);
    graph4->Draw("L");
  }

  if ( graph5 ) {
    graph5->SetLineColor(colors[4]);
    graph5->SetLineWidth(2);
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->SetMarkerSize(1);
    graph5->Draw("L");
  }

  if ( graph6 ) {
    graph6->SetLineColor(colors[5]);
    graph6->SetLineWidth(2);
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->SetMarkerSize(1);
    graph6->Draw("L");
  }
  
  TLegend* legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(legendTextSize);
  legend->AddEntry(graph1, legendEntry1.data(), "p");
  if ( graph2 ) legend->AddEntry(graph2, legendEntry2.data(), "p");
  if ( graph3 ) legend->AddEntry(graph3, legendEntry3.data(), "p");
  if ( graph4 ) legend->AddEntry(graph4, legendEntry4.data(), "p");
  if ( graph5 ) legend->AddEntry(graph5, legendEntry5.data(), "p");
  if ( graph6 ) legend->AddEntry(graph6, legendEntry6.data(), "p");
  legend->Draw();

  TPaveText* label = 0;
  if ( labelTextLines.size() > 0 ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "brNDC");
    for ( std::vector<std::string>::const_iterator labelTextLine = labelTextLines.begin();
	  labelTextLine != labelTextLines.end(); ++labelTextLine ) {
      label->AddText(labelTextLine->data());
    }
    label->SetFillColor(10);
    label->SetBorderSize(0);
    label->SetTextColor(1);
    label->SetTextAlign(12);
    label->SetTextSize(labelTextSize);
    label->Draw();
  }

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete dummyHistogram;
  delete label;
  delete legend;
  delete canvas;  
}

int main(int argc, char* argv[]) 
{
  std::cout << "<checkLEPlimits>:" << std::endl;

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  std::map<double, double> LEPlimit;
  LEPlimit[50.]      = 100.;
  LEPlimit[91.]      = 100; 
  LEPlimit[91.8]     = 30.02624; 
  LEPlimit[91.845]   = 22.07032; 
  LEPlimit[91.845]   = 17.12491; 
  LEPlimit[91.84523] = 13.64727; 
  LEPlimit[92.61388] = 11.94143; 
  LEPlimit[93.38253] = 10.03852; 
  LEPlimit[94.91982] = 9.021481; 
  LEPlimit[95.68846] = 8.107481; 
  LEPlimit[97.22578] = 7.141608; 
  LEPlimit[99.5317]  = 6.680381; 
  LEPlimit[103.375]  = 7.189448; 
  LEPlimit[104.1436] = 7.841313; 
  LEPlimit[106.4496] = 8.326916; 
  LEPlimit[109.5242] = 8.609568; 
  LEPlimit[112.5988] = 8.438845; 
  LEPlimit[115.6733] = 8.107481; 
  LEPlimit[118.748]  = 7.384029; 
  LEPlimit[122.5912] = 6.547911; 
  LEPlimit[126.4344] = 5.963618; 
  LEPlimit[131.815]  = 5.359424; 
  LEPlimit[138.7328] = 4.752558; 
  LEPlimit[144.1134] = 4.445624; 
  LEPlimit[149.4939] = 4.186368; 
  LEPlimit[156.4118] = 3.968637; 
  LEPlimit[164.8669] = 3.687628; 
  LEPlimit[177.1653] = 3.472575; 
  LEPlimit[187.9264] = 3.29197; 
  LEPlimit[203.2994] = 3.141663; 
  LEPlimit[221.7469] = 2.978266; 
  LEPlimit[241.7318] = 2.861322; 
  LEPlimit[261.7167] = 2.767383; 
  LEPlimit[283.2388] = 2.676528; 
  LEPlimit[304.761]  = 2.641027; 
  LEPlimit[334.7383] = 2.554322; 
  LEPlimit[357.0292] = 2.50367; 
  LEPlimit[383.9319] = 2.48701; 
  LEPlimit[420.8271] = 2.454023; 
  LEPlimit[452.3417] = 2.421473; 
  LEPlimit[487.6996] = 2.405361; 
  LEPlimit[500.]     = 2.405361; 

  mssm_xs_tools mssmTool_mhmax_boundByBelow;
  mssmTool_mhmax_boundByBelow.SetInput("../data_nocrab/out.mhmax-mu200-8TeV-tanbHigh-nnlo.root.1");
  mssm_xs_tools mssmTool_mhmax_boundByAbove;
  mssmTool_mhmax_boundByAbove.SetInput("../data_nocrab/out.mhmax-mu200-8TeV-tanbLow-nnlo.root");

  mssm_xs_tools mssmTool_mhmodPlus_boundByBelow;
  mssmTool_mhmodPlus_boundByBelow.SetInput("../data_nocrab/mhmodp_8TeV_tanbHigh.root");
  mssm_xs_tools mssmTool_mhmodPlus_boundByAbove;
  mssmTool_mhmodPlus_boundByAbove.SetInput("../data_nocrab/mhmodp_8TeV_tanbLow.root");

  mssm_xs_tools mssmTool_mhmodMinus_boundByBelow;
  mssmTool_mhmodMinus_boundByBelow.SetInput("../data_nocrab/mhmodm_8TeV_tanbHigh.root");
  mssm_xs_tools mssmTool_mhmodMinus_boundByAbove;
  mssmTool_mhmodMinus_boundByAbove.SetInput("../data_nocrab/mhmodm_8TeV_tanbLow.root");

  int numPoints = LEPlimit.size();

  TGraph* graph_LEPlimit = new TGraph(numPoints);
  
  int idxPoint = 0;
  for ( std::map<double, double>::const_iterator LEPlimit_i = LEPlimit.begin();
	LEPlimit_i != LEPlimit.end(); ++LEPlimit_i ) {
    double mA = LEPlimit_i->first;

    double LEPlimit_tanBeta = LEPlimit_i->second;
    graph_LEPlimit->SetPoint(idxPoint, mA, LEPlimit_tanBeta);

    ++idxPoint;
  }

  TGraph* graph_mhmaxLimit_boundByBelow = new TGraph(401);
  TGraph* graph_mhmaxLimit_boundByAbove = new TGraph(401);
  TGraph* graph_mhmodPlusLimit_boundByBelow = new TGraph(401);
  TGraph* graph_mhmodPlusLimit_boundByAbove = new TGraph(401);
  TGraph* graph_mhmodMinusLimit_boundByBelow = new TGraph(401);
  TGraph* graph_mhmodMinusLimit_boundByAbove = new TGraph(401);
      
  //const double mhlimit = 114.;
  const double mhlimit = 122.;

  idxPoint = 0;
  for ( double mA = 90.; mA <= 500.; mA += 1. ) {
    double mhmaxLimit_tanBeta_boundByBelow = getLimit_tanBeta_boundByBelow(mssmTool_mhmax_boundByBelow, mA, mhlimit, 1.e+3);
    graph_mhmaxLimit_boundByBelow->SetPoint(idxPoint, mA, mhmaxLimit_tanBeta_boundByBelow);
    double mhmaxLimit_tanBeta_boundByAbove = getLimit_tanBeta_boundByAbove(mssmTool_mhmax_boundByAbove, mA, mhlimit, 1.e+3);
    graph_mhmaxLimit_boundByAbove->SetPoint(idxPoint, mA, mhmaxLimit_tanBeta_boundByAbove);

    double mhmodPlusLimit_tanBeta_boundByBelow = getLimit_tanBeta_boundByBelow(mssmTool_mhmodPlus_boundByBelow, mA, mhlimit, 1.e+3);
    graph_mhmodPlusLimit_boundByBelow->SetPoint(idxPoint, mA, mhmodPlusLimit_tanBeta_boundByBelow);
    double mhmodPlusLimit_tanBeta_boundByAbove = getLimit_tanBeta_boundByAbove(mssmTool_mhmodPlus_boundByAbove, mA, mhlimit, 1.e+3);
    graph_mhmodPlusLimit_boundByAbove->SetPoint(idxPoint, mA, mhmodPlusLimit_tanBeta_boundByAbove);
    
    double mhmodMinusLimit_tanBeta_boundByBelow = getLimit_tanBeta_boundByBelow(mssmTool_mhmodMinus_boundByBelow, mA, mhlimit, 1.e+3);
    graph_mhmodMinusLimit_boundByBelow->SetPoint(idxPoint, mA, mhmodMinusLimit_tanBeta_boundByBelow);
    double mhmodMinusLimit_tanBeta_boundByAbove = getLimit_tanBeta_boundByAbove(mssmTool_mhmodMinus_boundByAbove, mA, mhlimit, 1.e+3);
    graph_mhmodMinusLimit_boundByAbove->SetPoint(idxPoint, mA, mhmodMinusLimit_tanBeta_boundByAbove);

    ++idxPoint;
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLogy(true);

  TH1* dummyHistogram = new TH1F("dummyHistogram", "dummyHistogram", 41, 90., 500.);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(0.5);
  dummyHistogram->SetMaximum(60.);
  dummyHistogram->GetXaxis()->SetTitle("m_{A} / GeV");
  dummyHistogram->GetYaxis()->SetTitle("tan#beta");

  dummyHistogram->Draw("axis");

  graph_LEPlimit->SetLineColor(8);
  graph_LEPlimit->SetLineWidth(-802);
  graph_LEPlimit->SetFillStyle(3004);
  graph_LEPlimit->SetFillColor(8);
  graph_LEPlimit->Draw("L");

  graph_mhmaxLimit_boundByBelow->SetLineColor(2);
  graph_mhmaxLimit_boundByBelow->SetLineWidth(-402);
  graph_mhmaxLimit_boundByBelow->SetFillStyle(3005);
  graph_mhmaxLimit_boundByBelow->SetFillColor(2);
  graph_mhmaxLimit_boundByBelow->Draw("L");

  if ( mhlimit < 120. ) {
    graph_mhmaxLimit_boundByAbove->SetLineColor(2);
    graph_mhmaxLimit_boundByAbove->SetLineWidth(402);
    graph_mhmaxLimit_boundByAbove->SetFillStyle(3005);
    graph_mhmaxLimit_boundByAbove->SetFillColor(2);
    graph_mhmaxLimit_boundByAbove->Draw("L");
  }

  graph_mhmodPlusLimit_boundByBelow->SetLineColor(4);
  graph_mhmodPlusLimit_boundByBelow->SetLineWidth(-402);
  graph_mhmodPlusLimit_boundByBelow->SetFillStyle(3004);
  graph_mhmodPlusLimit_boundByBelow->SetFillColor(4);
  graph_mhmodPlusLimit_boundByBelow->Draw("L");

/*
  // CV: buggy, disable plot for now
  graph_mhmodPlusLimit_boundByAbove->SetLineColor(4);
  graph_mhmodPlusLimit_boundByAbove->SetLineWidth(402);
  graph_mhmodPlusLimit_boundByAbove->SetFillStyle(3004);
  graph_mhmodPlusLimit_boundByAbove->SetFillColor(4);
  graph_mhmodPlusLimit_boundByAbove->Draw("L");
 */

  graph_mhmodMinusLimit_boundByBelow->SetLineColor(7);
  graph_mhmodMinusLimit_boundByBelow->SetLineWidth(-402);
  graph_mhmodMinusLimit_boundByBelow->SetFillStyle(3005);
  graph_mhmodMinusLimit_boundByBelow->SetFillColor(7);
  graph_mhmodMinusLimit_boundByBelow->Draw("L");

/*
  // CV: buggy, disable plot for now
    graph_mhmodMinusLimit_boundByAbove->SetLineColor(7);
    graph_mhmodMinusLimit_boundByAbove->SetLineWidth(402);
    graph_mhmodMinusLimit_boundByAbove->SetFillStyle(3005);
    graph_mhmodMinusLimit_boundByAbove->SetFillColor(7);
    graph_mhmodMinusLimit_boundByAbove->Draw("L"); 
 */

  TLegend* legend = new TLegend(0.67, 0.65, 0.87, 0.89, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(graph_LEPlimit, "LEP", "l");
  legend->AddEntry(graph_mhmaxLimit_boundByBelow, "m_{h}^{max}", "l");
  legend->AddEntry(graph_mhmodPlusLimit_boundByBelow, "m_{h}^{mod+}", "l");
  legend->AddEntry(graph_mhmodMinusLimit_boundByBelow, "m_{h}^{mod-}", "l");
  legend->Draw();

  canvas->Print("plots/checkLEPlimits.pdf");

  delete dummyHistogram;

  std::vector<double> tanBeta_values;
  tanBeta_values.push_back(1.);
  tanBeta_values.push_back(3.);
  tanBeta_values.push_back(5.);
  tanBeta_values.push_back(7.);
  tanBeta_values.push_back(10.);
  tanBeta_values.push_back(15.);
  std::map<double, TGraph*> tanBeta_graphs_mhmax;
  std::map<double, TGraph*> tanBeta_graphs_mhmodPlus;
  std::map<double, TGraph*> tanBeta_graphs_mhmodMinus;
  for ( std::vector<double>::const_iterator tanBeta = tanBeta_values.begin();
	tanBeta != tanBeta_values.end(); ++tanBeta ) {
    TGraph* graph_mhmax = new TGraph(411);
    TGraph* graph_mhmodPlus = new TGraph(411);
    TGraph* graph_mhmodMinus = new TGraph(411);
    int idxPoint = 0;
    for ( double mA = 90.; mA <= 500.; mA += 1. ) {
      double mh_mhmax = mssmTool_mhmax_boundByBelow.Give_Mass_h(mA, *tanBeta);
      //std::cout << "mhmax: mA = " << mA << ", tanBeta = " << (*tanBeta) << " -> mh = " << mh_mhmax << std::endl;
      graph_mhmax->SetPoint(idxPoint, mA, mh_mhmax);
      double mh_mhmodPlus = mssmTool_mhmodPlus_boundByBelow.Give_Mass_h(mA, *tanBeta);
      //std::cout << "mhmodPlus: mA = " << mA << ", tanBeta = " << (*tanBeta) << " -> mh = " << mh_mhmodPlus << std::endl;
      graph_mhmodPlus->SetPoint(idxPoint, mA, mh_mhmodPlus);
      double mh_mhmodMinus = mssmTool_mhmodMinus_boundByBelow.Give_Mass_h(mA, *tanBeta);
      //std::cout << "mhmodMinus: mA = " << mA << ", tanBeta = " << (*tanBeta) << " -> mh = " << mh_mhmodMinus << std::endl;
      graph_mhmodMinus->SetPoint(idxPoint, mA, mh_mhmodMinus);
      ++idxPoint;
    }
    tanBeta_graphs_mhmax[*tanBeta] = graph_mhmax;
    tanBeta_graphs_mhmodPlus[*tanBeta] = graph_mhmodPlus;
    tanBeta_graphs_mhmodMinus[*tanBeta] = graph_mhmodMinus;
  }

  int colors[6] = { 1, 2, 3, 4, 6, 8 };  
  //int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };
  int markerStyles[6] = { 8, 8, 8, 8, 8, 8 };

  std::vector<std::string> labelText_mhmax;
  labelText_mhmax.push_back("m_{h}^{max}");
  std::vector<std::string> labelText_mhmodPlus;
  labelText_mhmodPlus.push_back("m_{h}^{mod+}");
  std::vector<std::string> labelText_mhmodMinus;
  labelText_mhmodMinus.push_back("m_{h}^{mod-}");

  showGraphs(800, 600,
	     tanBeta_graphs_mhmax[1.],  "tan#beta = 1",
	     tanBeta_graphs_mhmax[3.],  "tan#beta = 3",
	     tanBeta_graphs_mhmax[5.],  "tan#beta = 5",
	     tanBeta_graphs_mhmax[7.],  "tan#beta = 7",
	     tanBeta_graphs_mhmax[10.], "tan#beta = 10",
	     tanBeta_graphs_mhmax[15.], "tan#beta = 15",
	     colors, markerStyles,
	     0.04, 0.61, 0.65, 0.28, 0.24,
	     labelText_mhmax, 0.06, 0.175, 0.725, 0.24, 0.165, 
	     90., 500., "m_{A} / GeV", 1.2,
	     70., 160., "m_{h} / GeV", 1.4,
	     "checkLEPlimits_mh_vs_mA_mhmax.pdf");

  showGraphs(800, 600,
	     tanBeta_graphs_mhmodPlus[1.],  "tan#beta = 1",
	     tanBeta_graphs_mhmodPlus[3.],  "tan#beta = 3",
	     tanBeta_graphs_mhmodPlus[5.],  "tan#beta = 5",
	     tanBeta_graphs_mhmodPlus[7.],  "tan#beta = 7",
	     tanBeta_graphs_mhmodPlus[10.], "tan#beta = 10",
	     tanBeta_graphs_mhmodPlus[15.], "tan#beta = 15",
	     colors, markerStyles,
	     0.04, 0.61, 0.65, 0.28, 0.24,
	     labelText_mhmodPlus, 0.06, 0.175, 0.725, 0.24, 0.165, 
	     90., 500., "m_{A} / GeV", 1.2,
	     70., 160., "m_{h} / GeV", 1.4,
	     "checkLEPlimits_mh_vs_mA_mhmodPlus.pdf");

  showGraphs(800, 600,
	     tanBeta_graphs_mhmodMinus[1.],  "tan#beta = 1",
	     tanBeta_graphs_mhmodMinus[3.],  "tan#beta = 3",
	     tanBeta_graphs_mhmodMinus[5.],  "tan#beta = 5",
	     tanBeta_graphs_mhmodMinus[7.],  "tan#beta = 7",
	     tanBeta_graphs_mhmodMinus[10.], "tan#beta = 10",
	     tanBeta_graphs_mhmodMinus[15.], "tan#beta = 15",
	     colors, markerStyles,
	     0.04, 0.61, 0.65, 0.28, 0.24,
	     labelText_mhmodMinus, 0.06, 0.175, 0.725, 0.24, 0.165, 
	     90., 500., "m_{A} / GeV", 1.2,
	     70., 160., "m_{h} / GeV", 1.4,
	     "checkLEPlimits_mh_vs_mA_mhmodMinus.pdf");

  std::vector<double> mA_values;
  mA_values.push_back(138.7328);
  mA_values.push_back(144.1134);
  mA_values.push_back(149.4939);
  mA_values.push_back(156.4118);
  mA_values.push_back(164.8669);
  mA_values.push_back(177.1653);
  std::map<double, TGraph*> mA_graphs_mhmax;
  std::map<double, TGraph*> mA_graphs_mhmodPlus;
  std::map<double, TGraph*> mA_graphs_mhmodMinus;
  for ( std::vector<double>::const_iterator mA = mA_values.begin();
	mA != mA_values.end(); ++mA ) {
    TGraph* graph_mhmax = new TGraph(1401);
    TGraph* graph_mhmodPlus = new TGraph(1401);
    TGraph* graph_mhmodMinus = new TGraph(1401);
    int idxPoint = 0;
    for ( double tanBeta = 1.; tanBeta <= 15.; tanBeta += 1.e-2 ) {
      double mh_mhmax = mssmTool_mhmax_boundByBelow.Give_Mass_h(*mA, tanBeta);
      graph_mhmax->SetPoint(idxPoint, tanBeta, mh_mhmax);
      double mh_mhmodPlus = mssmTool_mhmodPlus_boundByBelow.Give_Mass_h(*mA, tanBeta);
      graph_mhmodPlus->SetPoint(idxPoint, tanBeta, mh_mhmodPlus);      
      double mh_mhmodMinus = mssmTool_mhmodMinus_boundByBelow.Give_Mass_h(*mA, tanBeta);
      graph_mhmodMinus->SetPoint(idxPoint, tanBeta, mh_mhmodMinus);
      ++idxPoint;
    }
    mA_graphs_mhmax[*mA] = graph_mhmax;
    mA_graphs_mhmodPlus[*mA] = graph_mhmodPlus;
    mA_graphs_mhmodMinus[*mA] = graph_mhmodMinus;
  }

  showGraphs(800, 600,
	     mA_graphs_mhmax[138.7328], "m_{A} = 138.7 GeV",
	     mA_graphs_mhmax[144.1134], "m_{A} = 144.1 GeV",
	     mA_graphs_mhmax[149.4939], "m_{A} = 149.5 GeV",
	     mA_graphs_mhmax[156.4118], "m_{A} = 156.4 GeV",
	     mA_graphs_mhmax[164.8669], "m_{A} = 164.9 GeV",
	     mA_graphs_mhmax[177.1653], "m_{A} = 177.2 GeV",
	     colors, markerStyles,
	     0.04, 0.61, 0.65, 0.28, 0.24,
	     labelText_mhmax, 0.06, 0.175, 0.725, 0.24, 0.165, 
	     1., 15., "tan#beta", 1.2,
	     70., 160., "m_{h} / GeV", 1.4,
	     "checkLEPlimits_mh_vs_tanBeta_mhmax.pdf");

  showGraphs(800, 600,
	     mA_graphs_mhmodPlus[138.7328], "m_{A} = 138.7 GeV",
	     mA_graphs_mhmodPlus[144.1134], "m_{A} = 144.1 GeV",
	     mA_graphs_mhmodPlus[149.4939], "m_{A} = 149.5 GeV",
	     mA_graphs_mhmodPlus[156.4118], "m_{A} = 156.4 GeV",
	     mA_graphs_mhmodPlus[164.8669], "m_{A} = 164.9 GeV",
	     mA_graphs_mhmodPlus[177.1653], "m_{A} = 177.2 GeV",
	     colors, markerStyles,
	     0.04, 0.61, 0.65, 0.28, 0.24,
	     labelText_mhmodPlus, 0.06, 0.175, 0.725, 0.24, 0.165, 
	     1., 15., "tan#beta", 1.2,
	     70., 160., "m_{h} / GeV", 1.4,
	     "checkLEPlimits_mh_vs_tanBeta_mhmodPlus.pdf");

  showGraphs(800, 600,
	     mA_graphs_mhmodMinus[138.7328], "m_{A} = 138.7 GeV",
	     mA_graphs_mhmodMinus[144.1134], "m_{A} = 144.1 GeV",
	     mA_graphs_mhmodMinus[149.4939], "m_{A} = 149.5 GeV",
	     mA_graphs_mhmodMinus[156.4118], "m_{A} = 156.4 GeV",
	     mA_graphs_mhmodMinus[164.8669], "m_{A} = 164.9 GeV",
	     mA_graphs_mhmodMinus[177.1653], "m_{A} = 177.2 GeV",
	     colors, markerStyles,
	     0.04, 0.61, 0.65, 0.28, 0.24,
	     labelText_mhmodMinus, 0.06, 0.175, 0.725, 0.24, 0.165, 
	     1., 15., "tan#beta", 1.2,
	     70., 160., "m_{h} / GeV", 1.4,
	     "checkLEPlimits_mh_vs_tanBeta_mhmodMinus.pdf");
  
  return 0;
}
