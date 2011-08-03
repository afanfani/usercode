//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooVoigtian.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include <TROOT.h>
#include "TH1F.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include <string>
#include "/home/CMS/fanfani/tdrstyle.C"
using namespace RooFit ;

void FitMass_FromRooDataSet(TString name, TString cutname, float XMinPt, float XMinY, float piMinPt, float dimuVtx, float piSiHits, float JpsiMinY, float JpsiMaxY, TString PsiPrimePDF="Voigt", int framebins=100, int Log=0, TString nprompt=""){


// usage:
// .L FitMass_FromRooDataSet.C++
// FitMass_FromRooDataSet("RooNew_May10Promptv41Jul_JPsiRap_All","test",9.,1.25,0.6,0.01,7.,0.,1.25,"Voigt",100);

setTDRStyle();
gROOT->ForceStyle();

gROOT->Reset();
gROOT->GetList()->Delete();
gROOT->GetListOfCanvases()->Delete();

TFile *file = new TFile("RooDataset_"+name+".root");
RooDataSet* FulldataSet = (RooDataSet*)file->Get("Mydata");


//// Important: keep the same name as the RooRealVar used to fill the RooDataset
 float RangeMassMin = 3.6 , RangeMassMax =4;
 RooRealVar X_Mass("X_Mass", "m(J/#Psi + #pi^{+}#pi^{-}) [GeV]",  RangeMassMin,  RangeMassMax) ;
 RooRealVar HardPion_pt = RooRealVar("HardPion_pt","HardPion Pt",0.4,200.);
 RooRealVar SoftPion_pt = RooRealVar("SoftPion_pt","SoftPion Pt",0.4,200.);
 RooRealVar HardPion_Hits = RooRealVar("HardPion_Hits","HardPion Hits",0.4,200.);
 RooRealVar SoftPion_Hits = RooRealVar("SoftPion_Hits","SoftPion Hits",0.4,200.);
 RooRealVar X_Pt = RooRealVar("X_Pt","X(3872) Pt",0.,200.);
 RooRealVar X_Rap = RooRealVar("X_Rap","X(3872) Rap", -2.4,2.4);
 RooRealVar X_ct = RooRealVar("X_ct","X(3872) ctau",-2.,5.);
 RooRealVar X_pdl = RooRealVar("X_pdl","X(3872) pdl",-2.,2.);
 RooRealVar JPsi_Rap = RooRealVar("JPsi_Rap","J/Psi Rap",-2.4,2.4);
 RooCategory DiMuonType = RooCategory("DiMuonType","Type of Dimuon");
 DiMuonType.defineType("A",1);
 DiMuonType.defineType("B",-1);
 RooRealVar X_Vtx = RooRealVar("X_Vtx","X(3872) Vtx Probability",0.01,1.,"");
 RooRealVar X_QValue = RooRealVar("X_QValue","X(3872) Q Value",-3.,3.);
 RooRealVar X_QValuePDG = RooRealVar("X_QValuePDG","X(3872) Q Value",-3.,3.);
//// Reduce dataset only to Mass 
 const RooArgSet* rowFull = FulldataSet->get() ;
 RooRealVar* xmassFull = (RooRealVar*) rowFull->find("X_Mass");


  char cutstring[300];
  if (nprompt=="pdl"){
     sprintf(cutstring,"X_pdl>0.01 && X_Pt>%f && X_Rap>%f && X_Rap<%f && SoftPion_pt>%f && HardPion_pt>%f && DiMuonVtx>%f && SoftPion_Hits>=%f && HardPion_Hits>=%f && abs(JPsi_Rap)>=%f && abs(JPsi_Rap)<%f",XMinPt,-XMinY,XMinY,piMinPt,piMinPt,dimuVtx,piSiHits,piSiHits,JpsiMinY,JpsiMaxY);
  } else if (nprompt=="ctau"){
     sprintf(cutstring,"X_ct>0.1 && X_Pt>%f && X_Rap>%f && X_Rap<%f && SoftPion_pt>%f && HardPion_pt>%f && DiMuonVtx>%f && SoftPion_Hits>=%f && HardPion_Hits>=%f && abs(JPsi_Rap)>=%f && abs(JPsi_Rap)<%f",XMinPt,-XMinY,XMinY,piMinPt,piMinPt,dimuVtx,piSiHits,piSiHits,JpsiMinY,JpsiMaxY);
  } else if (nprompt=="cowboy"){
     sprintf(cutstring,"DiMuonType>0 && X_Pt>%f && X_Rap>%f && X_Rap<%f && SoftPion_pt>%f && HardPion_pt>%f && DiMuonVtx>%f && SoftPion_Hits>=%f && HardPion_Hits>=%f && abs(JPsi_Rap)>=%f && abs(JPsi_Rap)<%f",XMinPt,-XMinY,XMinY,piMinPt,piMinPt,dimuVtx,piSiHits,piSiHits,JpsiMinY,JpsiMaxY);
  } else if (nprompt=="seagull"){
     sprintf(cutstring,"DiMuonType<0 && X_Pt>%f && X_Rap>%f && X_Rap<%f && SoftPion_pt>%f && HardPion_pt>%f && DiMuonVtx>%f && SoftPion_Hits>=%f && HardPion_Hits>=%f && abs(JPsi_Rap)>=%f && abs(JPsi_Rap)<%f",XMinPt,-XMinY,XMinY,piMinPt,piMinPt,dimuVtx,piSiHits,piSiHits,JpsiMinY,JpsiMaxY);
  } else if (nprompt=="XVtx2"){
     sprintf(cutstring,"X_Vtx>0.02 && X_Pt>%f && X_Rap>%f && X_Rap<%f && SoftPion_pt>%f && HardPion_pt>%f && DiMuonVtx>%f && SoftPion_Hits>=%f && HardPion_Hits>=%f && abs(JPsi_Rap)>=%f && abs(JPsi_Rap)<%f",XMinPt,-XMinY,XMinY,piMinPt,piMinPt,dimuVtx,piSiHits,piSiHits,JpsiMinY,JpsiMaxY);
  } else {
     sprintf(cutstring,"X_Pt>%f && X_Rap>%f && X_Rap<%f && SoftPion_pt>%f && HardPion_pt>%f && DiMuonVtx>%f && SoftPion_Hits>=%f && HardPion_Hits>=%f && abs(JPsi_Rap)>=%f && abs(JPsi_Rap)<%f",XMinPt,-XMinY,XMinY,piMinPt,piMinPt,dimuVtx,piSiHits,piSiHits,JpsiMinY,JpsiMaxY);
 }

 RooDataSet* dataSet = (RooDataSet*) FulldataSet->reduce(RooArgSet(*xmassFull),cutstring) ;

  
/////////////////////////
//// Fit
/////////////////////////

//// Define PDF

  // PsiPrime PDF
  RooRealVar   Psi_peak  ("Psi_mass"    ,"Psi_mass"  ,3.686, 3.65, 3.75 );
  RooRealVar   Psi_width1 ("Psi_#sigma1"    ,"Psi_#sigma1"  , 0.007, 0.001, 0.08 );
  RooRealVar   Psi_width2 ("Psi_#sigma2"    ,"Psi_#sigma2"  , 0.003, 0.001, 0.1 );
  RooRealVar   fpsi ("fpsi"    ,"Psi component fraction"  , 0.5, 0., 1.0 );
  RooGaussian  pdfPsi1("pdfPsi1_A", "pdfPsi1",  X_Mass, Psi_peak ,Psi_width1 ) ;
  RooGaussian  pdfPsi2("pdfPsi2_B", "pdfPsi2",  X_Mass, Psi_peak ,Psi_width2 ) ;
  // Voigtian
  RooRealVar   Psi_gamma ("Psi_#gamma"    ,"Psi_#gamma"  , 0.003, 0.001, 1.0 );
  RooVoigtian pdfPsiV("pdfPsiV", "pdfPsiV",  X_Mass, Psi_peak ,Psi_width2,Psi_gamma ) ;
  // 2 gaussians with different sigma
  RooAddPdf    pdfPsi("pdfPsi","SumPsi", RooArgList(pdfPsi1,pdfPsi2), fpsi ) ;
  // CB + gaussian
  RooRealVar   alpha("alpha"    ,"alpha"  , 1.7, 0.5, 3. );
  RooRealVar   n("n"    ,"n"  , 9., 3., 12. );
  RooCBShape   pdfCBall("CBall","Crystal Ball shape", X_Mass, Psi_peak, Psi_width2, alpha, n); // CB with sigma2
  RooAddPdf    pdfPsiCB("pdfPsiCB","SumPsiCB", RooArgList(pdfPsi1,pdfCBall), fpsi ) ; // CB with sigma2 + gaussian with sigma1
  // CB + Voigt
  RooAddPdf    pdfPsiCBV("pdfPsiCBV","SumPsiCBV", RooArgList(pdfPsiV,pdfCBall), fpsi ) ;


  // X PDF
  RooRealVar   X_peak  ("X_mass"    ,"X_mass"  ,3.87, 3.83, 3.91 );
  RooRealVar   X_width ("X_#sigma"    ,"X_#sigma"  , 0.006, 0.001, 0.01 );
  RooGaussian pdfX("pdfX", "pdfX", X_Mass, X_peak ,X_width ) ;
  // BKG PDF with Chebycev values
  RooRealVar c1("c1", "c1", 0.3,  -1., 1.0);
  RooRealVar c2("c2", "c2", -0.09,  -1.0, 1.0);
  RooRealVar c3("c3", "c3", 0.02,  -1.0, 1.0);
  RooRealVar c4("c4", "c4", 0.01,  -1.0, 1.0);
  RooChebychev bkg ("bkg","Background pdf",X_Mass,RooArgList(c1,c2));

  RooRealVar nsig1("Number_Psi2S", "signal fraction 1",70000.,0.,300000.);
  RooRealVar nsig2("Number_X3872", "signal fraction 2",5000.,0.,20000.);
  RooRealVar nbkg("nbkg", "bkg fraction ",1300000.,0.,3000000.);



//// Define Model
  RooArgList varlist(pdfX);  // X PDF
  if ( PsiPrimePDF =="CB"){
    varlist.add(pdfPsiCB);   // with CB+Gauss PDF for Psi(2S)
  } else if ( PsiPrimePDF =="CBVoigt"){
    varlist.add(pdfPsiCBV); // with CB+Voigt  for Psi(2S)
  } else if ( PsiPrimePDF =="Voigt"){
    varlist.add(pdfPsiV);   // with  Voigt  for Psi(2S)
  } else{
    varlist.add(pdfPsi);     // with 2Gauss PDF for Psi(2S)
  }
  varlist.add(bkg);   // Background PDF

  RooAddPdf model("model","model", varlist , RooArgList(nsig2,nsig1, nbkg) ) ;
  model.printCompactTree();

  RooFitResult* fitres = model.fitTo(*dataSet,NumCPU(3),FitOptions("mer"),Save(),Extended(kTRUE)); 

  gROOT->SetStyle("Plain");  


  float binsize= (RangeMassMax-RangeMassMin)/framebins; 
  float binsizeMeV= binsize*1000.;
  char binstring[4];
  sprintf(binstring,"%1.0f ",binsizeMeV);
  TString bintstring = TString(binstring);
  TCanvas * c=new TCanvas("c","c",960,800);
  c->SetLogy(Log);
  TGaxis::SetMaxDigits(3);
  RooPlot *frame = X_Mass.frame(framebins);
 
  //float linew=1.; 
  //float msize=linew;
  float linew=0.9;
  float msize=0.8;
  dataSet->plotOn(frame,MarkerStyle(21),MarkerSize(msize));
  model.plotOn(frame, Components(RooArgSet(bkg)),LineStyle(kDashed),LineColor(kBlue),LineWidth(linew),Range(RangeMassMin,RangeMassMax) );
  model.plotOn(frame,LineWidth(linew));
  frame->SetTitle(""); 
  frame->GetYaxis()->SetTitle("Candidates/ "+bintstring+" MeV"); 
  frame->GetXaxis()->SetLabelOffset(0.013);
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->SetMarkerStyle(21);
  frame->Draw();

  float xmintxt=0.25;
  float ymintxt=0.17;
  TPaveText *txt = new TPaveText(xmintxt,ymintxt,xmintxt,ymintxt,"NDC");
  txt->SetTextSize(0.04);
  txt->SetFillColor(0);
  txt->SetFillStyle(0);
  txt->SetLineStyle(2);
  txt->SetLineColor(0);
  txt->SetTextAlign(12);
  txt->SetTextFont(42);
  txt->AddText("CMS Preliminary");
  txt->Draw("same");
  float xmintxt2=0.6;
  float ymintxt2=0.26;
  TPaveText *txtE = new TPaveText(xmintxt2+0.01,ymintxt2,xmintxt2+0.01,ymintxt2,"NDC");
  txtE->SetTextSize(0.04);
  txtE->SetFillColor(0);
  txtE->SetFillStyle(0);
  txtE->SetLineStyle(2);
  txtE->SetLineColor(0);
  txtE->SetTextAlign(12);
  txtE->SetTextFont(42);
  txtE->AddText("#sqrt{s} = 7 TeV");
  txtE->Draw("same");
  TPaveText *txtL = new TPaveText(xmintxt2,ymintxt,xmintxt2,ymintxt,"NDC");
  txtL->SetTextSize(0.04);
  txtL->SetFillColor(0);
  txtL->SetFillStyle(0);
  txtL->SetLineStyle(2);
  txtL->SetLineColor(0);
  txtL->SetTextAlign(12);
  txtL->SetTextFont(42);
  txtL->AddText("#int L dt = 896 pb^{ - 1}");
  txtL->Draw("same");

  c->Print("XMass_"+name+"_"+cutname+"_"+PsiPrimePDF+".png");


// Print Fit results
    cout<<" ==================================="<<endl;
    cout<<" ===> Used cuts: "<<cutstring<<endl;
    if (PsiPrimePDF =="CB") cout<<" ===> Psi(2S) with CB + Gauss"<<endl;
    if (PsiPrimePDF =="CBVoigt") cout<<" ===> Psi(2S) with CB + Voigt"<<endl;
    if (PsiPrimePDF =="Voigt") cout<<" ===> Psi(2S) with Voigt"<<endl;  
    cout<<" ==================================="<<endl;
    fitres->Print();
    //RooArgSet* params = model.getVariables() ;  
    //params->Print("v") ;

   double Nbkg = nbkg.getVal();
   double NSignal = nsig2.getVal();
   double NSignalPsi = nsig1.getVal();
   double minRangeX = X_peak.getVal() -  (X_width.getVal()*2.0) ;
   double maxRangeX = X_peak.getVal() +  (X_width.getVal()*2.0) ;

   X_Mass.setRange("integral",minRangeX, maxRangeX) ;
   RooAbsReal* N  = bkg.createIntegral(RooArgSet(X_Mass),NormSet(X_Mass),Range("integral"));
   RooAbsReal* S  = pdfX.createIntegral(RooArgSet(X_Mass),NormSet(X_Mass),Range("integral")); 
   double SN=(NSignal*S->getVal())/(Nbkg*N->getVal());


   ofstream fout("FitResult_"+name+"_"+cutname+"_"+PsiPrimePDF+".txt");
   fout<<"-------------------------------"<<endl;
   fout<<" ===> Used cuts: "<<cutstring<<endl;
   fout << "chi^2 = " << frame->chiSquare() << endl ;
   fout<<" N_X = "<<NSignal<<"+/-"<<nsig2.getError()<<"  stat error ="<<100.*(nsig2.getError()/NSignal)<<"%"<<endl;
   fout<<" N_Psi = "<<NSignalPsi<<"+/-"<<nsig1.getError()<<endl;
   fout<<" N_bkg = "<<Nbkg<<"+/-"<<nbkg.getError()<<endl;
   fout<<" X(3872) S/N ="<<SN <<endl;
   fout<<"-------------------------------"<<endl;

    double minFixX = 3.86;
    double maxFixX = 3.88;
    X_Mass.setRange("integral",minFixX, maxFixX) ;
    RooAbsReal* NFix  = bkg.createIntegral(RooArgSet(X_Mass),NormSet(X_Mass),Range("integral"));
    RooAbsReal* SFix  = pdfX.createIntegral(RooArgSet(X_Mass),NormSet(X_Mass),Range("integral"));

//    TCanvas* cX = new TCanvas("cX","cX",900,300) ;
    TH1F* hMass= new TH1F("XMass","",80,3.6,4.);
    for(unsigned int i=0; i < dataSet->numEntries(); i++)
    {
        const RooArgSet* row = dataSet->get(i);
        RooRealVar* xmass = (RooRealVar*) row->find("X_Mass");
        hMass->Fill(xmass->getVal());
    }
    int bmin =hMass->FindBin(minFixX);
    int bmax =hMass->FindBin(maxFixX);
    double TotalInRange = hMass->Integral(bmin+1,bmax);
    hMass->GetXaxis()->SetRange(bmin+1,bmax);
    hMass->Draw();
//    cX->Print("XRange_"+name+"_"+cutname+"_"+PsiPrimePDF+".png");

    double BkgInRange=Nbkg*NFix->getVal(); // bkg from fit
    double SignalInRange=TotalInRange-BkgInRange;
    double SignalInRangeError=nbkg.getError()*NFix->getVal()+sqrt(TotalInRange);

    fout<<"-------------------------------"<<endl;
    fout<<" In Fixed mass window "<<minFixX<<"-"<<maxFixX<<"  (bin "<<bmin+1<<"-"<<bmax<<")"<<endl;
    fout<<"   Bkg in mass window (Fit) ="<<BkgInRange<<endl;
    fout<<"   Total in mass window (counting)="<<TotalInRange<<endl;
    fout<<"   Signal in mass window (counting - Bkg Fit)="<<SignalInRange<<" +/- "<<SignalInRangeError<<endl;
    fout<<" X(3872) Signal (counting) = "<<SignalInRange<<" +/- "<<SignalInRangeError<<endl;
    fout<<" X(3872) Noise= "<<BkgInRange<<endl;
    fout<<" X(3872) S/N ="<<SignalInRange/BkgInRange <<endl;
    fout<<" X(3872) Signal(Fit)="<<nsig2.getVal()*SFix->getVal()<<" +/- "<<nsig2.getError()*SFix->getVal()<<endl;
    fout<<"-------------------------------"<<endl;


  ////
  //// Residual histogram
  ////
  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid = frame->residHist() ;
  // Create a new frame to draw the residual distribution and add the distribution to the frame
  RooPlot* frame2 = X_Mass.frame(Title("Residual Distribution")) ;
  frame2->addPlotable(hresid,"P") ;
  TCanvas* cresidual = new TCanvas("res","res",900,300) ;
  frame2->Draw() ;
  cresidual->Print("Residual_"+name+"_"+cutname+"_"+PsiPrimePDF+".png");

  ////
  //// Pull histogram
  ////
  RooHist* hpull = frame->pullHist() ;
  RooPlot* frame3 = X_Mass.frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"P") ;
  TCanvas* cpull = new TCanvas("pull","pull",900,300) ;
  frame3->Draw() ;
  cpull->Print("Pull_"+name+"_"+cutname+"_"+PsiPrimePDF+".png");



}
