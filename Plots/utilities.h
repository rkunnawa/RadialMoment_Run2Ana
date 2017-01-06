#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>
#include <TRandom.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TCanvas.h"
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TBox.h>
#include "TF1.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include <stdio.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>

using namespace std;


// TAA 

TGraphErrors *tTAAerr[6]={0};
TGraphErrors *tTAAerrNpart=0;

const int nbins_cent = 6;
double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};

double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
double npart[nbins_cent+1] = {381.829, 329.41, 224.28, 108.12, 42.04, 11.43, 112.9};

int style[nbins_cent+1] = {20, 24, 21, 25, 33, 27, 34};
int color[nbins_cent+1] = {2, 4, 6, 7, 8, 9, 1};
int msize = 2;
int shade[] = {3004, 3006, 3005};
int mstyle[] = {29, 33, 34};
int lstyle[] = {3, 1, 2};
int mcolor[] = {kRed, kBlack, kBlue};

/* const int nbins_pt = 30; */
/* const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12,   15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300,   330, 362, 395 }; */

const char * algo = (char*) "Pu";
const char *jet_type = (char*) "PF";
char* etaWidth = (char*) "20_eta_20";
double deltaEta = 4.0;
const char * ptbins = (char*)"finebins"; /* finebins, atlaspTbins, atlasRcpbins */

int kRegDraw_R2[nbins_cent] = {4, 5, 5, 5, 4, 4}; 
int kRegDraw_R3[nbins_cent] = {5, 5, 8, 5, 4, 4}; 
int kRegDraw_R4[nbins_cent] = {5, 4, 5, 6, 4, 4}; 
/* int kRegDraw_R2[nbins_cent] = {4, 5, 5, 5, 4, 4}; */
/* int kRegDraw_R3[nbins_cent] = {5, 7, 8, 5, 6, 5}; */
/* int kRegDraw_R4[nbins_cent] = {5, 7, 5, 8, 7, 7}; */

int kRegDraw_R2_PP = 6;
int kRegDraw_R3_PP = 6;
int kRegDraw_R4_PP = 6;
const int SVD_kReg_num = 5;
int SVD_kReg_values[SVD_kReg_num] = {4, 5, 6, 7, 8};
const int SVD_kReg_PbPb_R2[nbins_cent] = {0, 1, 1, 1, 0, 0};
const int SVD_kReg_PbPb_R3[nbins_cent] = {1, 1, 4, 1, 0, 0};
const int SVD_kReg_PbPb_R4[nbins_cent] = {1, 0, 1, 2, 0, 0};
const int SVD_kReg_PP_R2 = 2;
const int SVD_kReg_PP_R3 = 2;
const int SVD_kReg_PP_R4 = 2;
int ErrorkReg[nbins_cent] = {10, 10, 10, 10, 8, 7};


const int unfoldingCut_R2 = 55;
const int unfoldingCut_R3 = 55;
const int unfoldingCut_R4 = 55;

static const int nbins_short = 7;
static const double boundaries_short[nbins_short+1] = {
  43, 49, 56, 64, 74, 84, 97, 114
};

static const int nbins_yaxian_large = 29;
static const double ptbins_yaxian_large[nbins_yaxian_large+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638,790,967};


static const int nbins_yaxian = 21;
static const double boundaries_yaxian[nbins_yaxian+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395};

static const int nbins_nlo_comp = 12;
static const double boundaries_nlo_comp[nbins_nlo_comp+1] = {64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300};

static const int dir = 50;


const double ptbins_nlo[] = { 
  3, 4, 5, 7, 9, 12,  
  15, 18, 21, 24, 28, 
  32, 37, 43, 49, 56, 
  64, 74, 84, 97, 114, 
  133, 153, 174, 196, 
  220, 245, 272, 300,  
  330, 362, 395, 430, 
  468, 507, 548, 592, 
  638, 686, 1000  
}; 
const int nbins_nlo = sizeof(ptbins_nlo)/sizeof(double) - 1; 
/* const int nbins_pt_truncated = 12; */
/* double boundaries_pt_truncated[nbins_pt_truncated+1] = {64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300};  */
const double ptbins_ana[] = {50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
const int nbins_ana = sizeof(ptbins_ana)/sizeof(double) -1;

const double qg_sys[nbins_cent][nbins_ana] =
  //50, 60 ,  70 ,  80 ,  90 , 100 , 110 , 130 , 150 , 170 , 190 , 210 , 240 , 270 , 300 
  {{0.12, 0.12, 0.13, 0.12, 0.12, 0.12, 0.10, 0.10, 0.08, 0.08, 0.08, 0.07, 0.07, 0.07},
   {0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.10, 0.10, 0.08, 0.08, 0.08, 0.07, 0.07, 0.07},
   {0.12, 0.12, 0.10, 0.09, 0.09, 0.08, 0.08, 0.07, 0.07, 0.07, 0.07, 0.06, 0.06, 0.06},
   {0.12, 0.10, 0.08, 0.08, 0.07, 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05, 0.05},
   {0.12, 0.08, 0.07, 0.07, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05, 0.05},
   {0.12, 0.08, 0.07, 0.07, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05, 0.05}
  };


const double boundaries_powheg_Uncert[] = {56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300};
const int nbins_powheg_Uncert = sizeof(boundaries_powheg_Uncert)/sizeof(double)-1;


const double boundaries_powheg[] = {3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 1000};
const int nbins_powheg = sizeof(boundaries_powheg)/sizeof(double)-1;


const double ptbins_ALICEComp[] = {30, 40, 50, 60, 70, 80, 90, 100, 110, 120};
const int nbins_ALICEComp = sizeof(ptbins_ALICEComp)/sizeof(double) -1;

const double ptbins_ATLASComp[] = {30, 40, 50, 60, 80, 100, 120, 160, 200, 250, 320};
const int nbins_ATLASComp = sizeof(ptbins_ATLASComp)/sizeof(double) -1;

const double ptbins_long[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000}; 
/* const double ptbins_long[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168, 170, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248, 250, 252, 254, 256, 258, 260, 262, 264, 266, 268, 270, 272, 274, 276, 278, 280, 282, 284, 286, 288, 290, 292, 294, 296, 298, 300, 302, 304, 306, 308, 310, 312, 314, 316, 318, 320, 322, 324, 326, 328, 330, 332, 334, 336, 338, 340, 342, 344, 346, 348, 350, 352, 354, 356, 358, 360, 362, 364, 366, 368, 370, 372, 374, 376, 378, 380, 382, 384, 386, 388, 390, 392, 394, 396, 398, 400, 402, 404, 406, 408, 410, 412, 414, 416, 418, 420, 422, 424, 426, 428, 430, 432, 434, 436, 438, 440, 442, 444, 446, 448, 450, 452, 454, 456, 458, 460, 462, 464, 466, 468, 470, 472, 474, 476, 478, 480, 482, 484, 486, 488, 490, 492, 494, 496, 498, 500, 502, 504, 506, 508, 510, 512, 514, 516, 518, 520, 522, 524, 526, 528, 530, 532, 534, 536, 538, 540, 542, 544, 546, 548, 550, 552, 554, 556, 558, 560, 562, 564, 566, 568, 570, 572, 574, 576, 578, 580, 582, 584, 586, 588, 590, 592, 594, 596, 598, 600, 602, 604, 606, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 628, 630, 632, 634, 636, 638, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660, 662, 664, 666, 668, 670, 672, 674, 676, 678, 680, 682, 684, 686, 688, 690, 692, 694, 696, 698, 700, 702, 704, 706, 708, 710, 712, 714, 716, 718, 720, 722, 724, 726, 728, 730, 732, 734, 736, 738, 740, 742, 744, 746, 748, 750, 752, 754, 756, 758, 760, 762, 764, 766, 768, 770, 772, 774, 776, 778, 780, 782, 784, 786, 788, 790, 792, 794, 796, 798, 800, 802, 804, 806, 808, 810, 812, 814, 816, 818, 820, 822, 824, 826, 828, 830, 832, 834, 836, 838, 840, 842, 844, 846, 848, 850, 852, 854, 856, 858, 860, 862, 864, 866, 868, 870, 872, 874, 876, 878, 880, 882, 884, 886, 888, 890, 892, 894, 896, 898, 900, 902, 904, 906, 908, 910, 912, 914, 916, 918, 920, 922, 924, 926, 928, 930, 932, 934, 936, 938, 940, 942, 944, 946, 948, 950, 952, 954, 956, 958, 960, 962, 964, 966, 968, 970, 972, 974, 976, 978, 980, 982, 984, 986, 988, 990, 992, 994, 996, 998, 1000}; */
const int nbins = sizeof(ptbins_long)/sizeof(double) -1;

const double boundaries_pt_truncated[] = {60, 70, 80, 90, 110, 130, 150, 170, 190, 220, 240, 270, 300};
const int nbins_pt_truncated = sizeof(boundaries_pt_truncated)/sizeof(double) -1;

const double boundaries_pt[] = {60, 70, 80, 90, 110, 130, 150, 170, 190, 220, 240, 270, 300};
const int nbins_pt = sizeof(boundaries_pt)/sizeof(double) -1;

const int nbins_pt_truncated_gen = 23;   
double boundaries_pt_truncated_gen[nbins_pt_truncated_gen+1] = {21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395};

const int nbins_pt_good = 23;  
const double boundaries_pt_good[nbins_pt_good+1] = { 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 272, 300, 500};

const double ptcut_high[] = {300.0, 270.0, 300.0, 300.0, 240.0, 170.0, 300.0};
const double ptcut_low[] = {50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0};
const double draw_pt_high[] = {299.0, 269.0, 299.0, 299.0, 239.0, 169.0, 299.0};
const double draw_pt_low[] = {60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0};
// smaller pT bins for the MC 
/* static const int nbins_pt = 17; */
/* static const double boundaries_pt[nbins_pt+1] = { 32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330}; */

// this is the cms full analysis pT bins 
/* const int nbins_pt = 32;   */
/* const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 501}; */

/* const double boundaries_pt[] = {60, 70, 80, 90, } */
/* const int nbins_pt = ; */


/* const int nbins_pt_fine = 285; */
/* const double boundaries_pt_fine[nbins_pt_fine+1] = {15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300}; */

// truncated analysis bins 
/* const int nbins_pt = 20; */
/* double boundaries_pt[nbins_pt+1] = {32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395}; */

// this is the atlas pT binning for the spectra. 
static const int nbins_atlas_pt = 12; 
static const double boundaries_atlas_pt[nbins_atlas_pt+1] = {31., 39., 50., 63., 79., 100., 125., 158., 199., 251., 316., 398., 501}; 

static const int nbins_alice_pt = 9;
static const double boundaries_alice_pt[nbins_alice_pt+1] = {30, 40, 50, 60, 70, 80, 90, 100, 110, 120};


// this is the atlas pT binning for the Rcp plots. 
/* static const int nbins_pt = 12;   */
/* static const double boundaries_pt[nbins_pt+1] = {38.36, 44.21, 50.94, 58.7, 67.64 , 77.94 , 89.81, 103.5, 119.3, 137.4 , 158.3, 182.5,  210.3};   */

void doTrigCorr(TH1F *h1, TH1F *h2)
{
  for(int ix = 1; ix<=h2->GetNbinsX(); ++ix){
    int xbin = h1->FindBin(h2->GetBinCenter(ix));
    double val = h1->GetBinContent(xbin);
    val = val * 1./h2->GetBinContent(ix);
    h1->SetBinContent(xbin, val);
    double err = h1->GetBinError(xbin);
    err = err * 1./h2->GetBinContent(ix);
    h1->SetBinError(xbin, err);
  }
}

void doCorr(TH1F *h1, TH1F *h2, int doRecp = 0)
{
  for(int ix = 1; ix<=h2->GetNbinsX(); ++ix){
    int xbin = h1->FindBin(h2->GetBinCenter(ix));
    double val = h1->GetBinContent(xbin);
    if(doRecp == 0)val = val * h2->GetBinContent(ix);
    if(doRecp == 1)val = val * 1./h2->GetBinContent(ix);
    h1->SetBinContent(xbin, val);
    double err = h1->GetBinError(xbin);
    if(doRecp == 0)err = err * h2->GetBinContent(ix);
    if(doRecp == 1)err = err * 1./h2->GetBinContent(ix);
    h1->SetBinError(xbin, err);
  }
}


void SetUnfoldBins1D(TH1F *h, float minpt, float maxpt)
{
  int lbin = h->GetXaxis()->FindBin(minpt);
  int hbin = h->GetXaxis()->FindBin(maxpt);
  for(int ix=1; ix<=h->GetNbinsX(); ix++){
    if(ix > lbin && ix < hbin)continue;
    h->SetBinContent(ix,0);
    h->SetBinError(ix,0);
  }
}

void SetUnfoldBins2D(TH2F *h, float minGenPt, float maxGenPt, float minRecPt, float maxRecPt)
{
  int lgen_bin = h->GetXaxis()->FindBin(minGenPt);
  int hgen_bin = h->GetXaxis()->FindBin(maxGenPt);
  int lrec_bin = h->GetYaxis()->FindBin(minRecPt);
  int hrec_bin = h->GetYaxis()->FindBin(maxRecPt);

  for(int ix=1; ix<=h->GetNbinsX(); ix++){
    for(int iy=1; iy<=h->GetNbinsY(); iy++){
      if( ( ix > lgen_bin && ix < hgen_bin )
	  && ( iy > lrec_bin && iy < hrec_bin ) )continue;
      h->SetBinContent(ix, iy, 0);
      h->SetBinError(ix, iy, 0);
    }
  }
}

void Truncate1D(TH1 * hInput, TH1 * hOutput){

  for(int k = 1; k<=hOutput->GetNbinsX(); ++k){
    int binx = hInput->FindBin(hOutput->GetBinCenter(k));
    Double_t val = hInput->GetBinContent(binx);
    Double_t val_Err = hInput->GetBinError(binx);
    hOutput->SetBinContent(k, val);
    hOutput->SetBinError(k, val_Err);
  }
  
}

void Truncate2D(TH2 * hInput, TH2 * hOutput){

  for(int k = 1; k<=hOutput->GetNbinsX(); ++k){
    for(int l = 1; l<=hOutput->GetNbinsY(); ++l){      
      int binx = hInput->GetXaxis()->FindBin(hOutput->GetXaxis()->GetBinCenter(k));
      int biny = hInput->GetYaxis()->FindBin(hOutput->GetYaxis()->GetBinCenter(l));
      Double_t val = hInput->GetBinContent(binx, biny);
      Double_t val_Err = hInput->GetBinError(binx, biny);
      hOutput->SetBinContent(k, l, val);
      hOutput->SetBinError(k, l, val_Err);
    }
  }
  
}

// 
void transferErrorBars(TH1 *hinput, TH1 *houtput){

  if(hinput->GetNbinsX() != houtput->GetNbinsX())
    cout<<"ERROR, both histograms are not in the same bins"<<endl;
  
  for(int bin = 1; bin <= hinput->GetNbinsX(); ++bin)
    houtput->SetBinError(bin, hinput->GetBinError(bin));    
  
}


TH1 *Truncate1D(TH1 * hInput, Int_t nbins, Int_t xLow, Int_t xHigh){

  TH1F * hOutput = new TH1F("hOutput","",nbins,xLow, xHigh);
  for(int j = 1; j<=hInput->GetNbinsX(); j++){
    if( j >= xLow && j <= xHigh){ 
    hOutput->SetBinContent(j-xLow, hInput->GetBinContent(j));
    hOutput->SetBinError(j-xLow, hInput->GetBinError(j));
    }
  }
  return hOutput; 
}

TH2 *Truncate2D(TH2 * hInput, Int_t nbins_X, Int_t xLow, Int_t xHigh, Int_t nbins_Y, Int_t yLow, Int_t yHigh){

  TH2F * hOutput = new TH2F("hOutput","", nbins_X, xLow, xHigh, nbins_Y, yLow, yHigh);
  for(int k = 1; k<=hInput->GetNbinsY(); k++){
    for(int j = 1; j<=hInput->GetNbinsX(); j++){
      if( j >= xLow && j <= xHigh && k >= yLow && k <= yHigh){ 
	hOutput->SetBinContent(j-xLow,k-yLow, hInput->GetBinContent(j,k));
	hOutput->SetBinError(j-xLow,k-yLow, hInput->GetBinError(j,k));
      }
    }
  }

  return hOutput; 

}


class SysData
{
 public:
  SysData() {
    for (int i=0;i<=nbins_cent;i++) {
      hSys[i]     = new TH1F(Form("hSys_cent%d",i), Form("Totalsys_cent%d",i),nbins_ana, ptbins_ana);
      hSysGeneral[i]= new TH1F(Form("hSysGeneral_cent%d",i), Form("TotalsysGeneral_cent%d",i),nbins_ana, ptbins_ana);
      hSysJEC[i]  = new TH1F(Form("hSysJEC_cent%d",i), Form("JECsys_cent%d",i),nbins_ana, ptbins_ana);
      hSysBKG[i]  = new TH1F(Form("hSysBKG_cent%d",i), Form("BKGsys_cent%d",i),nbins_ana, ptbins_ana);
      hSysJER[i]  = new TH1F(Form("hSysJER_cent%d",i), Form("JERsys_cent%d",i),nbins_ana, ptbins_ana);
      //hSysEff[i]  = new TH1F(Form("hSysEff_cent%d",i), Form("Effsys_cent%d",i),nbins_ana, ptbins_ana);
      hSysSmear[i]  = new TH1F(Form("hSysSmear_cent%d",i), Form("Smearsys_cent%d",i),nbins_ana, ptbins_ana);
      hSysIter[i] = new TH1F(Form("hSysIter_cent%d",i), Form("Itersys_cent%d",i),nbins_ana, ptbins_ana);
      hSysPrior[i] = new TH1F(Form("hSysPrior_cent%d",i),Form("PriorUnfolding_sys_cent%d",i), nbins_ana, ptbins_ana);
      hSysJetID[i] = new TH1F(Form("hSysJetID_cent%d",i), Form("JetID_sys_cent%d",i), nbins_ana, ptbins_ana);
      hSys[i]->SetLineColor(kGray);
      hSysJEC[i]->SetLineColor(4);
      hSysBKG[i]->SetLineColor(5);      
      hSysJER[i]->SetLineColor(5);
      hSysSmear[i]->SetLineColor(kGreen+1);
      hSysIter[i]->SetLineColor(2);
      hSysJetID[i]->SetLineColor(kGreen+2);
      hSysPrior[i]->SetLineColor(kGreen+3);
    }  
  }
  TH1F *hSys[nbins_cent+1];
  TH1F *hSysGeneral[nbins_cent+1];
  TH1F *hSysJEC[nbins_cent+1];
  TH1F *hSysBKG[nbins_cent+1];
  TH1F *hSysJER[nbins_cent+1];
  //TH1F *hSysEff[nbins_cent+1];
  TH1F *hSysSmear[nbins_cent+1];
  TH1F *hSysIter[nbins_cent+1];
  TH1F *hSysNoise[nbins_cent+1];
  TH1F *hSysJetID[nbins_cent+1];
  TH1F *hSysPrior[nbins_cent+1];
  bool isBayesian = false; 
	
  void calcTotalSys(int i, int isRAA = 0) {
    TF1 *fNoise = new TF1("f","1+0.3*0.16*abs(1-([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x))");
    fNoise->SetParameters(0.9521,0.001105,-9.397e-6,3.32e-8,-5.618e-11);
    hSysNoise[i] = functionHist(fNoise,hSys[i],Form("hSysNoise_cent%d",i));
    hSysNoise[i]->SetName(Form("hSysNoise_cent%d",i));
    hSysNoise[i]->SetLineColor(6);
    for (int j=1;j<=hSys[i]->GetNbinsX();j++) {
      //double effSys = 0.025;
      double jetidSys = 0.02;
      if(i == 6) jetidSys = 0.0;
      //hSysEff[i]->SetBinContent(j,1+effSys);
      //hSysSmear[i]->SetBinContent(j,1.02);
      hSysJetID[i]->SetBinContent(j, 1+jetidSys);
      double JERSys = 0.0;
      JERSys = hSysJER[i]->GetBinContent(j)-1;

      double JECSys = hSysJEC[i]->GetBinContent(j)-1;
      // additional systematics in scale due to HIN tracking. // already added in the histograms and they will be plotted. 
      /* if(i!=6){ */
      /* 	if(hSys[i]->GetBinCenter(j) <= 80) JECSys = JECSys + 0.065; */
      /* 	if(hSys[i]->GetBinCenter(j) > 80 && hSys[i]->GetBinCenter(j) <= 110) JECSys = JECSys + 0.068; */
      /* 	if(hSys[i]->GetBinCenter(j) > 110) JECSys = JECSys + 0.04; */
      /* } */
      double SmearSys = hSysSmear[i]->GetBinContent(j)-1;
      if(i == 6) SmearSys = 0.0;
      double IterSys = hSysIter[i]->GetBinContent(j)-1;
      
      double BKGSys  = 0.0;
      if(i!=6) {
      	BKGSys  = hSysBKG[i]->GetBinContent(j)-1; 
      /* 	if(hSys[i]->GetBinCenter(j)<=70) BKGSys = BKGSys + 0.07; */
      /* 	if(hSys[i]->GetBinCenter(j)>70 && hSys[i]->GetBinCenter(j)<=80) BKGSys = BKGSys + 0.04; */
      /* 	if(hSys[i]->GetBinCenter(j)>80 && hSys[i]->GetBinCenter(j)<=90) BKGSys = BKGSys + 0.03; */
      /* 	if(hSys[i]->GetBinCenter(j)>90 && hSys[i]->GetBinCenter(j)<=100) BKGSys = BKGSys + 0.025; */
      /* 	if(hSys[i]->GetBinCenter(j)>100 && hSys[i]->GetBinCenter(j)<=110) BKGSys = BKGSys + 0.02; */
      } 
      double PriorSys = 0.0;
      if(isBayesian == true && i == 6) {
	PriorSys = 0.03;
	hSysPrior[i]->SetBinContent(j, 1+PriorSys);
      }else if(isBayesian == true && i !=6) {
	PriorSys = 0.015;
	hSysPrior[i]->SetBinContent(j, 1+PriorSys);
      }else PriorSys = hSysPrior[i]->GetBinContent(j)-1;
      
      double NoiseSys = hSysNoise[i]->GetBinContent(j)-1;
      if(isRAA) NoiseSys = 0;
      //cout <<effSys<<" "<<JECSys<<" "<<IterSys<<endl;
      double totalSys = sqrt( //effSys * effSys +
			      PriorSys * PriorSys +
			      JECSys * JECSys +
			      JERSys * JERSys +
			      SmearSys * SmearSys +
			      NoiseSys * NoiseSys +
			      IterSys* IterSys +
			      jetidSys * jetidSys +
			      BKGSys * BKGSys
			      );
      hSys[i]->SetBinContent(j,totalSys+1);	
      hSys[i]->SetLineWidth(2);				
    }
  }
  
  void calcTotalSysNoUnfolding(int i) {
    TF1 *fNoise = new TF1("f","1+0.3*0.16*abs(1-([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x))");
    fNoise->SetParameters(0.9521,0.001105,-9.397e-6,3.32e-8,-5.618e-11);
    hSysNoise[i] = functionHist(fNoise,hSys[i],Form("hSysNoise_cent%d",i));
    hSysNoise[i]->SetName(Form("hSysNoise_cent%d",i));
    hSysNoise[i]->SetLineColor(6);
    for (int j=1;j<=hSysGeneral[i]->GetNbinsX();j++) {
      //double effSys = 0.01;
      double jetidSys = 0.02;
      //hSysEff[i]->SetBinContent(j,1+effSys);
      hSysSmear[i]->SetBinContent(j,1.02);
      hSysJetID[i]->SetBinContent(j, 1+jetidSys);
      /* if(i == 0) hSysJEC[i]->SetBinContent(j, 1.12); */
      /* if(i == 1) hSysJEC[i]->SetBinContent(j, 1.10); */
      /* if(i == 2) hSysJEC[i]->SetBinContent(j, 1.08); */
      /* if(i == 3) hSysJEC[i]->SetBinContent(j, 1.05); */
      /* if(i == 4) hSysJEC[i]->SetBinContent(j, 1.03); */
      /* if(i == 5) hSysJEC[i]->SetBinContent(j, 1.02); */
      double JECSys = hSysJEC[i]->GetBinContent(j)-1;
      double SmearSys = hSysSmear[i]->GetBinContent(j)-1;
      double NoiseSys = hSysNoise[i]->GetBinContent(j)-1; 
      double totalSys = sqrt( //effSys * effSys +
			      JECSys * JECSys +
			      SmearSys * SmearSys +
			      NoiseSys * NoiseSys +
			      jetidSys * jetidSys
			      );
      hSysGeneral[i]->SetBinContent(j,totalSys+1);	
      hSysGeneral[i]->SetLineWidth(2);				
    }
  }
  
  void Draw(TH1F *h,int i, int color, int style = 0) {

    Int_t beginning = h->FindBin(ptcut_low[i])+2; Int_t end = h->FindBin(ptcut_high[i]);
    
    //if(i==0){ beginning = h->FindBin(100); }
    //if(i==5) { end = h->FindBin(180);}
    //if(i==4) { end = h->FindBin(240);}

    const Int_t N = end-beginning;
    cout<<"centrality bin = "<<i<<endl;
    cout<<"beginning = "<<beginning<<", end = "<<end<<" and N = "<<N<<endl;
    double x[N];
    double y[N];
    double ex[N];
    double ey[N];

    int counter = 0;
    for (int j=beginning;j<end;j++) {
      double val = h->GetBinContent(j);
      double err = hSys[i]->GetBinContent(j)-1;
      cout << "Sys Value Check" <<val<<" "<<err<<" "<<h->GetBinLowEdge(j)<<" "<<val*(1-err)<<" "<<h->GetBinLowEdge(j+1)<<" "<<val*(1+err)<<endl;
      x[counter] = h->GetBinCenter(j);
      y[counter] = h->GetBinContent(j);
      ex[counter] = h->GetBinLowEdge(j);
      ey[counter] = val*(1+err) - val;
      cout<<"systematic tgraph check: N = "<<counter<<", x = "<<x[counter]<<", y = "<<y[counter]<<", ex = "<<ex[counter]<<", ey = "<<ey[counter]<<endl;
      
      /* TBox *b = new TBox(h->GetBinLowEdge(j),val*(1-err),h->GetBinLowEdge(j+1),val*(1+err)); */
      /* //b->SetFillColor(kGray); */
      /* b->SetFillStyle(style); */
      /* //b->SetLineColor(kGray); */
		             
      /* //\***********For Gunther's Color Systematics Band Peference */
      /* b->SetFillColor(color); */
      /* b->SetLineColor(color); */
      /* b->Draw(); */
      counter++;
      
      if(j == beginning){
	TLine * lineup = new TLine(h->GetBinCenter(j), val*(1+err), h->GetBinLowEdge(j+1), val*(1+err));
	lineup->SetLineColor(color);
	lineup->SetLineStyle(style);
	lineup->Draw();
	TLine * linedown = new TLine(h->GetBinCenter(j), val*(1-err), h->GetBinLowEdge(j+1), val*(1-err));
	linedown->SetLineColor(color);
	linedown->SetLineStyle(style);
	linedown->Draw();
      }
      else if(j == end-1){
	TLine * lineup = new TLine(h->GetBinLowEdge(j), val*(1+err), h->GetBinCenter(j), val*(1+err));
	lineup->SetLineColor(color);
	lineup->SetLineStyle(style);
	lineup->Draw();
	TLine * linedown = new TLine(h->GetBinLowEdge(j), val*(1-err), h->GetBinCenter(j), val*(1-err));
	linedown->SetLineColor(color);
	linedown->SetLineStyle(style);
	linedown->Draw();

      }
      else {
	TLine * lineup = new TLine(h->GetBinLowEdge(j), val*(1+err), h->GetBinLowEdge(j+1), val*(1+err));
	lineup->SetLineColor(color);
	lineup->SetLineStyle(style);
	lineup->Draw();
	TLine * linedown = new TLine(h->GetBinLowEdge(j), val*(1-err), h->GetBinLowEdge(j+1), val*(1-err));
	linedown->SetLineColor(color);
	linedown->SetLineStyle(style);
	linedown->Draw();
      }
      
    }

    TGraphErrors* ge = new TGraphErrors(N, x, y, ex, ey);
    ge->SetFillColor(color);
    ge->SetFillStyle(style);
    /* ge->Draw("3"); */
    
  }
  
  void DrawTGraph(TGraphErrors *h,int i) {
    double xv;
    double val;
    for(int j=0;j<h->GetN();j++){	
      h->GetPoint(j,xv,val);
      double err = hSysGeneral[i]->GetBinContent(j+1)-1;
      //cout <<"value" <<val<<" "<<err<<" "<<hSysGeneral[i]->GetBinLowEdge(j+1)<<" "<<val*(1-err)<<" "<<hSysGeneral[i]->GetBinLowEdge(j+2)<<" "<<val*(1+err)<<endl;
      TBox *b_ = new TBox(hSysGeneral[i]->GetBinLowEdge(j+1),val*(1-err),hSysGeneral[i]->GetBinLowEdge(j+2),val*(1+err));
      b_->SetFillColor(kGray);
      b_->SetFillStyle(1001);
      b_->SetLineColor(kGray);
      b_->Draw();
    }
  }
	
  void DrawUnfoErr(TH1F *h,int i) {
    for (int j=1;j<=hSysIter[i]->GetNbinsX();j++) {
      double val = h->GetBinContent(j);
      double err = hSysIter[i]->GetBinContent(j)-1;
      //cout <<"value" << val<<" "<<err<<" "<<h->GetBinLowEdge(j)<<" "<<val*(1-err)<<" "<<h->GetBinLowEdge(j+1)<<" "<<val*(1+err)<<endl;
      TBox *b = new TBox(h->GetBinLowEdge(j),val*(1-err),h->GetBinLowEdge(j+1),val*(1+err));
      b->SetFillColor(29);
      b->SetFillStyle(1001);
      b->SetLineColor(29);
      b->Draw();
    }
  }
  
  
  void DrawNpartSys(double yvalue,int i,double xvalue, int binno, int style = 1001, int color = kGray) {
    double yerrorNpart[6]= {0.0409, 0.0459,0.0578,0.0944, 0.143, 0.176 };
    double err = hSys[i]->GetBinContent(hSys[i]->FindBin(binno))-1;
    TBox *b = new TBox(xvalue-10.,yvalue*(1-err),xvalue+10.,yvalue*(1+err));
    //cout << "value " << yvalue<<" err   "<<err<<" xvalue  "<<xvalue<<" "<<yvalue*(1-err)<<" "<<yvalue*(1+err)<<endl;
    b->SetFillColor(color);
    b->SetFillStyle(style);
    b->SetLineColor(color);
    /* b->Draw(); */
    
    TLine * lineup = new TLine(xvalue-10.,yvalue*(1+err),xvalue+10.,yvalue*(1+err));
    lineup->SetLineColor(color);
    lineup->SetLineStyle(style);
    lineup->Draw();
    TLine * linedown = new TLine(xvalue-10.,yvalue*(1-err),xvalue+10.,yvalue*(1-err));
    linedown->SetLineColor(color);
    linedown->SetLineStyle(style);
    linedown->Draw();
    		
  }
  
  
  void DrawComponent(int i, int isPP = 0, int isRAA = 0) {
    calcTotalSys(i);
    TH1D *h = new TH1D(Form("hSysTmp_cent%d",i),"",nbins_ana, ptbins_ana);
    makeHistTitle(h,"","Jet p_{T} (GeV/c)","Systematic uncertainty");
    h->SetAxisRange(-0.7,0.7,"Y");
    if(isPP==1) h->SetAxisRange(-0.3,0.5,"Y");
    h->SetAxisRange(draw_pt_low[i], draw_pt_high[i], "X");
    TH1F * sysJER;
    TH1F * sysJetID;
    TH1F * sysSmear;
    TH1F * sysBKG;
    TH1F * sysNoise;
    h->Draw();
    TH1F* sys = drawEnvelope(hSys[i],"same",hSys[i]->GetLineColor(),1001,hSys[i]->GetLineColor(),-1);
    TH1F* sysIter = drawEnvelope(hSysIter[i],"same",hSysIter[i]->GetLineColor(),3004,hSysIter[i]->GetLineColor(),-1);
    TH1F* sysPrior = drawEnvelope(hSysPrior[i],"same",hSysPrior[i]->GetLineColor(),3004,hSysPrior[i]->GetLineColor(),-1);
    TH1F* sysJEC = drawEnvelope(hSysJEC[i],"same",hSysJEC[i]->GetLineColor(),3005,hSysJEC[i]->GetLineColor(),-1);
    //for(int bin = 1; bin<=nbins_ana; ++bin)
    //hSysJER[i]->SetBinContent(bin, 1.046);
    sysJER = drawEnvelope(hSysJER[i],"same",hSysJER[i]->GetLineColor(),3003,hSysJER[i]->GetLineColor(),-1);
    if(i != 6) sysJetID = drawEnvelope(hSysJetID[i],"same",hSysJetID[i]->GetLineColor(),3005,hSysJetID[i]->GetLineColor(),-1);
    if(i != 6) sysSmear =  drawEnvelope(hSysSmear[i],"same",hSysSmear[i]->GetLineColor(),3001,hSysSmear[i]->GetLineColor(),-1);
    if(i != 6) sysBKG =  drawEnvelope(hSysBKG[i],"same",hSysBKG[i]->GetLineColor(),3001,hSysBKG[i]->GetLineColor(),-1);
    //TH1F* sysEff = drawEnvelope(hSysEff[i],"same",hSysEff[i]->GetLineColor(),3002,hSysEff[i]->GetLineColor(),-1);
    if(!isRAA) sysNoise = drawEnvelope(hSysNoise[i],"same",hSysNoise[i]->GetLineColor(),3001,hSysNoise[i]->GetLineColor(),-1);

    for(int k = 0; k<sys->GetNbinsX(); ++k){
      if(k < sys->FindBin(ptcut_low[i]) || k > sys->FindBin(ptcut_high[i])){
	sys->SetBinContent(k, 0);
	sysIter->SetBinContent(k, 0);
	sysPrior->SetBinContent(k, 0);
	sysJEC->SetBinContent(k, 0);
	sysJER->SetBinContent(k, 0);
	if(i != 6)sysJetID->SetBinContent(k, 0);
	if(i != 6)sysSmear->SetBinContent(k, 0);
	if(i != 6)sysBKG->SetBinContent(k, 0);
	//sysEff->SetBinContent(k, 0);
	if(!isRAA) sysNoise->SetBinContent(k, 0);
      }
    }
    
    TLine *l = new TLine(h->GetBinLowEdge(1),0,h->GetBinLowEdge(h->GetNbinsX()+1),0);
    //l->Draw();
    TLine *l2 = new TLine(h->GetBinLowEdge(1),-0.25,h->GetBinLowEdge(1),0.4);
    //l2->Draw();
    TLegend *leg;
    if(isPP == 0)leg = myLegend(0.22,0.6,0.65,0.93);
    if(isPP == 1)leg = myLegend(0.16,0.55,0.65,0.85);
    leg->SetTextSize(0.035);
    leg->AddEntry(sys,"Total Systematics","f");
    leg->AddEntry(sysIter,"Unfolding Iterations","f");
    leg->AddEntry(sysPrior,"Unfolding Prior Smearing","f");
    leg->AddEntry(sysJEC,"Jet Energy Scale","f");
    leg->AddEntry(sysJER,"Jet Energy Resolution","f");
    //leg->AddEntry(sysEff,"Jet Trigger Efficiency","f");
    if(i != 6)leg->AddEntry(sysJetID, "Jet ID efficiency","f");
    if(i != 6)leg->AddEntry(sysSmear,"UE Diff between Data and MC","f");
    if(i != 6)leg->AddEntry(sysBKG,"Combinatorial BKG Subtraction","f");
    if(!isRAA)leg->AddEntry(sysNoise,"HCAL Noise","f");
    if (i==nbins_cent-1 && isPP == 0)leg->Draw();
    else if(isPP == 1) leg->Draw();
  }	
 
};

void getRAA_Uncert(TH1F *hRAA, TH1F * hPbPb, TH1F * hPP){
  hRAA->Print("base");
  for(int bin = 1; bin<=hPbPb->GetNbinsX(); ++bin){
    double pbpb_uncert = 1.0 - hPbPb->GetBinContent(bin);
    double pp_uncert = 1.0 - hPP->GetBinContent(bin);

    // when they are fully un-correlated
    //double raa_uncert =  1.0 + TMath::Sqrt((pbpb_uncert)*(pbpb_uncert) + (pp_uncert)*(pp_uncert));

    // when the JER uncertainty cancels between the two for each centrality bin where peripheral PbPb ends up similar to PP
    double raa_uncert = 1.0 + TMath::Abs(pbpb_uncert - pp_uncert);

    // when the JER uncertainty cancels between the two for each centrality bin where peripheral PbPb ends up similar to PP
    //double raa_uncert = 0.0;
    //if(pp_uncert!=0) raa_uncert =  1.0 + TMath::Abs((float)pbpb_uncert/pp_uncert);
    //else raa_uncert = 1.0 + pbpb_uncert;
    
    hRAA->SetBinContent(bin, raa_uncert);
    cout<<"ptbin = "<<hRAA->GetBinCenter(bin)<<"; value = "<<hRAA->GetBinContent(bin)<<endl;
  }
}

// Remove error 
void removeError(TH1F *h)
{
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      h->SetBinError(i,0);
    }   	
}

void getTotalUncert(TH1F * hUncert, TH1F * hSys, TH1F * hSpectra, int mode = 0){
  if(mode == 0) hSpectra->Print("base");
  for(int bin = 1; bin<=hUncert->GetNbinsX(); ++bin){
    if(hSpectra->GetBinContent(bin) <= 0.0) continue;
    double sys = 0.0;
    if(mode == 0) sys = hSys->GetBinContent(bin)-1;
    double stat = (double)hSpectra->GetBinError(bin)/hSpectra->GetBinContent(bin);
    double uncert = (double)TMath::Sqrt(sys*sys + stat*stat);
    
    if(mode == 0) cout<<endl<<"Value = "<<hSpectra->GetBinContent(bin)<<endl;
    if(mode == 0) cout<<"Sys Error = "<<sys<<", stat error = "<<stat<<endl;
    if(mode == 0) cout<<"Total Uncert = "<<uncert<<endl;
   
    hUncert->SetBinContent(bin, 1.0+uncert);    
  }
}
  
// Remove Zero
void removeZero(TH1 *h)
{
  double min = 0;
  for(int i = 1;i<h->GetNbinsX();i++){
    if(h->GetBinContent(i)>min&&h->GetBinContent(i)>0)
      min = h->GetBinContent(i);
  }

  for(int i = 1;i<h->GetNbinsX();i++){
    if(h->GetBinContent(i) == 0){
      h->SetBinContent(i,min/10.);
      h->SetBinError(i,min/10.);
    }
  }
}


// make systematic histogram
void checkMaximumSys(TH1F *hSys, TH1F *h, int opt=0,double minVal = 1)
{
  if (h->GetNbinsX()!=hSys->GetNbinsX()) {
    cout <<"ERROR! Different NBins in subroutine checkMaximumSys!"<<endl;
  } else {
    double val = minVal;
    for (int i=1;i<=h->GetNbinsX();i++) {
      //cout <<i<<" "<<val<<" "<<hSys->GetBinContent(i)<<" "<<h->GetBinContent(i)<<endl;
      if (h->GetBinContent(i)==0) continue;
      if (opt==0) val=minVal;
      if (fabs(hSys->GetBinContent(i))>val) val = fabs(hSys->GetBinContent(i));
      if (fabs(h->GetBinContent(i)-1)+1>val) val=fabs(h->GetBinContent(i)-1)+1;
      hSys->SetBinContent(i,val);
    }
  }
}



void prepareNcollUnc(int nbins, float maxpt=299.){
	
  int fillsty = 1001;
	
  const int n = nbins;
	
  double xvalue[n];
  double yvalue[n];
  double xerror[n];
  double yerror1[n], yerror2[n], yerror3[n], yerror4[n], yerror5[n], yerror6[n];

  for(int i=0;i<n;i++){

    xvalue[i] = ptbins_ana[i];
    yvalue[i] = 1.0;
    xerror[i] = 0.0;  
    
    // TAA
    yerror1[i]=0.0409, yerror2[i]=0.0459, yerror3[i]=0.0578, yerror4[i]=0.0944, yerror5[i]=0.143, yerror6[i]=0.176;
    
    // add 6% error 
    yerror1[i]=TMath::Sqrt(yerror1[i]*yerror1[i]+0.06*0.06);
    yerror2[i]=TMath::Sqrt(yerror2[i]*yerror2[i]+0.06*0.06);
    yerror3[i]=TMath::Sqrt(yerror3[i]*yerror3[i]+0.06*0.06);
    yerror4[i]=TMath::Sqrt(yerror4[i]*yerror4[i]+0.06*0.06);
    yerror5[i]=TMath::Sqrt(yerror5[i]*yerror5[i]+0.06*0.06);
    yerror6[i]=TMath::Sqrt(yerror6[i]*yerror6[i]+0.06*0.06);
    //cout<<"TAA + Lumi uncertainty = "<<yerror1[i]<<endl;   
  
  }
  
  // int ci = 29;
  int ci = 15;

  tTAAerr[0] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror1);
  tTAAerr[0]->SetFillColor(ci);
  tTAAerr[0]->SetLineColor(ci);
  tTAAerr[0]->SetFillStyle(fillsty);

  tTAAerr[1] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror2);
  tTAAerr[1]->SetFillColor(ci);
  tTAAerr[1]->SetFillStyle(fillsty);

  tTAAerr[2] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror3);
  tTAAerr[2] ->SetFillColor(ci);
  tTAAerr[2] ->SetFillStyle(fillsty);

  tTAAerr[3]  = new TGraphErrors(n,xvalue,yvalue,xerror,yerror4);
  tTAAerr[3]->SetFillColor(ci);
  tTAAerr[3]->SetFillStyle(fillsty);

  tTAAerr[4] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror5);
  tTAAerr[4]->SetFillColor(ci);
  tTAAerr[4]->SetFillStyle(fillsty);

  tTAAerr[5] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror6);
  tTAAerr[5]->SetFillColor(ci);
  tTAAerr[5]->SetFillStyle(fillsty);
  


}

void DrawNpartTAABand(){
  double xvalueNpart[6];
  double yerrorNpart[6];

  xvalueNpart[0] = 381.29; xvalueNpart[1] = 329.41; xvalueNpart[2] = 224.28;
  xvalueNpart[3] = 108.12; xvalueNpart[4] = 42.04;  xvalueNpart[5] = 11.43;

  yerrorNpart[0]=0.0409, yerrorNpart[1]=0.0459, yerrorNpart[2]=0.0578, yerrorNpart[3]=0.0944, yerrorNpart[4]=0.143, yerrorNpart[5]=0.176;
  	
  int ci = 30;
  	
  for (int i=0;i<6;i++) {
		
    TBox *b = new TBox(xvalueNpart[i]-5,1.-yerrorNpart[i]/2,xvalueNpart[i]+5,1.+yerrorNpart[i]/2);
    b->SetFillColor(ci);
    b->SetFillStyle(3001);
    b->SetLineColor(ci);
    b->Draw();
  }
 
 
}

void drawNpartTAAlegend(double x1 = 0.2, double y1 = 0.2, double x2 = 0.22, double y2 = 0.22){

  int ci = 30;
  TBox *b = new TBox(x1, y1, x2, y2);
  b->SetFillColor(ci);
  b->SetFillStyle(3001);
  b->SetLineColor(ci);
  b->Draw();
  
}

void drawSysBoxlegend(double x1 = 0.2, double y1 = 0.2, double x2 = 0.22, double y2 = 0.22, int color = 2, int style = 0){

  TBox *b = new TBox(x1, y1, x2, y2);
  b->SetFillColor(color);
  b->SetFillStyle(style);
  b->SetLineColor(color);
  b->Draw(); 

  TLine * lineup = new TLine(x1, y1, x2, y1);
  lineup->SetLineColor(color);
  lineup->SetLineStyle(style);
  /* lineup->Draw(); */
  TLine * linedown = new TLine(x1, y2, x2, y2);
  linedown->SetLineColor(color);
  linedown->SetLineStyle(style);
  /* linedown->Draw(); */

}

void drawSysBoxNpartlegend(double x1 = 0.2, double y1 = 0.2, double x2 = 0.22, double y2 = 0.22, int style = 1001, int color = kGray){

  TBox *b = new TBox(x1, y1, x2, y2);
  b->SetFillColor(color);
  b->SetFillStyle(style);
  b->SetLineColor(color);
  /* b->Draw(); */
  
  TLine * lineup = new TLine(x1, y1, x2, y1);
  lineup->SetLineColor(color);
  lineup->SetLineStyle(style);
  lineup->Draw();
  TLine * linedown = new TLine(x1, y2, x2, y2);
  linedown->SetLineColor(color);
  linedown->SetLineStyle(style);
  linedown->Draw();
  
}

void DrawTAALumiUncert(int bin){

  double xvalue, yerrorNpart[6];

  xvalue = 294.0;

  yerrorNpart[0]=0.0409, yerrorNpart[1]=0.0459, yerrorNpart[2]=0.0578, yerrorNpart[3]=0.0944, yerrorNpart[4]=0.143, yerrorNpart[5]=0.176;
  	  
  // add 6% error 
  yerrorNpart[0]=TMath::Sqrt(yerrorNpart[0]*yerrorNpart[0]+0.037*0.037);
  yerrorNpart[1]=TMath::Sqrt(yerrorNpart[1]*yerrorNpart[1]+0.037*0.037);
  yerrorNpart[2]=TMath::Sqrt(yerrorNpart[2]*yerrorNpart[2]+0.037*0.037);
  yerrorNpart[3]=TMath::Sqrt(yerrorNpart[3]*yerrorNpart[3]+0.037*0.037);
  yerrorNpart[4]=TMath::Sqrt(yerrorNpart[4]*yerrorNpart[4]+0.037*0.037);
  yerrorNpart[5]=TMath::Sqrt(yerrorNpart[5]*yerrorNpart[5]+0.037*0.037);

  int ci = 30;

  TBox * b = new TBox(xvalue-5, 1.-yerrorNpart[bin]/2, xvalue+5, 1.+yerrorNpart[bin]/2);
  b->SetFillColor(ci);
  b->SetFillStyle(3001);
  b->SetLineColor(ci);
  b->Draw();
  
}

 

void dumpDatatoTxt(const char *centbin,TH1F *h, TH1F *hsys, TH1F *htotStat, const char *txtfile)
{
  ofstream outf(txtfile,ios::out);
  for(int ix=1;ix<=h->GetNbinsX();ix++){
    double pt = h->GetBinCenter(ix);
    double val = h->GetBinContent(ix);
    double Uncorerr = h->GetBinError(ix);
    double syserr = hsys->GetBinContent(ix)-1;
    double totStaterr = htotStat->GetBinError(ix);

    outf<<setprecision(0) << fixed <<pt<<"\t" << setprecision(3) << fixed <<val<<"\t" << setprecision(5) << fixed << totStaterr<<"\t"<< setprecision(4) << fixed << syserr*val << endl;
  }
  outf.close();
}
