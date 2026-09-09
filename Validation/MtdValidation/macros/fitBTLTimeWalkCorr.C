// This ROOT macro determines the BTL time-walk corrections fitting the TProfile DeltaTsimvsE in
// the BTL/LocalReco validation with a power-law + constant function: f(x) = [0]*x^[1] + [2]
//
// CMSSW instructions to calculate new time-walk corrections for each aging scenario:
//  - run the step 3 switching off the current time-walk corrections:
//
//       process.mtdUncalibratedRecHits.barrel.timeWalkCorrection = cms.string('0.')
//
//  - run the MTD validation enabling the flag:
//
//       process.btlLocalRecoValid.FillTimeWalkPlots = True
//
// Macro usage:
//   root -l 'fitBTLTimeWalkCorr.C+("DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root")'

#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include <cstdio>

TProfile* p_deltaT_vs_E = nullptr;
TF1* fitFunc = nullptr;

TCanvas* c1 = nullptr;

void fitBTLTimeWalkCorr(const char* fileName = "DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root")
{
    const char* histPath = "DQMData/Run 1/MTD/Run summary/BTL/LocalReco/DeltaTsimvsE";

    TFile* f = TFile::Open(fileName);
    if (!f || f->IsZombie()) {
        printf("Error: could not open file %s\n", fileName);
        return;
    }

    f->GetObject(histPath, p_deltaT_vs_E);

    if (!p_deltaT_vs_E) {
        printf("Error: could not find TProfile 'DeltaTsimvsE' in %s\n", fileName);
        f->Close();
        return;
    }

    // Determine the fit range from the first and last bins with non-zero content
    int firstBin = -1, lastBin = -1;
    for (int i = 1; i <= p_deltaT_vs_E->GetNbinsX(); ++i) {
        if (p_deltaT_vs_E->GetBinContent(i) != 0) {
            if (firstBin < 0) firstBin = i;
            lastBin = i;
        }
    }

    if (firstBin < 0) {
        printf("Error: all bins of 'DeltaTsimvsE' are empty\n");
        f->Close();
        return;
    }

    double xMin = p_deltaT_vs_E->GetXaxis()->GetBinLowEdge(firstBin);
    double xMax = p_deltaT_vs_E->GetXaxis()->GetBinUpEdge(lastBin);

    // Define fit function: power law + constant
    fitFunc = new TF1("fitFunc", "[0]*pow(x,[1])+[2]", xMin, xMax);
    fitFunc->SetParameters(1.0, -1.0, 0.0);
    fitFunc->SetParNames("Norm", "Power", "Const");

    // "R" restricts the fit to the function range; bins with zero content
    // (no entries) are automatically excluded from the chi2 calculation
    p_deltaT_vs_E->Fit(fitFunc, "R");

    c1 = new TCanvas("c1", "DeltaTsimvsE fit", 800, 600);
    p_deltaT_vs_E->Draw();
    fitFunc->Draw("same");

    printf("Fit results:\n");
    printf("  Norm  = %.4g +/- %.4g\n", fitFunc->GetParameter(0), fitFunc->GetParError(0));
    printf("  Power = %.4g +/- %.4g\n", fitFunc->GetParameter(1), fitFunc->GetParError(1));
    printf("  Const = %.4g +/- %.4g\n", fitFunc->GetParameter(2), fitFunc->GetParError(2));
    printf("  Chi2/NDF = %.4g / %d = %.4g\n",
           fitFunc->GetChisquare(), fitFunc->GetNDF(),
           fitFunc->GetChisquare() / fitFunc->GetNDF());
}
