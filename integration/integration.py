from __future__ import division
import lhapdf
import math
import random
import argparse
from subprocess import call

import sys
import os

def parsing():
    parser = argparse.ArgumentParser(description='Parameters Passed to the STRINGS Generator')
    parser.add_argument('-v', '--version', action='version', version='STRINGS-2.00')
    parser.add_argument('-RandGenSEED', '--RandGenSEEDvalue', default=123456, help='Seed for Random Number Generator',type=float)
    parser.add_argument('-Number', '--Numbervalue', default=10000, help='Number of Events',type=int)
    parser.add_argument('-Ms', '--Msvalue', default=7000, help='String Scale of the String Theory (GeV)',type=float)
    parser.add_argument('-COME', '--COMEvalue', default=13000, help='Centre of Mass Energy (GeV)',type=float)
    parser.add_argument('-MinMass', '--MinMassvalue', default=6000, help='Lower Bound of the Invariant Mass (GeV)',type=float)
    parser.add_argument('-MaxMass', '--MaxMassvalue', default=8000, help='Upper Bound of the Invariant Mass (GeV)',type=float)
    parser.add_argument('-yMax', '--yMaxvalue', default=2.5, help='Upper Limit for the Rapidity of the Outgoing Partons',type=float)
    parser.add_argument('-PDFSet', '--PDFSetvalue', default="cteq6l1", help='Parton Distribution Function Set of the LHAPDF',type=str)
    parser.add_argument('-PDFScale', '--PDFScalevalue', default=-22, help='Scale at which the Parton Distribution Function is Evaluated (GeV)',type=float)
    parser.add_argument('-Coupling', '--Couplingvalue', default=-1, help='Coupling Constant',type=float)
    parser.add_argument('-CouplingScale', '--CouplingScalevalue', default=-22, help='Scale at Which the Running Coupling Constant is Calculated(GeV)',type=float)
    parser.add_argument('-FirstStringCoeff', '--FirstStringCoeffvalue', default=False, type=lambda x: (str(x).lower() == 'true'), help='First String Resonance')
    parser.add_argument('-SecondStringCoeff', '--SecondStringCoeffvalue', default=False, type=lambda x: (str(x).lower() == 'true'), help='Second String Resonance')
    parser.add_argument('-QCDCoeff', '--QCDCoeffvalue', default=False, type=lambda x: (str(x).lower() == 'true'), help='QCD Tree-Level')
    parser.add_argument('-dMass', '--dMassvalue', default=5e-3, help='Mass of the Down Quark (GeV)',type=float)
    parser.add_argument('-uMass', '--uMassvalue', default=2e-3, help='Mass of the Up Quark (GeV)',type=float)
    parser.add_argument('-sMass', '--sMassvalue', default=1e-3, help='Mass of the Strange Quark (GeV)',type=float)
    parser.add_argument('-cMass', '--cMassvalue', default=1.25, help='Mass of the Charm Quark (GeV)',type=float)
    parser.add_argument('-bMass', '--bMassvalue', default=4.2, help='Mass of the Bottom Quark (GeV)',type=float)
    parser.add_argument('-tMass', '--tMassvalue', default=172.5, help='Mass of the Top Quark (GeV)',type=float)
    parser.add_argument('-gg2gg', '--gg2ggvalue', default=False, type=lambda x: (str(x).lower() == 'true'), help='gg -> gg Subprocess')
    parser.add_argument('-gg2qqbar', '--gg2qqbarvalue',default=False, type=lambda x: (str(x).lower() == 'true'), help='gg -> qqbar Subprocess')
    parser.add_argument('-gq2gq', '--gq2gqvalue', default=False, type=lambda x: (str(x).lower() == 'true'), help='gq -> gq Subprocess')
    parser.add_argument('-gqbar2gqbar', '--gqbar2gqbarvalue', default=False, type=lambda x: (str(x).lower() == 'true'), help='gqbar -> gqbar Subprocess')
    parser.add_argument('-qqbar2gg', '--qqbar2ggvalue',default=False, type=lambda x: (str(x).lower() == 'true'), help='qqbar -> gg Subprocess')
    parser.add_argument('-gg2gGamma', '--gg2gGammavalue',default=False, type=lambda x: (str(x).lower() == 'true'), help='gg -> gGamma Subprocess')
    parser.add_argument('-gq2qGamma', '--gq2qGammavalue', default=False, type=lambda x: (str(x).lower() == 'true'), help='gq -> qGamma Subprocess')
    
    
    
    parser.add_argument('-lowerBound', '--lowerBoundvalue', default=1000, help='Lower Bound of Mass',type=int)
    parser.add_argument('-upperBound', '--upperBoundvalue', default=9000, help='Upper Bound of Mass',type=int)
    parser.add_argument('-stepSize', '--stepSizevalue', default=20, help='Step size',type=int)

    


    args = parser.parse_args()
    
    return args

def checker():
    if(args.Numbervalue <= 0):
        print("Error: Invalid Value for the Number of Generated Events\nIt Should be Positive\n")
        return False
    if(args.Msvalue <= 0):
        print("Error: Invalid Value for the String Scale\nIt Should be Positive\n")
        return False
    if(args.COMEvalue <= 0):
        print("Error: Invalid Value for Centre of Mass Energy\nIt Should be Positive\n")
        return False
    if(args.MinMassvalue <= 0):
        print("Error: Invalid Value for the Lower Bound of the Invariant Mass\nIt Should be Positive\n")
        return False
    if(args.MaxMassvalue <= 0):
        print("Error: Invalid Value for the Upper Bound of the Invariant Mass\nIt Should be Positive\n")
        return False
    if(args.MaxMassvalue > args.COMEvalue):
        print("Error: Invalid Value for the Upper Bound of the Invariant Mass\nIt Should Not be larger than the Centre of Mass Energy\n")
        return False
    elif(args.MaxMassvalue <= args.MinMassvalue):
        print("Error: Invalid Value for the Upper Bound of the Invariant Mass\nIt Should be greater than the Lower Bound of the Invariant Mass\n")
        return False
    elif(args.MaxMassvalue > args.COMEvalue):
        print("Error: Invalid Value for the Upper Bound of the Invariant Mass\nIt Should not be greater than the Centre of Mass Energy\n")
        return False
    if(args.yMaxvalue == -1):
        args.yMaxvalue = float("inf")
    elif(args.yMaxvalue < 0):
        print("Error: Invalid Value for the Upper Bound of the Rapidity\nIt Should be Positive or for the Upper Bound to be Infinity, Use -1\n")
        return False
    if(args.PDFScalevalue <= 0 and args.PDFScalevalue != -22):
        print("Error: Invalid Value for the Scale at Which PDFs are Determined\nIt Should be Positive\n")
        return False
    if(args.Couplingvalue == -1):
        args.Couplingvalue = -1
    elif(args.Couplingvalue < 0):
        print("Error: Invalid Value for the Coupling Constant\nIt Should be positive or for running coupling use -1\n")
        return False
    if(args.CouplingScalevalue <= 0 and args.CouplingScalevalue != -22):
        print("Error: Invalid Value for the Scale at Which Running Coupling is Measured\nIt Should be Positive\n")
        return False
    if(args.QCDCoeffvalue == False):
        if(args.FirstStringCoeffvalue == False):
            if(args.SecondStringCoeffvalue == False):
                return False
                print("Error: Type of Event Generation is not Specified\n")
    if(args.gg2ggvalue == False):
        if(args.gg2qqbarvalue == False):
            if(args.gq2gqvalue == False):
                if(args.gqbar2gqbarvalue == False):
                    if(args.qqbar2ggvalue == False):
                        if(args.gg2gGammavalue == False):
                            if(args.gq2qGammavalue == False):
                                return False
                                print("Error: No Subprocess is Turned On\n")
    if(args.gg2ggvalue == False):
        if(args.gg2qqbarvalue == False):
            if(args.gq2gqvalue == False):
                if(args.gqbar2gqbarvalue == False):
                    if(args.qqbar2ggvalue == False):
                        if(args.gg2gGammavalue == True):
                            if(args.QCDCoeffvalue == True):
                                return False
                                print("Error: QCD Diparton Cannot be Produced Through gg->gGamma Scattering\nPlease Turn On Either of the 2->2 Partonic Scatterings")
                            if(args.SecondStringCoeffvalue == True):
                                return False
                                print("Error: Second string cannot be produced through gg->gGamma scattering\nPlease Turn On Either of the 2->2 Partonic Scatterings")
                        if(args.gq2qGammavalue == True):
                            if(args.QCDCoeffvalue == True):
                                return False
                                print("Error: QCD Diparton Cannot be Produced Through gq->qGamma Scattering\nPlease Turn On Either of the 2->2 Partonic Scatterings")
                            if(args.SecondStringCoeffvalue == True):
                                return False
                                print("Error: Second string cannot be produced through gq->qGamma scattering\nPlease Turn On Either of the 2->2 Partonic Scatterings")
    if(args.dMassvalue < 0):
        print("Error: Invalid Value for the Mass of the Down Quark\nIt Should be Positive\n")
        return False
    if(args.uMassvalue < 0):
        print("Error: Invalid Value for the Mass of the UP Quark\nIt Should be Positive\n")
        return False
    if(args.sMassvalue < 0):
        print("Error: Invalid Value for the Mass of the Strange Quark\nIt Should be Positive\n")
        return False
    if(args.cMassvalue < 0):
        print("Error: Invalid Value for the Mass of the Charm Quark\nIt Should be Positive\n")
        return False
    if(args.bMassvalue < 0):
        print("Error: Invalid Value for the Mass of the Bottom Quark\nIt Should be Positive\n")
        return False
    if(args.tMassvalue < 0):
        print("Error: Invalid Value for the Mass of the Top Quark\nIt Should be Positive\n")
        return False
    if(args.gg2ggvalue != True and args.gg2ggvalue != False):
        print("Error: Invalid Value; In order to enable or disable gg --> gg subprocess, set gg2gg to true or false, respectively.\n")
        return False
    if(args.gg2qqbarvalue != True and args.gg2qqbarvalue != False):
        print("Error: Invalid Value; In order to enable or disable gg --> qqbar subprocess, set gg2qqbar to true or false, respectively.\n")
        return False
    if(args.gq2gqvalue != True and args.gq2gqvalue != False):
        print("Error: Invalid Value; In order to enable or disable gq --> gq subprocess, set gq2gq to true or false, respectively.\n")
        return False
    if(args.gqbar2gqbarvalue != True and args.gqbar2gqbarvalue != False):
        print("Error: Invalid Value; In order to enable or disable gqbar --> gqbar subprocess, set gqbar2gqbar to true or false, respectively.\n")
        return False
    if(args.qqbar2ggvalue != True and args.qqbar2ggvalue != False):
        print("Error: Invalid Value; In order to enable or disable qqbar --> gg subprocess, set qqbar2gg to true or false, respectively.\n")
        return False
    if(args.gg2gGammavalue != True and args.gg2gGammavalue != False):
        print("Error: Invalid Value; In order to enable or disable gg --> gGamma subprocess, set gg2gGamma to true or false, respectively.\n")
        return False
    if(args.gq2qGammavalue != True and args.gq2qGammavalue != False):
        print("Error: Invalid Value; In order to enable or disable gq --> qGamma subprocess, set gq2qGamma to true or false, respectively.\n")
        return False
    
    return True

def main(RandGenSEED, Number, Ms, COME, MinMass, MaxMass, yMax, PDFSet, PDFScale, Coupling, CouplingScale, StringCoeff, SecondStringCoeff, QCDCoeff, dMass, uMass, sMass, cMass, bMass, tMass, gg2gg, gg2qqbar, gq2gq, gqbar2gqbar, qqbar2gg, gg2gGamma, gq2qGamma, mass):
    if(PDFScale == -22):
        PDFScale = Ms
        
    if(Coupling == -1):
        if(CouplingScale == -22):
            CouplingScale = Ms
            
    else:
        pass
        
    if(StringCoeff):
        StringCoeff11 = 1
        
    else:
        StringCoeff11 = 0
        
    if(SecondStringCoeff):
        SecondStringCoeff11 = 1
       
    else:
        SecondStringCoeff11 = 0
        
    if(QCDCoeff):
        QCDCoeff11 = 1
        
    else:
        QCDCoeff11 = 0

    NPRUP = 0

    if(gg2gg):
        gg2gg11 = 1
        NPRUP = NPRUP + 1
        
    else:
        gg2gg11 = 0

    if(gg2qqbar):
        gg2qqbar11 = 1
        NPRUP = NPRUP + 1
        
    else:
        gg2qqbar11 = 0

    if(gq2gq):
        gq2gq11 = 1
        NPRUP = NPRUP + 1
        
    else:
        gq2gq11 = 0

    if(gqbar2gqbar):
        gqbar2gqbar11 = 1
        NPRUP = NPRUP + 1
        
    else:
        gqbar2gqbar11 = 0

    if(qqbar2gg):
        qqbar2gg11 = 1
        NPRUP = NPRUP + 1
        
    else:
        qqbar2gg11 = 0

    if(gg2gGamma):
        gg2gGamma11 = 1
        NPRUP = NPRUP + 1
        
    else:
        gg2gGamma11 = 0

    if(gq2qGamma):
        gq2qGamma11 = 1
        NPRUP = NPRUP + 1
        
    else:
        gq2qGamma11 = 0


    DistFunc = lhapdf.mkPDF(PDFSet, 0)

    # Quark Masses
    QuarkMasses = [dMass, uMass, sMass, cMass, bMass, tMass]

    # Centre of mass energy
    s = COME**2

    # Number of colours and flavours
    N = 3
    Nf = 6

    # Coupling constants
    if(Coupling == -1):
        # Running Coupling Constant, calculated at the CouplingScale
        alpha3 = 1/(1/0.118 + 7/2/math.pi*math.log(CouplingScale/91.2))
    else:
        # Fixed coupling constant
        alpha3 = Coupling

    g3 = math.sqrt(4*math.pi*alpha3) # Strong coupling constant
    QED = 1/128  # QED coupling constant
    Q = 1/6      # Fundamental charge

    # Decay widths that will be used in the scattering amplitudes
    # Decay Widths of First Resonances
    GJ0gS   = g3**2/4/math.pi*Ms*N/4
    GJ0cS   = g3**2/4/math.pi*Ms*N/2
    GJhfqS  = g3**2/4/math.pi*Ms*N/8
    GJ2gS   = g3**2/4/math.pi*Ms*(N/10 + Nf/40)
    GJ2cS   = g3**2/4/math.pi*Ms*(N/5  + Nf/40)
    GJ3hfqS = g3**2/4/math.pi*Ms*N/16

    # Decay Widths of Second Resonances
    GssJhfqS  = g3**2/4/math.pi*math.sqrt(2)*Ms*N/24
    GssJ3hfqS = g3**2/4/math.pi*math.sqrt(2)*Ms*19*N/240
    GssJ5hfqS = g3**2/4/math.pi*math.sqrt(2)*Ms*N/60



    def MonteCarlogggg(parameters, *args):

        Y, y                           = parameters
        M, FstpartonID, ScndpartonID   = args
        
        tau = M**2/s

        xa  = math.sqrt(tau)*math.exp( Y)
        xb  = math.sqrt(tau)*math.exp(-Y)

        shat = M**2
        that = -0.5*M**2*math.exp(-y)/math.cosh(y)
        uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

        # Production of the first string resonance
        MSquareString = StringCoeff11*(8/N**2*g3**4/Ms**4*((N**2-4)**2/4/(N**2-1)*(Ms**8/((shat-Ms**2)**2+(Ms*GJ0gS)**2) + (uhat**4+that**4)/((shat-Ms**2)**2 + (Ms*GJ2gS)**2)) + Ms**8/((shat-Ms**2)**2+(Ms*GJ0cS)**2) + (uhat**4+that**4)/((shat-Ms**2)**2 + (Ms*GJ2cS)**2)))

        # Production of the QCD diparton
        MSquareQCD    = QCDCoeff11*(g3**4*(1/shat**2+1/that**2+1/uhat**2)*(2*N**2/(N**2-1)*(shat**2+that**2+uhat**2)+4*(3-N**2)/N**2/(N**2-1)*(shat+that+uhat)**2))

        # Total scattering amplitude
        MSquareTotal  = MSquareString + MSquareQCD

        #Convolution of the scattering amplitudes with parton distribution function
        return(gg2gg11*DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale)*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))


    def MonteCarloggqqbar(parameters, *args):
        
        Y, y                           = parameters
        M, FstpartonID, ScndpartonID   = args
        
        tau = M**2/s

        xa  = math.sqrt(tau)*math.exp( Y)
        xb  = math.sqrt(tau)*math.exp(-Y)

        shat = M**2
        that = -0.5*M**2*math.exp(-y)/math.cosh(y)
        uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

        # Production of the first string resonance
        MSquareString = StringCoeff11*(2/N/(N**2-1)*Nf*g3**4/Ms**4*((N**2-4)/2*uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2gS)**2) + uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2cS)**2)))

        # Production of the QCD diparton
        MSquareQCD    = QCDCoeff11*(2*g3**4*Nf*(that**2+uhat**2)/shat**2*(1/2/N/uhat/that*(that+uhat)**2-N/(N**2-1)))

        # Total scattering amplitude
        MSquareTotal  = MSquareString + MSquareQCD

        #Convolution of the scattering amplitudes with parton distribution function
        return(gg2qqbar11*DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale)*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

    #gq->gq
    def MonteCarlogqgq(parameters, *args):
        
        Y, y                           = parameters
        M, FstpartonID, ScndpartonID   = args
        
        tau = M**2/s

        xa  = math.sqrt(tau)*math.exp( Y)
        xb  = math.sqrt(tau)*math.exp(-Y)

        shat = M**2
        that = -0.5*M**2*math.exp(-y)/math.cosh(y)
        uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

        # Production of the first string resonance
        MSquareString       = StringCoeff11*(0.5*(N**2-1)/2/N**2*g3**4/Ms**2*(-Ms**4*that/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - that**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2) - Ms**4*uhat/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - uhat**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2)))

        # Production of the QCD diparton
        MSquareQCD          = QCDCoeff11*(g3**4*(shat**2+uhat**2)/that**2*(1-(N**2-1)/2/N**2/uhat/shat*(shat+uhat)**2) + g3**4*(shat**2+that**2)/uhat**2*(1-(N**2-1)/2/N**2/that/shat*(shat+that)**2))

        # Production of the second string resonance
        MSquareSecondString = SecondStringCoeff11*(2*(N**2-1)/N**2*(g3**4/2/Ms**2*(1/9*Ms**4*(-uhat)/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJhfqS)**2) + (-1/9)*uhat*(3*that+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2)) + g3**4/8/Ms**6*(9/25*Ms**4*(-uhat)**3/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2) + (-1/25)*uhat**3*(5*that+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ5hfqS)**2))))

        #Total scattering amplitude
        MSquareTotal        = MSquareString + MSquareQCD + MSquareSecondString

        #Convolution of the scattering amplitudes with parton distribution function
        return(gq2gq11*(DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale) + DistFunc.xfxQ(FstpartonID,xb,PDFScale)*DistFunc.xfxQ(ScndpartonID,xa,PDFScale))*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

    #gqbar->gqbar
    def MonteCarlogqbargqbar(parameters, *args):
        
        Y, y                           = parameters
        M, FstpartonID, ScndpartonID   = args
        
        tau = M**2/s

        xa  = math.sqrt(tau)*math.exp( Y)
        xb  = math.sqrt(tau)*math.exp(-Y)

        shat = M**2
        that = -0.5*M**2*math.exp(-y)/math.cosh(y)
        uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

        # Production of the first string resonance
        MSquareString       = StringCoeff11*(0.5*(N**2-1)/2/N**2*g3**4/Ms**2*(-Ms**4*that/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - that**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2) - Ms**4*uhat/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - uhat**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2)))

        # Production of the QCD diparton
        MSquareQCD          = QCDCoeff11*(g3**4*(shat**2+uhat**2)/that**2*(1-(N**2-1)/2/N**2/uhat/shat*(shat+uhat)**2) + g3**4*(shat**2+that**2)/uhat**2*(1-(N**2-1)/2/N**2/that/shat*(shat+that)**2))

        # Production of the second string resonance
        MSquareSecondString = SecondStringCoeff11*(2*(N**2-1)/N**2*(g3**4/2/Ms**2*(1/9*Ms**4*(-that)/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJhfqS)**2) + (-1/9)*that*(3*uhat+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2)) + g3**4/8/Ms**6*(9/25*Ms**4*(-that)**3/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2) + (-1/25)*that**3*(5*uhat+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ5hfqS)**2))))

        # Total scattering amplitude
        MSquareTotal        = MSquareString + MSquareQCD + MSquareSecondString

        
        #Convolution of the scattering amplitudes with parton distribution function
        return(gqbar2gqbar11*(DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale) + DistFunc.xfxQ(FstpartonID,xb,PDFScale)*DistFunc.xfxQ(ScndpartonID,xa,PDFScale))*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

    #qqbar->gg
    def MonteCarloqqbargg(parameters, *args):

        Y, y                           = parameters
        M, FstpartonID, ScndpartonID   = args
        
        tau = M**2/s
        
        xa  = math.sqrt(tau)*math.exp( Y)
        xb  = math.sqrt(tau)*math.exp(-Y)

        shat = M**2
        that = -0.5*M**2*math.exp(-y)/math.cosh(y)
        uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

        # Production of the first string resonance
        MSquareString = StringCoeff11*(2*(N**2-1)/N**3*g3**4/Ms**4*((N**2-4)/2*uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2gS)**2) + uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2cS)**2)))

        # Production of the QCD diparton
        MSquareQCD    = QCDCoeff11*(2*g3**4*(that**2+uhat**2)/shat**2*((N**2-1)**2/2/N**3/uhat/that*(that+uhat)**2-(N**2-1)/N))

        # Total scattering amplitude
        MSquareTotal  = MSquareString + MSquareQCD

        #Convolution of the scattering amplitudes with parton distribution function
        return(qqbar2gg11*(DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale) + DistFunc.xfxQ(FstpartonID,xb,PDFScale)*DistFunc.xfxQ(ScndpartonID,xa,PDFScale))*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

    #gg->gGamma
    def MonteCarlogggGamma(parameters, *args):

        Y, y                           = parameters
        M, FstpartonID, ScndpartonID   = args
        
        tau = M**2/s
        
        xa  = math.sqrt(tau)*math.exp( Y)
        xb  = math.sqrt(tau)*math.exp(-Y)

        shat = M**2
        that = -0.5*M**2*math.exp(-y)/math.cosh(y)
        uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

        # Production of the first string resonance
        MSquareTotal  = StringCoeff11*5*g3**4*Q**2/3/Ms**4*(Ms**8/((shat-Ms**2)**2 + (GJ0gS*Ms)**2) + (that**4 + uhat**4)/((shat-Ms**2)**2 + (GJ2gS*Ms)**2))

        #Convolution of the scattering amplitudes with parton distribution function
        return(2*gg2gGamma11*DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale)*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

    #gq->qGamma
    def MonteCarlogqqGamma(parameters, *args):

        Y, y                           = parameters
        M, FstpartonID, ScndpartonID   = args
        
        tau = M**2/s
        
        xa  = math.sqrt(tau)*math.exp( Y)
        xb  = math.sqrt(tau)*math.exp(-Y)

        shat = M**2
        that = -0.5*M**2*math.exp(-y)/math.cosh(y)
        uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

        # Production of the first string resonance
        MSquareTotal  = -StringCoeff11*g3**4*Q**2/3/Ms**2*(Ms**4*uhat/((shat-Ms**2)**2 + (GJhfqS*Ms)**2) + uhat**3/((shat-Ms**2)**2 + (GJ3hfqS*Ms)**2))

        #Convolution of the scattering amplitudes with parton distribution function
        return(gq2qGamma11*(DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale) + DistFunc.xfxQ(FstpartonID,xb,PDFScale)*DistFunc.xfxQ(ScndpartonID,xa,PDFScale))*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

        def MonteCarloggqqbar(parameters, *args):
            
            Y, y                           = parameters
            M, FstpartonID, ScndpartonID   = args
            
            tau = M**2/s

            xa  = math.sqrt(tau)*math.exp( Y)
            xb  = math.sqrt(tau)*math.exp(-Y)

            shat = M**2
            that = -0.5*M**2*math.exp(-y)/math.cosh(y)
            uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

            # Production of the first string resonance
            MSquareString = StringCoeff11*(2/N/(N**2-1)*Nf*g3**4/Ms**4*((N**2-4)/2*uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2gS)**2) + uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2cS)**2)))

            # Production of the QCD diparton
            MSquareQCD    = QCDCoeff11*(2*g3**4*Nf*(that**2+uhat**2)/shat**2*(1/2/N/uhat/that*(that+uhat)**2-N/(N**2-1)))

            # Total scattering amplitude
            MSquareTotal  = MSquareString + MSquareQCD

            #Convolution of the scattering amplitudes with parton distribution function
            return(gg2qqbar11*DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale)*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

        #gq->gq
        def MonteCarlogqgq(parameters, *args):
            
            Y, y                           = parameters
            M, FstpartonID, ScndpartonID   = args
            
            tau = M**2/s

            xa  = math.sqrt(tau)*math.exp( Y)
            xb  = math.sqrt(tau)*math.exp(-Y)

            shat = M**2
            that = -0.5*M**2*math.exp(-y)/math.cosh(y)
            uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

            # Production of the first string resonance
            MSquareString       = StringCoeff11*(0.5*(N**2-1)/2/N**2*g3**4/Ms**2*(-Ms**4*that/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - that**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2) - Ms**4*uhat/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - uhat**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2)))

            # Production of the QCD diparton
            MSquareQCD          = QCDCoeff11*(g3**4*(shat**2+uhat**2)/that**2*(1-(N**2-1)/2/N**2/uhat/shat*(shat+uhat)**2) + g3**4*(shat**2+that**2)/uhat**2*(1-(N**2-1)/2/N**2/that/shat*(shat+that)**2))

            # Production of the second string resonance
            MSquareSecondString = SecondStringCoeff11*(2*(N**2-1)/N**2*(g3**4/2/Ms**2*(1/9*Ms**4*(-uhat)/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJhfqS)**2) + (-1/9)*uhat*(3*that+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2)) + g3**4/8/Ms**6*(9/25*Ms**4*(-uhat)**3/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2) + (-1/25)*uhat**3*(5*that+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ5hfqS)**2))))

            #Total scattering amplitude
            MSquareTotal        = MSquareString + MSquareQCD + MSquareSecondString

            #Convolution of the scattering amplitudes with parton distribution function
            return(gq2gq11*(DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale) + DistFunc.xfxQ(FstpartonID,xb,PDFScale)*DistFunc.xfxQ(ScndpartonID,xa,PDFScale))*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

        #gqbar->gqbar
        def MonteCarlogqbargqbar(parameters, *args):
            
            Y, y                           = parameters
            M, FstpartonID, ScndpartonID   = args
            
            tau = M**2/s

            xa  = math.sqrt(tau)*math.exp( Y)
            xb  = math.sqrt(tau)*math.exp(-Y)

            shat = M**2
            that = -0.5*M**2*math.exp(-y)/math.cosh(y)
            uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

            # Production of the first string resonance
            MSquareString       = StringCoeff11*(0.5*(N**2-1)/2/N**2*g3**4/Ms**2*(-Ms**4*that/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - that**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2) - Ms**4*uhat/((shat-Ms**2)**2+(Ms*GJhfqS)**2) - uhat**3/((shat-Ms**2)**2+(Ms*GJ3hfqS)**2)))

            # Production of the QCD diparton
            MSquareQCD          = QCDCoeff11*(g3**4*(shat**2+uhat**2)/that**2*(1-(N**2-1)/2/N**2/uhat/shat*(shat+uhat)**2) + g3**4*(shat**2+that**2)/uhat**2*(1-(N**2-1)/2/N**2/that/shat*(shat+that)**2))

            # Production of the second string resonance
            MSquareSecondString = SecondStringCoeff11*(2*(N**2-1)/N**2*(g3**4/2/Ms**2*(1/9*Ms**4*(-that)/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJhfqS)**2) + (-1/9)*that*(3*uhat+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2)) + g3**4/8/Ms**6*(9/25*Ms**4*(-that)**3/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ3hfqS)**2) + (-1/25)*that**3*(5*uhat+shat)**2/((shat-2*Ms**2)**2 + (math.sqrt(2)*Ms*GssJ5hfqS)**2))))

            # Total scattering amplitude
            MSquareTotal        = MSquareString + MSquareQCD + MSquareSecondString

            
            #Convolution of the scattering amplitudes with parton distribution function
            return(gqbar2gqbar11*(DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale) + DistFunc.xfxQ(FstpartonID,xb,PDFScale)*DistFunc.xfxQ(ScndpartonID,xa,PDFScale))*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

        #qqbar->gg
        def MonteCarloqqbargg(parameters, *args):

            Y, y                           = parameters
            M, FstpartonID, ScndpartonID   = args
            
            tau = M**2/s
            
            xa  = math.sqrt(tau)*math.exp( Y)
            xb  = math.sqrt(tau)*math.exp(-Y)

            shat = M**2
            that = -0.5*M**2*math.exp(-y)/math.cosh(y)
            uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

            # Production of the first string resonance
            MSquareString = StringCoeff11*(2*(N**2-1)/N**3*g3**4/Ms**4*((N**2-4)/2*uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2gS)**2) + uhat*that*(uhat**2+that**2)/((shat-Ms**2)**2+(Ms*GJ2cS)**2)))

            # Production of the QCD diparton
            MSquareQCD    = QCDCoeff11*(2*g3**4*(that**2+uhat**2)/shat**2*((N**2-1)**2/2/N**3/uhat/that*(that+uhat)**2-(N**2-1)/N))

            # Total scattering amplitude
            MSquareTotal  = MSquareString + MSquareQCD

            #Convolution of the scattering amplitudes with parton distribution function
            return(qqbar2gg11*(DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale) + DistFunc.xfxQ(FstpartonID,xb,PDFScale)*DistFunc.xfxQ(ScndpartonID,xa,PDFScale))*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

        #gg->gGamma
        def MonteCarlogggGamma(parameters, *args):

            Y, y                           = parameters
            M, FstpartonID, ScndpartonID   = args
            
            tau = M**2/s
            
            xa  = math.sqrt(tau)*math.exp( Y)
            xb  = math.sqrt(tau)*math.exp(-Y)

            shat = M**2
            that = -0.5*M**2*math.exp(-y)/math.cosh(y)
            uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

            # Production of the first string resonance
            MSquareTotal  = StringCoeff11*5*g3**4*Q**2/3/Ms**4*(Ms**8/((shat-Ms**2)**2 + (GJ0gS*Ms)**2) + (that**4 + uhat**4)/((shat-Ms**2)**2 + (GJ2gS*Ms)**2))

            #Convolution of the scattering amplitudes with parton distribution function
            return(2*gg2gGamma11*DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale)*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

        #gq->qGamma
        def MonteCarlogqqGamma(parameters, *args):

            Y, y                           = parameters
            M, FstpartonID, ScndpartonID   = args
            
            tau = M**2/s
            
            xa  = math.sqrt(tau)*math.exp( Y)
            xb  = math.sqrt(tau)*math.exp(-Y)

            shat = M**2
            that = -0.5*M**2*math.exp(-y)/math.cosh(y)
            uhat = -0.5*M**2*math.exp( y)/math.cosh(y)

            # Production of the first string resonance
            MSquareTotal  = -StringCoeff11*g3**4*Q**2/3/Ms**2*(Ms**4*uhat/((shat-Ms**2)**2 + (GJhfqS*Ms)**2) + uhat**3/((shat-Ms**2)**2 + (GJ3hfqS*Ms)**2))

            #Convolution of the scattering amplitudes with parton distribution function
            return(gq2qGamma11*(DistFunc.xfxQ(FstpartonID,xa,PDFScale)*DistFunc.xfxQ(ScndpartonID,xb,PDFScale) + DistFunc.xfxQ(FstpartonID,xb,PDFScale)*DistFunc.xfxQ(ScndpartonID,xa,PDFScale))*MSquareTotal*M/(math.cosh(y)**2*16*math.pi*shat**2))

    # Function to numerically integrate
    def Integrate(Func, args, N, a, b, x):
        value = 0
        for i in range(1, N+1):
            Y = a+((i-(1/2))*((b-a)/N))
            c = -(x - math.fabs(Y))
            d = x - math.fabs(Y)
            for j in range(1, N+1):
                y = c+((j-(1/2))*((d-c)/N))
                params = [Y, y]
                value += ((b-a)/N)*((d-c)/N)*Func(params, *args)
        return(value)
    
    
    # Number of integration slices
    MC = 5
    
    # Mass weighting
    DeltaMass = (MaxMass - MinMass)/Number
    
    # Assigning the parameters that are used in the LHE file
    IDBMUP = 2212              # Proton beams
    EBMUP  = COME/2            # Energy of each proton beam
    PDFGUP = 0
    PDFSUP = DistFunc.lhapdfID # ID of the PDF set
    IDWTUP = 3                 # Unweighted events
    XWGTUP = 1                 # Unweighted events
    NUP = 4                    # Number of particles
    
    # Finding the maximum cross sections
    # Max and min parameters
    tauMin = MinMass**2/s
    YmaxMin = min([math.log(1/math.sqrt(tauMin)), yMax])
    
    if(Ms < MinMass or Ms > MaxMass):
        Msss = MinMass
    else:
        Msss = Ms
       
    tauMax = Msss**2/s
    YmaxMax = min([math.log(1/math.sqrt(tauMax)), yMax])
	

    
    
    tau = mass**2/s
    Ymax = min([math.log(1/math.sqrt(tau)), yMax])
    
    total_list = []
    params_list = []
   
    
    # yMax_range = [i*0.01 for i in range(0, 250, 1)]
    
    # Ymax_range = [min([math.log(1/math.sqrt(tau)), i]) for i in yMax_range]
    
    # params_list=list(zip(Ymax_range, yMax_range))
    
    
    N = 100
    a, b, x = -Ymax, Ymax, yMax
    for i in range(1, N+1):
        Y = a+((i-(1/2))*((b-a)/N))
        c = -(x - math.fabs(Y))
        d = x - math.fabs(Y)
        for j in range(1, N+1):
            y = c+((j-(1/2))*((d-c)/N))
            if math.fabs(Y) < Ymax and math.fabs(y) < yMax:
                params_list.append([Y, y])
            
    # print(params_list)
    
    valuesDict = {}
          
        
      
    for params_k in params_list:
    
        sumgggg      = MonteCarlogggg(params_k, mass, 21, 21) 

        sumggqqbar   = MonteCarloggqqbar(params_k, mass, 21, 21)

        sumgdgd = MonteCarlogqgq(params_k, mass, 21, 1) 
        sumgugu = MonteCarlogqgq(params_k, mass, 21, 2) 
        sumgsgs = MonteCarlogqgq(params_k, mass, 21, 3) 
        sumgcgc = MonteCarlogqgq(params_k, mass, 21, 4) 
        sumgbgb = MonteCarlogqgq(params_k, mass, 21, 5) 
        sumgtgt = MonteCarlogqgq(params_k, mass, 21, 6) 

        sumgdbargdbar = MonteCarlogqbargqbar(params_k, mass, 21, -1) 
        sumgubargubar = MonteCarlogqbargqbar(params_k, mass, 21, -2) 
        sumgsbargsbar = MonteCarlogqbargqbar(params_k, mass, 21, -3) 
        sumgcbargcbar = MonteCarlogqbargqbar(params_k, mass, 21, -4) 
        sumgbbargbbar = MonteCarlogqbargqbar(params_k, mass, 21, -5) 
        sumgtbargtbar = MonteCarlogqbargqbar(params_k, mass, 21, -6) 

        sumddbargg = MonteCarloqqbargg(params_k, mass, 1, -1)
        sumuubargg = MonteCarloqqbargg(params_k, mass, 2, -2)
        sumssbargg = MonteCarloqqbargg(params_k, mass, 3, -3)
        sumccbargg = MonteCarloqqbargg(params_k, mass, 4, -4)
        sumbbbargg = MonteCarloqqbargg(params_k, mass, 5, -5)
        sumttbargg = MonteCarloqqbargg(params_k, mass, 6, -6)
        
        value = sumgggg+sumggqqbar+sumgdgd+sumgugu+sumgsgs+sumgcgc+sumgbgb+sumgtgt+sumgdbargdbar+sumgubargubar+sumgsbargsbar+sumgcbargcbar+sumgbbargbbar+sumgtbargtbar+sumddbargg+sumuubargg+sumssbargg+sumccbargg+sumbbbargg+sumttbargg
        if params_k[0] not in list(valuesDict.keys()):
            valuesDict[params_k[0]] = {}
        valuesDict[params_k[0]][params_k[1]] = value
    return valuesDict
    
args = parsing()

lowerBound = args.lowerBoundvalue
upperBound = args.upperBoundvalue
stepSize = args.stepSizevalue

if checker() == True:

    mass = 7000
    valuesDict = main(args.RandGenSEEDvalue, args.Numbervalue, args.Msvalue, args.COMEvalue, args.MinMassvalue, args.MaxMassvalue, args.yMaxvalue, args.PDFSetvalue, args.PDFScalevalue, args.Couplingvalue, args.CouplingScalevalue, args.FirstStringCoeffvalue, args.SecondStringCoeffvalue, args.QCDCoeffvalue, args.dMassvalue, args.uMassvalue, args.sMassvalue, args.cMassvalue, args.bMassvalue, args.tMassvalue, args.gg2ggvalue, args.gg2qqbarvalue, args.gq2gqvalue, args.gqbar2gqbarvalue, args.qqbar2ggvalue, args.gg2gGammavalue, args.gq2qGammavalue, mass)
    
    print(valuesDict)
    
    # with open("./list.txt", "a") as file:
        # file.write(" ".join(str(x) for x in valuesDict))

else:
    print('Error')
    

