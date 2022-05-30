export PATH=$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH


lowerBound=1000 #Lower bound (GeV) for mass
upperBound=9000 #Upper bound (GeV) for mass 
stepSize=20 #step size

RandGenSEED=$(shuf -i 1-100000000000 -n 1 ) #### A Seed for the Random Number Generator
COME=13000            #### Centre of Mass Energy of the Proton-Proton Collision (GeV).
Ms=7000               #### String Scale (GeV) ==> This is the Same Scale for QCD.
MinMass=8050          #### Lower Bound for the Invariant Mass (GeV).
MaxMass=13000          #### Upper Bound for the Invariant Mass (GeV). The Maximum Value for this Variable is COME. The Default Value is the COME
yMax=2.5              #### Upper Bound for the Rapidity of the Outgoing Partons  |y1|,|y2| < ymax.
Number=1          #### Number of Events to be Generated.
PDFSet=cteq6l1        #### Parton Distribution Function (PDF) Set. Other PDF Sets that Can be Set Here: "CT14lo", "CT10"
Coupling=-1           #### The Running Coupling Constant, alpha_s (= g^2/4*pi). In Order to Use a Fixed Coupling Constant, Simply Use the Intended Value.
CouplingScale=$Ms     #### Scale at Which Coupling Constant is Calculated (Only if the Running Coupling is Being Used). The Default Value for this Variable is the String Scale Ms.
PDFScale=$Ms          #### Scale at which PDFs are Determined. The Default Value for this Variable is the QCD Scale "CouplingScale".

################### SubProcessess #################

gg2gg=true       #### gg -> gg       Subprocess (ID = 1)
gg2qqbar=True    #### gg -> qqbar    Subprocess (ID = 2)
gq2gq=True       #### gq -> gq       Subprocess (ID = 3)
gqbar2gqbar=True #### gqbar -> gqbar Subprocess (ID = 4)
qqbar2gg=True    #### qqbar -> gg    Subprocess (ID = 5)
gg2gGamma=False  #### gg -> gGamma   Subprocess (ID = 6)
gq2qGamma=False  #### gq -> qGamma   Subprocess (ID = 7)

###################################################

QCDCoeff=false           #### Production of QCD tree-level diparton
FirstStringCoeff=true    #### Production of First  String Resonance ( 2 --> 2 Partonic Scattering and 2-Parton --> Parton-Gamma Scattering )
SecondStringCoeff=False  #### Production of Second String Resonance ( 2 --> 2 Partonic Scattering )

################## Qurak Masses ###################

dMass=5e-3    #### GeV
uMass=2e-3    #### GeV
sMass=1e-3    #### GeV
cMass=1.25    #### GeV
bMass=4.2     #### GeV
tMass=172.5   #### GeV

####################################################

python ./differential.py -RandGenSEED $RandGenSEED \
			  -Number $Number \
			  -Ms $Ms \
			  -COME $COME \
			  -MinMass $MinMass \
			  -MaxMass $MaxMass \
			  -yMax $yMax \
			  -PDFSet $PDFSet \
			  -PDFScale $PDFScale \
			  -Coupling $Coupling \
			  -CouplingScale $CouplingScale \
			  -FirstStringCoeff $FirstStringCoeff  \
			  -SecondStringCoeff $SecondStringCoeff \
			  -QCDCoeff $QCDCoeff \
			  -dMass $dMass \
			  -uMass $uMass \
			  -sMass $sMass \
			  -cMass $cMass \
			  -bMass $bMass \
			  -tMass $tMass \
			  -gg2gg $gg2gg \
			  -gg2qqbar $gg2qqbar \
			  -gq2gq $gq2gq \
			  -gqbar2gqbar $gqbar2gqbar \
			  -qqbar2gg $qqbar2gg \
			  -gg2gGamma $gg2gGamma \
			  -gq2qGamma $gq2qGamma \
			  -lowerBound $lowerBound \
			  -upperBound $upperBound \
			  -stepSize $stepSize \
