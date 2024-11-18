//Voriconazole PBPK model for a typical adult male
[CMT]

GUT 
ADIPOSE 
BRAIN 
HEART 
BONE 
KIDNEY 
LIVER 
LUNG 
MUSCLE 
SPLEEN 
REST 
ART 
VEN
URINE 
AUC // these are dummy variables for model validation. 


[PARAM] 
//Tissue volumes (L); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981
Vad      = 18.2   //adipose
Vbo      = 10.5   //bone
Vbr      = 1.45   //brain
VguWall  = 0.65   //gut wall
Vhe      = 0.33   //heart
Vki      = 0.31   //kidneys
Vli      = 1.8    //liver
Vlu      = 0.5    //lungs
Vmu      = 29     //muscle
Vsp      = 0.15   //spleen
Vbl      = 5.6    //blood
Vss 	   = 0.99 // (L/kg)

//Tissue blood flows (L/h); Cardiac output = 6.5 (L/min); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981
Qad = 0.05*6.5*60
Qbo = 0.05*6.5*60
Qbr = 0.12*6.5*60
Qgu = 0.15*6.5*60 
Qhe = 0.04*6.5*60
Qki = 0.19*6.5*60
Qmu = 0.17*6.5*60
Qsp = 0.03*6.5*60
Qha = 0.065*6.5*60  //hepatic artery
Qlu = 6.5*60        //same as cardiac output

//partition coefficients
Kpad = 0.193  //adipose:plasma
Kpbo = 0.467 // bone:plasma
Kpbr = 0.858  //brain:plasma
Kpgu = 0.773  //gut:plasma
Kphe = 0.814  //heart:plasma
Kpki = 0.84   //kidney:plasma
Kpli = 0.814  //liver:plasma
Kplu = 0.864  //lungs:plasma
Kpmu = 0.809  //muscle:plasma; optimized
Kpsp = 0.848  //spleen:plasma
BP = 1.15       //blood:plasma ratio; 

//other parameters
WEIGHT  = 70          //(kg)
fup     = 0.25         //fraction of unbound drug in plasma
GFR 	= 6 		// glomerular filtration rate (L/hr)

// dummy clearance parameter to scale up secretion
renal_secretion = 1.04  // https://pubmed.ncbi.nlm.nih.gov/35335919/

// hepatic clearance related rates; from https://pubmed.ncbi.nlm.nih.gov/35335919/ Table 2
biliary_excretion = 24.8 // [uL/min/1E6 cells]
// hepatic metabolic clearance 
HEP = 30.38  // intrinsic metabolic clearance; [uL/min/1E6 cells]
HLM = 86.29  // intrinsic metabolic clearance calculate per mg of miscrosomal protein; [uL/min/mg protein]
HLC = 38.74  // intrinsic metabolic clearance calculated per mg of cytosolic protein; [uL/min/mg protein]

// hepatic clearance scaling factor
hepatocyte_per_gram_liver = 139 // [1E6 cells per gram liver] https://pubmed.ncbi.nlm.nih.gov/16930941/
liver_weight = 1561 // [g]; https://pubmed.ncbi.nlm.nih.gov/22182984/
MPPGL = 30.3   // adult mg microsomal protein per g liver; [mg/g];  http://dmd.aspetjournals.org/content/38/1/25.long
CPPGL = 41     // cytosolic protein per gram of liver; [mg/g]; https://pubmed.ncbi.nlm.nih.gov/37735064/

CLintHep = 6  // [L/h]

$MAIN
// clearance 
double CLrenal = GFR + renal_secretion ; //(L/hr) 

// double CLintHep = ( ((biliary_excretion+HEP)*hepatocyte_per_gram_liver + HLM*MPPGL + HLC*CPPGL)*60/1E6*liver_weight )  ; 

// parition coefficient of the rest of the tissue; average of all non-adipose tissue Kps
double Kpre = ( Kpbo + Kpbr + Kpgu + Kphe + Kpki + Kpli + Kplu + Kpmu + Kpsp )/9; 

//additional volume derivations
double Vve = 0.705 * Vbl;         //venous blood
double Var = 0.295 * Vbl;         //arterial blood
double Vre = Vss * WEIGHT - (Vli + Vki + Vsp + Vhe + Vlu + Vbo + Vbr + Vmu + Vad + VguWall + Vbl); //volume of rest of the body compartment

//additional blood flow derivation
double Qli  = Qgu + Qsp + Qha;
double Qtot = Qli + Qki + Qbo + Qhe + Qmu + Qad + Qbr;
double Qre  = Qlu - Qtot;

$ODE
//Calculation of tissue drug concentrations (mg/L)
double Cadipose  = ADIPOSE / Vad;
double Cbone     = BONE / Vbo;
double Cbrain    = BRAIN / Vbr; 
double Cheart    = HEART / Vhe; 
double Ckidney   = KIDNEY / Vki;
double Cliver    = LIVER / Vli; 
double Clung     = LUNG / Vlu; 
double Cmuscle   = MUSCLE / Vmu;
double Cspleen   = SPLEEN / Vsp;
double Crest     = REST / Vre;
double Carterial = ART / Var;
double Cvenous   = VEN / Vve;
double Cgut      = GUT / VguWall;

//Free Concentration Calculations
double Cliverfree  = fup * Cliver; 
double Ckidneyfree = fup * Ckidney;

//ODEs

dxdt_GUT =  Qgu*(Carterial - Cgut/(Kpgu/BP)) ;

dxdt_ADIPOSE  = Qad * (Carterial - Cadipose / (Kpad / BP)); 

dxdt_BRAIN    = Qbr * (Carterial - Cbrain / (Kpbr / BP));

dxdt_HEART    = Qhe * (Carterial - Cheart / (Kphe / BP));

dxdt_KIDNEY   = Qki * (Carterial - Ckidney / (Kpki / BP)) - CLrenal * (Ckidneyfree / (Kpki / BP));

dxdt_LIVER    = Qgu * (Cgut / (Kpgu / BP)) 
				+ Qsp * (Cspleen / (Kpsp / BP)) 
				+ Qha * (Carterial)
				- Qli * (Cliver / (Kpli / BP))
				// - CLintHep * (Cliverfree / (Kpli / BP)
				- Qli * (Cliverfree * CLintHep)/( Qli + Cliverfree * CLintHep ) ; 

dxdt_LUNG     = Qlu * (Cvenous - Clung / (Kplu / BP));

dxdt_MUSCLE   = Qmu * (Carterial - Cmuscle / (Kpmu / BP));

dxdt_SPLEEN   = Qsp * (Carterial - Cspleen / (Kpsp / BP));

dxdt_BONE     = Qbo * (Carterial - Cbone / (Kpbo / BP));

dxdt_REST     = Qre * (Carterial - Crest / (Kpre / BP));

dxdt_VEN      = Qad * (Cadipose / (Kpad / BP)) + 
  Qbr * (Cbrain / (Kpbr / BP)) +
  Qhe * (Cheart / (Kphe / BP)) + 
  Qki * (Ckidney / (Kpki / BP)) + 
  Qli * (Cliver / (Kpli / BP)) + 
  Qmu * (Cmuscle / (Kpmu / BP)) + 
  Qbo * (Cbone / (Kpbo / BP)) + 
  Qre * (Crest / (Kpre / BP)) - 
  Qlu * Cvenous;

dxdt_ART      = Qlu * (Clung / (Kplu / BP) - Carterial);

// calculate for the dummy variables
dxdt_AUC = Cvenous/BP; 

dxdt_URINE = CLrenal * (Ckidneyfree / (Kpki / BP)); 


$CAPTURE
Cvenous 