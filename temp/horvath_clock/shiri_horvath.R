# using
# https://horvath.genetics.ucla.edu/html/dnamage/

install.packages("WGCNA")
install.packages("sqldf")
#install the Bioconductor installer
install.packages("BiocInstaller",repos="http://www.bioconductor.org/packages/2.13/bioc")

#source("http://bioconductor.org/biocLite.R") # doesn't work R 3.5 need to do with Manager
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("biocLite")

#biocLite("impute") # doesn't work
#install “impute” from Bioconductor
BiocManager::install("impute")

#############

# LOAD FILE
# Copy and pase the following R software code
# Use forward slashes /”, as R will misread filepaths with backslashes
#setwd("C:/Users/SHorvath/Documents/DNAmAge/Example55")
setwd("./horvath_clock")

library(WGCNA)
library(sqldf)

source("NORMALIZATION.R")

#Age transformation and probe annotation functions
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("datMiniAnnotation27k.csv")
datClock=read.csv("AdditionalFile3.csv")

#Read in the DNA methylation data (beta values)
# For a small file, e.g. measured on the 27k platform you could just use read.csv.
# But for large files, e.g. those measured on the 450K platform, I recommend you use read.csv.sql.
dat0=read.csv.sql("beta1_padded.csv") ; # misha here we changed
nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]

# the following command may not be needed. But it is sometimes useful when you use read.csv.sql
dat0[,1]= gsub(x=dat0 [,1],pattern="\"",replacement="")

#Create a log file which will be output into your directory
# The code looks a bit complicated because it serves to create a log file (for error checks etc).
# It will automatically create a log file.
file.remove("LogFile.txt")
file.create("LogFile.txt")
DoNotProceed=FALSE
cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ", nProbes, " probes."),file="LogFile.txt")
if (nSamples==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be no
samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql .
Samples correspond to columns in that file ."), file="LogFile.txt",append=TRUE) }
if (nProbes==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be zero
probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql
CpGs correspond to rows.") , file="LogFile.txt",append=TRUE) }
if ( nSamples > nProbes ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG
probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose
the data and then resubmit them? In any event, I will proceed with the analysis."),file="LogFile.txt",append=TRUE) }
if ( is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe
identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file
contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1] ),file="LogFile.txt",append=TRUE) }
if ( !is.character(dat0[,1]) ) { cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg
numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe
identifiers such as cg00000292. Instead it contains ", dat0[1:3,1] ),file="LogFile.txt",append=TRUE) }
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."),
                  Comment=c("", "email Steve Horvath."))
if ( ! DoNotProceed ) {
  nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
  for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
  if ( sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain nonnumeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. 
4
Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure
this makes sense.\n" ),file="LogFile.txt",append=TRUE) }
  
  XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
  selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
  selectXchromosome[is.na(selectXchromosome)]=FALSE
  meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
  if ( sum(selectXchromosome) >=500 ) {
    meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
  if ( sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal
probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these
samples.\n " ),file="LogFile.txt",append=TRUE) }
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  if ( sum( is.na(match1))>0 ) {
    missingProbes= probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]
    DoNotProceed=TRUE; cat(paste( "\n \n Input error: You forgot to include the following ", length(missingProbes), " CpG probes
(or probe names):\n ", paste( missingProbes, sep="",collapse=", ")),file="LogFile.txt",append=TRUE) 
  }
  
  #STEP 2: Restrict the data to 21k probes and ensure they are numeric
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  if ( sum( is.na(match1))>0 ) stop(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
  dat1= dat0[match1,]
  asnumeric1=function(x) {as.numeric(as.character(x))}
  dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)
  
  #STEP 3: Create the output file called datout
  set.seed(1)
  # Do you want to normalize the data (recommended)?
  normalizeData=TRUE
  source("StepwiseAnalysis.txt")
  
  if ( sum( datout$Comment != "" ) ==0 ) { 
    cat(paste( "\n The individual samples appear to be fine."),file="LogFile.txt",append=TRUE) 
  }
  if ( sum( datout$Comment != "" ) >0 ) { 
    cat(paste( "\n Warnings were generated for the following samples.\n",
               datout[,1][datout$Comment != ""], "\n Hint: Check the output file for more details."),file="LogFile.txt",append=TRUE) 
  }
}
# output the results into the directory
write.table(datout,"Output.csv", row.names=F, sep="," )

#The methylation data set contains 8 samples (e.g. arrays) and  866091  probes.

#Input error: You forgot to include the following  1059  CpG probes
#(or probe names):
#  cg00005847, cg00012792, cg00019495, cg00020533, cg00033516, cg00061059, cg00071250, cg00089071, cg00090147, cg00101227, cg00109274, cg00113020, cg00145118, cg00209951, cg00226923, cg00230502, cg00234961, cg00253200, cg00263760, cg00298357, cg00308665, cg00328227, cg00338893, cg00396163, cg00402366, cg00433406, cg00445824, cg00463577, cg00485380, cg00514895, cg00556408, cg00559473, cg00567749, cg00613255, cg00613344, cg00638514, cg00702729, cg00731459, cg00736326, cg00762512, cg00777121, cg00818693, cg00821764, cg00826384, cg00840516, cg00849368, cg00930194, cg00940278, cg00948524, cg00953277, cg01062029, cg01069256, cg01072821, cg01130192, cg01139059, cg01148741, cg01161597, cg01172972, cg01197831, cg01204271, cg01240931, cg01265637, cg01278696, cg01297972, cg01299496, cg01301664, cg01320040, cg01354473, cg01372689, cg01373706, cg01413314, cg01413516, cg01422136, cg01429391, cg01431114, cg01447498, cg01476044, cg01484156, cg01485797, cg01485998, cg01522975, cg01585703, cg01594214, cg01605783, cg01638735, cg01684579, cg01703884, cg01708236, cg01722994, cg01743370, cg01785339, cg01802635, cg01805540, cg01819572, cg01856970, cg01868128, cg01888566, cg01888601, cg01900364, cg01971122, cg01988129, cg01991150, cg01995513, cg02008544, cg02066681, cg02077702, cg02092466, cg02110963, cg02196655, cg02200132, cg02227605, cg02235880, cg02248486, cg02254649, cg02266731, cg02301754, cg02314308, cg02320454, cg02377083, cg02441831, cg02564523, cg02582754, cg02611282, cg02643667, cg02647265, cg02654291, cg02706881, cg02776251, cg02780849, cg02841886, cg02882813, cg02885771, cg02905245, cg02910574, cg02915261, cg02917064, cg02932689, cg02936872, cg02941816, cg02949544, cg02972551, cg02975142, cg03016571, cg03029145, cg03032025, cg03070588, cg03082060, cg03098643, cg03109066, cg03112433, cg03147542, cg03149130, cg03159329, cg03166779, cg03190825, cg03201604, cg03214212, cg03266904, cg03271651, cg03332970, cg03343942, cg03364781, cg03366896, cg03373080, cg03410718, cg03415518, cg03494005, cg03573747, cg03582451, cg03591238, cg03611555, cg03642518, cg03648823, cg03680758, cg03724463, cg03835296, cg03841029, cg03849074, cg03856786, cg03875195, cg03890877, cg03940103, cg03954150, cg03958344, cg03986640, cg04009666, cg04020816, cg04101379, cg04136480, cg04137128, cg04180460, cg04293726, cg04300115, cg04304130, cg04394967, cg04405541, cg04422896, cg04431054, cg04435420, cg04437590, cg04451770, cg04466870, cg04473302, cg04520084, cg04531710, cg04542938, cg04551925, cg04583874, cg04600618, cg04612566, cg04668164, cg04678793, cg04700925, cg04711050, cg04720330, cg04726200, cg04754011, cg04785461, cg04790874, cg04826422, cg04881903, cg04894993, cg04951204, cg04970352, cg05028467, cg05055720, cg05103623, cg05157725, cg05284582, cg05352668, cg05461841, cg05462826, cg05538387, cg05554936, cg05588972, cg05590257, cg05636175, cg05764607, cg05768141, cg05788526, cg05793077, cg05801648, cg05824484, cg05882691, cg05896682, cg05927159, cg05940463, cg05942970, cg05949173, cg05950276, cg05957518, cg05958582, cg05959508, cg05976074, cg05990214, cg06037693, cg06055013, cg06088032, cg06117855, cg06131936, cg06141025, cg06161738, cg06162003, cg06194186, cg06200697, cg06211616, cg06216677, cg06236061, cg06236748, cg06263943, cg06394897, cg06434451, cg06438404, cg06462703, cg06494782, cg06540941, cg06564900, cg06614002, cg06625811, cg06655089, cg06661765, cg06707886, cg06750167, cg06754047, cg06763078, cg06768707, cg06824727, cg06834261, cg06851207, cg06862644, cg06866657, cg06885524, cg06911084, cg06942110, cg06948294, cg06963844, cg06981910, cg06994747, cg07007400, cg07016325, cg07018708, cg07042144, cg07050831, cg07054095, cg07122126, cg07154408, cg07175007, cg07181881, cg07197823, cg07233761, cg07242565, cg07380959, cg07383972, cg07408740, cg07447922, cg07457944, cg07465609, cg07473175, cg07476030, cg07478100, cg07490776, cg07505695, cg07512517, cg07536847, cg07550362, cg07611177, cg07671949, cg07675169, cg07685221, cg07749808, cg07750111, cg07845390, cg07895149, cg07896225, cg07899425, cg07912181, cg07947930, cg07959896, cg08007665, cg08015496, cg08023751, cg08082692, cg08096946, cg08097882, cg08097973, cg08146221, cg08151470, cg08156349, cg08164315, cg08222514, cg08244028, cg08341978, cg08390021, cg08393541, cg08395899, cg08397758, cg08398233, cg08412305, cg08422599, cg08492707, cg08551088, cg08568512, cg08574156, cg08604885, cg08651674, cg08660285, cg08661290, cg08687825, cg08708323, cg08724636, cg08728865, cg08738040, cg08748415, cg08797471, cg08811349, cg08860143, cg08907850, cg08914623, cg08918749, cg08986416, cg09014354, cg09018862, cg09059945, cg09072120, cg09144196, cg09164559, cg09173897, cg09199317, cg09209002, cg09234859, cg09242541, cg09251995, cg09272869, cg09305224, cg09325711, cg09350141, cg09384159, cg09473585, cg09499698, cg09504196, cg09508556, cg09518735, cg09558850, cg09572106, cg09588653, cg09660171, cg09671005, cg09674215, cg09715672, cg09729918, cg09757277, cg09785172, cg09787254, cg09788866, cg09793866, cg09863772, cg09869858, cg09879797, cg09906488, cg09978128, cg09991975, cg10065130, cg10082991, cg10095719, cg10108208, cg10134939, cg10174864, cg10180697, cg10216105, cg10368981, cg10415235, cg10467098, cg10493739, cg10503138, cg10563284, cg10569616, cg10575841, cg10581256, cg10588377, cg10617171, cg10624583, cg10643489, cg10645113, cg10660136, cg10730174, cg10784341, cg10792326, cg10829134, cg10844844, cg10845200, cg10866709, cg10894512, cg10910525, cg10912077, cg10936763, cg10957680, cg10962407, cg10982775, cg10984852, cg11034784, cg11053574, cg11054882, cg11054936, cg11062095, cg11075556, cg11078738, cg11104347, cg11149033, cg11157872, cg11196870, cg11219178, cg11375622, cg11418559, cg11419863, cg11422541, cg11460364, cg11532513, cg11547724, cg11592951, cg11601662, cg11653500, cg11654620, cg11822932, cg11825652, cg11879480, cg11899895, cg11973688, cg12024292, cg12038710, cg12040555, cg12061127, cg12076344, cg12080675, cg12085225, cg12095491, cg12124516, cg12128839, cg12149391, cg12287813, cg12291552, cg12292686, cg12335708, cg12346881, cg12374431, cg12376277, cg12393697, cg12422450, cg12453778, cg12460105, cg12504957, cg12515638, cg12530270, cg12569516, cg12585282, cg12598198, cg12601757, cg12648583, cg12706983, cg12717668, cg12718440, cg12792367, cg12818699, cg12867448, cg12879297, cg12936220, cg12940430, cg12949975, cg12971694, cg12977318, cg13051968, cg13098382, cg13132204, cg13191049, cg13283942, cg13294439, cg13326338, cg13351583, cg13405161, cg13439299, cg13488557, cg13493001, cg13500388, cg13500819, cg13559495, cg13560436, cg13581941, cg13585240, cg13599025, cg13678049, cg13679077, cg13682722, cg13724813, cg13745870, cg13754949, cg13760546, cg13945069, cg13979182, cg13998293, cg14023558, cg14031452, cg14035045, cg14036856, cg14059475, cg14060836, cg14063008, cg14155416, cg14155482, cg14265075, cg14265670, cg14279899, cg14288464, cg14297876, cg14303330, cg14329157, cg14334099, cg14391855, cg14405589, cg14414534, cg14419187, cg14427173, cg14541950, cg14542238, cg14605021, cg14614901, cg14623805, cg14643520, cg14646244, cg14651992, cg14700577, cg14718680, cg14776416, cg14809191, cg14832904, cg14869028, cg14893163, cg14961550, cg14982472, cg15070798, cg15128898, cg15146752, cg15158783, cg15188491, cg15201291, cg15214137, cg15280964, cg15301694, cg15334028, cg15341760, cg15361590, cg15372098, cg15377649, cg15379633, cg15383087, cg15383120, cg15384717, cg15394409, cg15410903, cg15413566, cg15427774, cg15439196, cg15489294, cg15492003, cg15538820, cg15562197, cg15563952, cg15654964, cg15692239, cg15699524, cg15733507, cg15741583, cg15767406, cg15777781, cg15796819, cg15806518, cg15840658, cg15865742, cg15869022, cg15898840, cg15973234, cg15983538, cg15984612, cg16008138, cg16041660, cg16046376, cg16097079, cg16101800, cg16176301, cg16176600, cg16196812, cg16199381, cg16267266, cg16267491, cg16282101, cg16284168, cg16293656, cg16321029, cg16326979, cg16330965, cg16348491, cg16381688, cg16382256, cg16408565, cg16417937, cg16425577, cg16443152, cg16444968, cg16449972, cg16467694, cg16483466, cg16490390, cg16494477, cg16495265, cg16515820, cg16555388, cg16558348, cg16579354, cg16580737, cg16597172, cg16607065, cg16636571, cg16648297, cg16715722, cg16721845, cg16743781, cg16752778, cg16762386, cg16763895, cg16833551, cg16908782, cg16927136, cg16954341, cg16957569, cg16979445, cg17036441, cg17095731, cg17119688, cg17155167, cg17171916, cg17185543, cg17186163, cg17189058, cg17240836, cg17264618, cg17272094, cg17296078, cg17306637, cg17391474, cg17391877, cg17408647, cg17465827, cg17470697, cg17498321, cg17505428, cg17524886, cg17525003, cg17530977, cg17575811, cg17607973, cg17648620, cg17710288, cg17718210, cg17793621, cg17797815, cg17814481, cg17824393, cg17838765, cg17854440, cg17864730, cg17869167, cg17897879, cg17918700, cg17947053, cg18023080, cg18053505, cg18089852, cg18159180, cg18189601, cg18250832, cg18312807, cg18313702, cg18329036, cg18354594, cg18374517, cg18396533, cg18433086, cg18436172, cg18463686, cg18474934, cg18502522, cg18503260, cg18519564, cg18638496, cg18652852, cg18705776, cg18718102, cg18728606, cg18773223, cg18783796, cg18805835, cg18824775, cg18855178, cg18870231, cg18878891, cg18896687, cg18906795, cg18936744, cg18959422, cg18987220, cg18995088, cg19046959, cg19096540, cg19109431, cg19118533, cg19118812, cg19167673, cg19216044, cg19220586, cg19221369, cg19250891, cg19273182, cg19278165, cg19280776, cg19286986, cg19292008, cg19370284, cg19434936, cg19524009, cg19564877, cg19569684, cg19591881, cg19649173, cg19721889, cg19728223, cg19744122, cg19746675, cg19793718, cg19795898, cg19855618, cg19945840, cg19947621, cg19948397, cg19951443, cg19955521, cg20050826, cg20070077, cg20083871, cg20092728, cg20118424, cg20158826, cg20326853, cg20339230, cg20358834, cg20366831, cg20394284, cg20427879, cg20443314, cg20452583, cg20461541, cg20483763, cg20509080, cg20537992, cg20566118, cg20573420, cg20649991, cg20655558, cg20695562, cg20715219, cg20718816, cg20723355, cg20744464, cg20752438, cg20757758, cg20790056, cg20811607, cg20854748, cg20903926, cg20904010, cg20909686, cg20918903, cg21017752, cg21032203, cg21032583, cg21033494, cg21068152, cg21094669, cg21179457, cg21222430, cg21229859, cg21240812, cg21255605, cg21272774, cg21304187, cg21365602, cg21546057, cg21549904, cg21610516, cg21613754, cg21626086, cg21639968, cg21665774, cg21684385, cg21688707, cg21692936, cg21716693, cg21794225, cg21896452, cg21912241, cg21926138, cg21943652, cg21968169, cg21988465, cg22063056, cg22085335, cg22131691, cg22156632, cg22188495, cg22201387, cg22209624, cg22213042, cg22220722, cg22318304, cg22333868, cg22380033, cg22409383, cg22467567, cg22546318, cg22552168, cg22560190, cg22572159, cg22580905, cg22601917, cg22606869, cg22610305, cg22680204, cg22778981, cg22805308, cg22807700, cg22874560, cg22874695, cg22880820, cg22902083, cg22909609, cg22925639, cg22933847, cg22940789, cg22994720, cg23001650, cg23032316, cg23033845, cg23036555, cg23043119, cg23074453, cg23088430, cg23128146, cg23141333, cg23211240, cg23234154, cg23297477, cg23301687, cg23378495, cg23423382, cg23430664, cg23442323, cg23472215, cg23495733, cg23527067, cg23580945, cg23583739, cg23640701, cg23683201, cg23727583, cg23807009, cg23851011, cg23867494, cg23889093, cg23897067, cg23916751, cg24022301, cg24076830, cg24107665, cg24127989, cg24164563, cg24292612, cg24346429, cg24362726, cg24399106, cg24422489, cg24441911, cg24471894, cg24473385, cg24497877, cg24576425, cg24579667, cg24641737, cg24659201, cg24685926, cg24695828, cg24777710, cg24793265, cg24801210, cg24816455, cg24826867, cg24868525, cg24884084, cg24901474, cg24903794, cg24914244, cg24970539, cg24974599, cg25017304, cg25025866, cg25060020, cg25097676, cg25104030, cg25115460, cg25140571, cg25186143, cg25201980, cg25221625, cg25300941, cg25306927, cg25336579, cg25374854, cg25391023, cg25411534, cg25421002, cg25438054, cg25443975, cg25476145, cg25477181, cg25482967, cg25522312, cg25589890, cg25596297, cg25764464, cg25778479, cg25820693, cg25842356, cg25885771, cg25902460, cg25915982, cg25937832, cg25946952, cg25957124, cg26057752, cg26079992, cg26093687, cg26124016, cg26174387, cg26267341, cg26279025, cg26286036, cg26289824, cg26304237, cg26331192, cg26353877, cg26390526, cg26391080, cg26411702, cg26413355, cg26422832, cg26471445, cg26519339, cg26523389, cg26590537, cg26605086, cg26626042, cg26643856, cg26656135, cg26677448, cg26682717, cg26783856, cg26789453, cg26817382, cg26840318, cg26869604, cg26887625, cg26974738, cg26990660, cg27016307, cg27016863, cg27038197, cg27043141, cg27138018, cg27147306, cg27164762, cg27180443, cg27190239, cg27190537, cg27196745, cg27210447, cg27265637, cg27287498, cg27319898, cg27329371, cg27486427, cg27554769, cg27557143, cg27588902






































#-OLD OLD OLD OLD OLD OLD ------------------------ now on our sample

install.packages("WGCNA")
install.packages("sqldf")
setwd("./horvath_clock")
library(WGCNA)
library(sqldf)
install.packages("BiocInstaller",repos="http://www.bioconductor.org/packages/2.13/bioc")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("biocLite")
BiocManager::install("impute")
source("NORMALIZATION.R")

#Age transformation and probe annotation functions
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("datMiniAnnotation27k.csv")
datClock=read.csv("AdditionalFile3.csv")

#Read in the DNA methylation data (beta values)
# For a small file, e.g. measured on the 27k platform you could just use read.csv.
# But for large files, e.g. those measured on the 450K platform, I recommend you use read.csv.sql.
#dat0=read.csv.sql("MethylationDataExample55.csv") ; # here we use the beta calculated in minif.R
dat0=beta
# original: nSamples=dim(dat0)[[2]]-1
nSamples=dim(dat0)[[2]]
nProbes= dim(dat0)[[1]]

# the following command may not be needed. But it is sometimes useful when you use read.csv.sql
dat0[,1]= gsub(x=dat0 [,1],pattern="\"",replacement="")

#Create a log file which will be output into your directory
# The code looks a bit complicated because it serves to create a log file (for error checks etc).
# It will automatically create a log file.
file.remove("LogFile.txt")
file.create("LogFile.txt")
DoNotProceed=FALSE
cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ", nProbes, " probes."),file="LogFile.txt")
if (nSamples==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be no
samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql .
Samples correspond to columns in that file ."), file="LogFile.txt",append=TRUE) }
if (nProbes==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be zero
probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql
CpGs correspond to rows.") , file="LogFile.txt",append=TRUE) }
if ( nSamples > nProbes ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG
probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose
the data and then resubmit them? In any event, I will proceed with the analysis."),file="LogFile.txt",append=TRUE) }
if ( is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe
identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file
contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1] ),file="LogFile.txt",append=TRUE) }
if ( !is.character(dat0[,1]) ) { cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg
numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe
identifiers such as cg00000292. Instead it contains ", dat0[1:3,1] ),file="LogFile.txt",append=TRUE) }
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."),
                  Comment=c("", "email Steve Horvath."))
if ( ! DoNotProceed ) {
  nonNumericColumn=rep(FALSE, nSamples)
  # original   for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
  for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(as.numeric(dat0[,i])) }
  if ( sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain nonnumeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. 
4
Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure
this makes sense.\n" ),file="LogFile.txt",append=TRUE) }
  XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
  #ORIGINAL  selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
  misha_names = dat0[,0]
  selectXchromosome=is.element(dat0[,0], XchromosomalCpGs )
  selectXchromosome[is.na(selectXchromosome)]=FALSE
  meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
  if ( sum(selectXchromosome) >=500 ) {
    meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
  if ( sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal
probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these
samples.\n " ),file="LogFile.txt",append=TRUE) } # it seems we don't have that so we can't predict age
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  
  
  #####  #####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####  here we have an error match1 is NA
  if ( sum( is.na(match1))>0 ) {
    missingProbes= probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]
    DoNotProceed=TRUE; cat(paste( "\n \n Input error: You forgot to include the following ", length(missingProbes), " CpG probes
(or probe names):\n ", paste( missingProbes, sep="",collapse=", ")),file="LogFile.txt",append=TRUE) 
  }
  
  #STEP 2: Restrict the data to 21k probes and ensure they are numeric
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  if ( sum( is.na(match1))>0 ) stop(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
  dat1= dat0[match1,]
  asnumeric1=function(x) {as.numeric(as.character(x))}
  dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)
  
  #STEP 3: Create the output file called datout
  set.seed(1)
  # Do you want to normalize the data (recommended)?
  normalizeData=TRUE
  source("StepwiseAnalysis.txt")
  
  if ( sum( datout$Comment != "" ) ==0 ) { 
    cat(paste( "\n The individual samples appear to be fine."),file="LogFile.txt",append=TRUE) 
  }
  if ( sum( datout$Comment != "" ) >0 ) { 
    cat(paste( "\n Warnings were generated for the following samples.\n",
               datout[,1][datout$Comment != ""], "\n Hint: Check the output file for more details."),file="LogFile.txt",append=TRUE) 
  }
}
# output the results into the directory
write.table(datout,"Output_Tirosh.csv", row.names=F, sep="," )


write.csv(beta, "beta.csv", quote=FALSE, row.names=TRUE)

