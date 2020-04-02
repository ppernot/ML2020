# Install packages if necessary ####
libs = c(
  "xtable","mixtools","inlmisc","rlist","boot",
  "repmis","sessioninfo","moments"
)
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(
      lib,
      dependencies = TRUE,
      repos = "https://cran.univ-paris1.fr"
    )
  }
}
# Load packages and generate biblio
repmis::LoadandCite(libs,file='../article/packages.bib')

# Parallel options for bootstrap
options(boot.parallel = "multicore")
options(boot.ncpus = 4)

# Load misc. functions
source('../analysis/functions.R')

## Github package
lib = "ErrViewLib"
if(!require(lib,character.only = TRUE))
  devtools::install_github(paste0("ppernot/",lib))
library(lib,character.only = TRUE)


# Set graphical params ####
gPars = list(
  cols     = rev(inlmisc::GetColors(8))[1:7],
  cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.2))[1:7],
  cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.75))[1:7],
  pty      = 's',
  mar      = c(3,3,1.6,.2),
  mgp      = c(2,.75,0),
  tcl      = -0.5,
  lwd      = 4.0,
  cex      = 4.0,
  cex.leg  = 0.7,
  reso     = 1200  # (px) base resolution for png figs
)
# Expose gPars list
for (n in names(gPars))
  assign(n, rlist::list.extract(gPars, n))

# Define Data and Results repositories
dataRepo = '../data/'
figRepo  = '../results/figs/'
tabRepo  = '../results/tables/'

sink(file ='./sessionInfo.txt')
print(sessioninfo::session_info())
sink()

# Load data ####
composTab = read.csv(
  file = paste0(dataRepo, 'ComposTab.csv'),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
rownames(composTab) = composTab[,1]
composTab = composTab[,-1]
composTabLearn = composTab[1:1000,]
composTabTest = composTab[1001:nrow(composTab),]
data = read.csv(
  file = paste0(dataRepo, 'Data.csv'),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Calculate error sets ####
idTest = data[,1]
rownames(data) = idTest
Ref = data[,2]
Data = data[,-c(1,2)]
methList = colnames(Data)
errors = Ref - Data

stats = ErrViewLib::estBS1(
  errors,
  props = c('mue','mse','rmsd','q95hd','P1','skew','kurt','W'),
  eps = 1,
  do.sip = FALSE)
df1 = ErrViewLib::genTabStat(
  stats,
  comp = FALSE,
  units = 'kcal/mol'
)

# Linear correction ####
errors_lc = errors
for(m in methList) {
  reg = lm(Ref ~ Data[,m])
  errors_lc[,m] = residuals(reg)
}

stats_lc = ErrViewLib::estBS1(
  errors_lc,
  props = c('mue','mse','rmsd','q95hd','P1','skew','kurt','W'),
  eps = 1,
  do.sip = FALSE)
df2 = ErrViewLib::genTabStat(
  stats_lc,
  comp = FALSE,
  units = 'kcal/mol'
)

# Outliers by mixture analysis ####

## Fit SLATM-L2 distribution as bi-normal mixture
err = errors[,'SLATM-L2']
mixmdl = mixtools::normalmixEM(err, k=2)

sink(paste0(tabRepo,'biNormal.tex'))
df= data.frame(
  mu    = signif(mixmdl$mu,2),
  sigma = signif(mixmdl$sigma,2),
  weight= signif(mixmdl$lambda,2)
)
print(
  xtable::xtable(
    df,
    type = 'latex',
    caption = 'Parameters of the bi-normal fit of SALTM-L2 errors',
    label = "tab:MLBiNorm"
  ),
  comment = FALSE,
  include.rownames = FALSE,
  caption.placement ='bottom'
)
sink()

# Detect outliers ####
i1 = which.min(mixmdl$sigma) # id of most concentrated compnt
mu1    = mixmdl$mu[i1]
sigma1 = mixmdl$sigma[i1]

### 5*sigma bounds on larger dist
lims5s    = mu1 + 5*sigma1*c(-1,1)
sel5s     = (err-lims5s[1])*(err-lims5s[2]) > 0
outList5s = idTest[sel5s]
errors5s  = err[sel5s]
outCompos = composTabTest[sel5s,]

## Info table of largest outliers ####
uniqueComp = outCompos[!duplicated(outCompos),]
testIso = outIso = isoMass = learnIso = MAE = MaxE = c()
for (i in 1:nrow(uniqueComp)) {
  ref = uniqueComp[i,]
  testIso[i]  = sum(apply(composTabTest ,1,function(x) all(x == ref)))
  outIso[i]   = sum(apply(outCompos ,1,function(x) all(x == ref)))
  learnIso[i] = sum(apply(composTabLearn,1,function(x) all(x == ref)))
  selo    = which(apply(outCompos ,1,function(x) all(x == ref)))
  MAE[i]  = median(abs(errors5s[selo]))
  MaxE[i] = max(abs(errors5s[selo]))
}

dbe = as.integer(apply(uniqueComp,1,dbeCalc))
io = order(MaxE, decreasing = TRUE)
molF = c()
for (i in 1:nrow(uniqueComp)) {
  x = uniqueComp[i,]
  sel = which(x != 0)
  molF[i] = paste0('\\ce{',paste0(names(x)[sel],x[sel],collapse=''),'}')
}

# Save Table II ####
sink(paste0(tabRepo,'5sOutliersCompo.tex'))
df= data.frame(
  name = molF,
  DBE = dbe,
  medAE = MAE,
  maxAE = MaxE,
  'nOutl/nTest'   = paste0(outIso,'/',testIso),
  nLearn  = learnIso
)[io,]
print(
  xtable::xtable(
    df,
    type = 'latex',
    caption = 'Composition of the SLATM-L2 5*sigma outliers',
    label = "tab:ML5soutl",
    digits = c(0,0,0,2,2,0,0)
  ),
  comment = FALSE,
  include.rownames = FALSE,
  caption.placement ='bottom'
)
sink()

# Load data for augmented learning sets ####
## Outliers augmented LS
data_wo = read.csv(
  file = paste0(dataRepo, 'Data_wo.csv'),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
idTest_wo = data_wo[,1]
rownames(data_wo) = idTest_wo
Ref_wo = data_wo[,2]
Data_wo = data.frame(data_wo[,-c(1,2)])
colnames(Data_wo) = colnames(data_wo)[3]
rownames(Data_wo) = idTest_wo
errors_wo = Ref_wo - Data_wo

## Randomly augmented LS
data_wor = read.csv(
  file = paste0(dataRepo, 'Data_wor.csv'),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
idTest_wor = data_wor[,1]
rownames(data_wor) = idTest_wor
Ref_wor = data_wor[,2]
Data_wor = data.frame(data_wor[,-c(1,2)])
colnames(Data_wor) = colnames(data_wor)[3]
rownames(Data_wor) = idTest_wor
errors_wor = Ref_wor - Data_wor

# Original error set pruned from outliers
errors_no = matrix(err[!sel5s],ncol=1)
colnames(errors_no) = 'SLATM-L2(-o)'

stats_out = ErrViewLib::estBS1(
  cbind(errors_no,errors_wo,errors_wor),
  props = c('mue','mse','rmsd','q95hd','P1','skew','kurt','W'),
  eps = 1,
  do.sip = FALSE)
df_out = ErrViewLib::genTabStat(
  stats_out,
  comp = FALSE,
  units = 'kcal/mol'
)

df = rbind(df1,df2,df_out)
df = df[-c(6,11),] # Remove extra rows with units

# Save Table I ####
sink(paste0(tabRepo,'tabStatsML_BS.tex'))
print(
  xtable::xtable(
    df,
    type = 'latex',
    caption = 'Error summary statistics',
    label = "tab:MLStats"
  ),
  comment = FALSE,
  include.rownames = FALSE,
  caption.placement ='bottom'
)
sink()

# Fig. 1(a-f) ####
# Scatter plots and histograms
ifig = 0
for (key in names(errors)) {
  png(
    filename = paste0(figRepo, 'histDist_', key, '.png'),
    width = reso,
    height = reso
  )
  ErrViewLib::plotDistHist(
    Data[,key],
    errors[,key],
    xlab = expression(paste(italic(E)^{'*'},"(kcal/mol)")),
    ylab = 'Errors (kcal/mol)',
    main = key,
    xlim = range(Data),
    ylim = c(-30, 30),
    plotGauss = TRUE,
    plotConf = TRUE,
    scalePoints = 0.1,
    nclass = 55,
    gPars = gPars
  )
  ifig = ifig + 1
  mtext(
    text = paste0('(', letters[ifig], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.5
  )
  dev.off()
}

key = colnames(errors_wo)
png(
  filename = paste0(figRepo, 'histDist_', key, '.png'),
  width = reso,
  height = reso
)
ErrViewLib::plotDistHist(
  Data_wo[,key],
  errors_wo[,key],
  xlab = expression(paste(italic(E)^{'*'},"(kcal/mol)")),
  ylab = 'Errors (kcal/mol)',
  main = key,
  xlim = range(Data),
  ylim = c(-30, 30),
  plotGauss = TRUE,
  scalePoints = 0.1,
  nclass = 55,
  gPars = gPars
)
ifig = ifig + 1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

key = colnames(errors_wor)
png(
  filename = paste0(figRepo, 'histDist_', key, '.png'),
  width = reso,
  height = reso
)
ErrViewLib::plotDistHist(
  Data_wor[,key],
  errors_wor[,key],
  xlab = expression(paste(italic(E)^{'*'},"(kcal/mol)")),
  ylab = 'Errors (kcal/mol)',
  main = key,
  xlim = range(Data),
  ylim = c(-30, 30),
  plotGauss = TRUE,
  scalePoints = 0.1,
  nclass = 55,
  gPars = gPars
)
ifig = ifig + 1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

# Fig. 2(a-b)####
ifig = 1
png(
  filename = paste0(figRepo, 'compareECDF.png'),
  width = reso,
  height = reso
)
ErrViewLib::plotUncEcdf(
  abs(errors),
  xmax = 10.0,
  xlab = '|Errors| (kcal/mol)',
  title = '',
  show.leg = TRUE,
  show.MAE = TRUE,
  gPars = gPars
)
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

ifig = 2
png(
  filename = paste0(figRepo, 'compareECDF_wo.png'),
  width = reso,
  height = reso
)
lt0 = nrow(errors)
lt  = nrow(errors_wo)
fill = rep(NA, lt0 - lt)
X = data.frame(
  abs(errors[,c(2,4)]),
  c(abs(errors_wo[,1]),fill),
  c(abs(errors_wor[,1]),fill))
colnames(X) = c(
  methList[c(2,4)],
  colnames(errors_wo),
  colnames (errors_wor))
ErrViewLib::plotUncEcdf(
  X,
  xmax = 10,
  xlab = '|Errors| (kcal/mol)',
  title = '',
  show.leg = TRUE,
  show.MAE = TRUE,
  col.index = c(2, 4, 5, 7),
  gPars = gPars
)
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

# Fig. 4(a) ####
ifig = 1
meth = 'SLATM-L2'
png(
  file = paste0(figRepo, 'QQplot_', meth, '.png'),
  width = reso,
  height = reso
)
y = errors[,meth]
ErrViewLib::plotZscoreQqnorm(
  scale(y, center = TRUE, scale = TRUE),
  1,
  lim = 4,
  title = meth,
  gPars = gPars)
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

# Fig. 4(b) ####
## Bi-normal analysis
png(
  filename = paste0(figRepo, 'biNormalML_',meth, '.png'),
  width = reso,
  height = reso
)
par(
  mfrow = c(1, 1),
  mar = mar,
  mgp = mgp,
  pty = pty,
  tcl = tcl,
  cex = cex,
  lwd = lwd,
  lend = 2,
  yaxs = 'i'
)
hist(
  errors[,meth],
  freq = FALSE,
  col = cols_tr2[5],
  border = NA,
  nclass = 99,
  xlim = 12 * c(-1, 1),
  xlab = 'SLATM-L2 errors',
  ylim = c(0, 0.45),
  main = ''
)
grid()
abline(
  v = 5 * min(mixmdl$sigma) * c(-1, 1),
  lty = 2,
  col = cols[1]
)
x = seq(-15, 15, length.out = 1000)
for (ic in 1:length(mixmdl$lambda))
  lines(
    x,
    mixmdl$lambda[ic] * dnorm(x, mixmdl$mu[ic], mixmdl$sigma[ic]),
    col = cols[c(7, 6)][ic],
    lty = 2,
    lwd = 1.5 * lwd
  )
lines(
  x,
  mixmdl$lambda[1] * dnorm(x, mixmdl$mu[1], mixmdl$sigma[1]) +
    mixmdl$lambda[2] * dnorm(x, mixmdl$mu[2], mixmdl$sigma[2]),
  col = cols[2],
  lty = 1,
  lwd = 1.5 * lwd
)
box()
ifig = 2
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()


## nH/DBE analysis

scale = max(abs(errors[,meth]))

dbeAll   = as.integer(apply(composTab,      1, dbeCalc))
dbeLearn = as.integer(apply(composTabLearn, 1, dbeCalc))
dbeTest  = as.integer(apply(composTabTest,  1, dbeCalc))
dbeOut   = as.integer(apply(outCompos,      1, dbeCalc))

# Fig. 3(a-b)####
ifig = 0
for (meth in c('SLATM-L2','MP2')) {
  png(
    filename = paste0(figRepo, 'errStat_H_DBE_',meth, '.png'),
    width = reso,
    height = reso
  )
  plotnHDBE(
    x = composTabTest[, 'H'],
    y = dbeTest,
    z = abs(errors[,meth]),
    scale = scale,
    main = paste0(meth, ' abs. errors'),
    gPars = gPars
  )
  ifig = ifig + 1
  mtext(
    text = paste0('(', letters[ifig], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.5
  )
  dev.off()
}

# Fig. 3(c) ####
ifig = 3
png(
  filename = paste0(figRepo, 'errStat_H_DBE_SLATM-L2_Outliers.png'),
  width = reso,
  height = reso
)
plotnHDBE(
  x = composTabTest[, 'H'][sel5s],
  y = dbeTest[sel5s],
  z = abs(errors[,'SLATM-L2'])[sel5s],
  scale = scale,
  xlim = range(composTabTest[, 'H']),
  ylim = range(dbeTest),
  main = 'SLATM-L2 outliers',
  gPars = gPars
)
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

# Fig. 3(d) ####
## Errors per DBE class

ifig = 4
png(
  filename = paste0(figRepo, 'errDist_DBE.png'),
  width = reso,
  height = reso
)
par(
  mfrow = c(1, 1),
  mar = mar,
  mgp = mgp,
  pty = pty,
  tcl = tcl,
  cex = cex,
  lwd = lwd,
  xaxs = 'i',
  yaxs = 'i',
  lend = 2
)
x = dbeTest
y = abs(errors[,'SLATM-L2'])
plot(
  x,
  y,
  type = 'n',
  xlab = 'DBE',
  xlim = c(min(x - 0.5), max(x) + 0.5),
  ylab = '|Errors| (kcal/mol)',
  ylim = c(0, 1.1 * max(y)),
  main = 'SLATM-L2 abs. errors'
)
grid()
points(
  x + rnorm(length(x), 0, 0.1),
  y,
  pch = 16,
  cex = 0.5,
  col = cols_tr2[5]
)
abline(h = 5 * sigma1,
       col = cols[2],
       lty = 2)
box()
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()


# Fig. 5(a-d) ####

## DBE distributions
ifig = 0
png(
  filename = paste0(figRepo, 'DBELearnML.png'),
  width = reso,
  height = reso
)
plotBar(
  X = table(dbeLearn)[1:7],
  xlab = 'DBE',
  ylab = 'Distribution in Learning Set',
  mean = NULL,
  gPars
)
ifig = ifig + 1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

png(
  filename = paste0(figRepo, 'DBEValidML.png'),
  width = reso,
  height = reso
)
plotBar(
  X = table(dbeTest)[1:7],
  xlab = 'DBE',
  ylab = 'Distribution in Test Set',
  mean = NULL,
  gPars
)
ifig = ifig + 1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

png(
  filename = paste0(figRepo, 'ratioDBELearnML.png'),
  width = reso,
  height = reso
)
plotBar(
  X = table(dbeLearn)[1:7] / table(c(dbeTest, dbeLearn))[1:7],
  xlab = 'DBE',
  ylab = 'Proportion in Learning Set',
  mean = length(dbeLearn) / length(c(dbeTest, dbeLearn)),
  gPars
)
ifig = ifig + 1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

## Excess of outliers in DBE classes
png(
  filename = paste0(figRepo, 'excessDBEOutML.png'),
  width = reso,
  height = reso
)
X = 1 / table(dbeTest)[1:7] * table(dbeOut)[1:7]
X[!is.finite(X)] = 0
plotBar(
  X =   X,
  xlab = 'DBE',
  ylab = 'Proportion in Outliers',
  mean = length(dbeOut) / length(dbeTest),
  gPars
)
ifig = ifig + 1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

# Fig. 6(a,b) ####
ifig = 0
y = errors_wo[,1]
key = colnames(errors_wo)[1]
png(
  file = paste0(figRepo, 'QQplot_', key, '.png'),
  width = reso,
  height = reso
)
ErrViewLib::plotZscoreQqnorm(
  scale(y, center = TRUE, scale = TRUE),
  1,
  title = key,
  lim = 4,
  gPars = gPars)
ifig = ifig +1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

y = errors_wor[,1]
key = colnames(errors_wor)[1]
png(
  file = paste0(figRepo, 'QQplot_', key, '.png'),
  width = reso,
  height = reso
)
ErrViewLib::plotZscoreQqnorm(
  scale(y, center = TRUE, scale = TRUE),
  1,
  title = key,
  lim = 4,
  gPars = gPars)
ifig = ifig +1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

# Fig. 6(c-d) ####
ifig = 2
key = colnames(errors_wo)[1]
id = rownames(errors_wo)
X = composTab[paste0(id),]
x = X[, 'H']
y = apply(X, 1, dbeCalc)
png(
  filename = paste0(figRepo, 'errStat_H_DBE_', key, '.png'),
  width = reso,
  height = reso
)
plotnHDBE(
  x = x,
  y = y,
  z = abs(errors_wo[,1]),
  scale = scale,
  main = paste0(key, ' abs. errors'),
  gPars = gPars
)
ifig = ifig + 1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

key = colnames(errors_wor)[1]
id = rownames(errors_wor)
X = composTab[paste0(id),]
x = X[, 'H']
y = apply(X, 1, dbeCalc)
png(
  filename = paste0(figRepo, 'errStat_H_DBE_', key, '.png'),
  width = reso,
  height = reso
)
plotnHDBE(
  x = x,
  y = y,
  z = abs(errors_wor[,1]),
  scale = scale,
  main = paste0(key, ' abs. errors'),
  gPars = gPars
)
ifig = ifig + 1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.5
)
dev.off()

