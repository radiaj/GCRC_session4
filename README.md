# GCRC_session4

Make sure to use the most recent version of R (version >= 3.3.2) and the following Bioconductor packages are installed:
```R
source("https://bioconductor.org/biocLite.R")
biocLite("gage")
biocLite("gageData")
biocLite("pathview")
biocLite("clusterProfiler")
biocLite("ReactomePA")
biocLite("fgsea")
biocLite("DESeq2")
biocLite("biomaRt")

```
You will need to make sure you install clusterProfiler (version 3.2.11). You can install from clusterProfiler_3.2.11.tgz for MAC users or the clusterProfiler_3.2.11.zip for Windows users. 

```R
# To check the version installed load the library and run the sessionInfo() function
library("clusterProfiler")
sessionInfo()
```

# Windows Users
To update R version from within R terminal in Windows:
```R
# installing/loading the package:
if(!require(installr)) {
install.packages("installr"); require(installr)} #load / install+load installr
 
# using the package:
updateR() # this will start the updating process of your R installation.  It will check for newer versions, and if one is available, will guide you through the decisions you'd need to make.

```
see https://www.r-statistics.com/2013/03/updating-r-from-r-on-windows-using-the-installr-package/ for more information
