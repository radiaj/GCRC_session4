# GCRC_session4

Make sure to use the most recent version of R (version )
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
