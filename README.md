# modSaRa2
Although it has been shown that the widely used change-point based methods can increase statistical power to identify variants, it remains challenging to effectively identify CNVs with weak signals due to the noisy nature of genotyping intensity data. modSaRa2 is a novel improvement of our previously developed method modified Screening and Ranking algorithm (modSaRa) by integrating the relative allelic intensity with prior information of statistics. modSaRa2 markedly improved both sensitivity and specificity over existing methods. The improvement for weak CNV signals is the most substantial, while simultaneously improving stability when CNV size varies. 
## Getting Started
The source package needs to be compiled first with the Rtools on Windows. 

For Mac, users need to have Xcode installed on the computer before installing modSaRa. To that end, the following steps are also recommended for installing the .tar.gz file.

1. Go to http://hpc.sourceforge.net

2. Download the newest version of gcc-x.x-bin.tar.gz

3. Run “gunzip gcc-x.x-bin.tar.gz” in the Mac terminal (change the directory to where this package was saved first)

4. Run “sudo tar -xvf gcc-x.x-bin.tar -C /” in the Mac terminal 
## Installing
```
install.packages("devtools")
library(devtools)
install_github("FeifeiXiaoUSC/modSaRa2",subdir="Package")
```
## Reference Manual
*[Reference Manual](https://github.com/FeifeiXiaoUSC/modSaRa2/blob/master/modSaRa2-manual.pdf)
