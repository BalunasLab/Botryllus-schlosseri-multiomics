# Multi-omics Analysis Reveals Important Role for Microbial-derived Metabolites from Botryllus schlosseri in Metal Interactions

## Description
All the analyses done for this manuscript were done using open source R code.
The description and how to use of each package can be found in their respective Githubs.

##installation
library("qiime2R", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("vegan", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("dplyr")
library("scales", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
install.packages(c("UpSetR","mltools","data.table","VennDiagram","grid","futile.logger","vegan","geosphere","mixOmics", "IMIFA", "compositions","phyloseq","microbiome","igraph", "associationsubgraphs","tidyr","stringr))

## References:
Bolyen, E., Rideout, J. R., Dillon, M. R., Bokulich, N. A., Abnet, C., Al-Ghalith, G. A., Alexander, H., Alm, E. J., Arumugam, M., Asnicar, F., Bai, Y., Bisanz, J. E., Bittinger, K., Brejnrod, A., Brislawn, C. J., Brown, C. T., Callahan, B. J., Caraballo-Rodríguez, A. M., Chase, J., … Caporaso, J. G. (2018). QIIME 2: Reproducible, interactive, scalable, and extensible microbiome data science. PeerJ Preprints, 6, e27295v1. https://doi.org/10.7287/peerj.preprints.27295v1
Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R., & Pfister, H. (2014). UpSet: Visualization of Intersecting Sets. IEEE Transactions on Visualization and Computer Graphics, 20(12), 1983–1992. IEEE Transactions on Visualization and Computer Graphics. https://doi.org/10.1109/TVCG.2014.2346248
Rohart, F., Gautier, B., Singh, A., & Cao, K.-A. L. (2017). mixOmics: An R package for ‘omics feature selection and multiple data integration. PLOS Computational Biology, 13(11), e1005752. https://doi.org/10.1371/journal.pcbi.1005752
Strayer, N., Zhang, S., Yao, L., Vessels, T., Bejan, C. A., Hsi, R. S., Shirey-Rice, J. K., Balko, J. M., Johnson, D. B., Phillips, E. J., Bick, A., Edwards, T. L., Velez Edwards, D. R., Pulley, J. M., Wells, Q. S., Savona, M. R., Cox, N. J., Roden, D. M., Ruderfer, D. M., & Xu, Y. (2023). Interactive network-based clustering and investigation of multimorbidity association matrices with associationSubgraphs. Bioinformatics, 39(1), btac780. https://doi.org/10.1093/bioinformatics/btac780

