##### **What is GR-shiny?**
##### GR-Shiny is a shiny app that allows users to benchmark various genomic-regions enrichment tools against their data. This comparison can be parametrized via metrics as sensitivity, specificity, precision, and prioritization.

##### **How is GR-shiny built?**
##### For a reference workflow, kindly follow this [workflow](https://nbviewer.jupyter.org/github/shauryajauhari/GSABenchmarkTestAnalysis/blob/master/testProtocol.ipynb), or visit the [repository](https://github.com/shauryajauhari/GSABenchmarkTestAnalysis) for function definitions.

##### **What is the scope of the current set of methods available in this application?**
##### Due to their availability as functions in R, this application is limited to three tools currently for the end user, viz. [ChIPEnrich](http://chip-enrich.med.umich.edu/), [BroadEnrich](http://broad-enrich.med.umich.edu/), and [Seq2Pathway](https://www.bioconductor.org/packages/release/bioc/html/seq2pathway.html). However, in our analysis we have also compiled results explicitly from [GREAT](http://bejerano.stanford.edu/great/public/html/) and [Enrichr](http://amp.pharm.mssm.edu/Enrichr/) as well, as the two are prominent methods in use. These provide for contemporary options for testing enrichment in genomic regions. In the future, there is a distinct poosibility of including novel tools as part of an extended anslysis.

##### **Are the gold standard datasets extendible?**
##### No. In it's current implementation, the gold standard datasets available for processing are *Colorectal Cancer*, *Gastric Cancer*, *Prostate Cancer*, and *Alzheimer's Disease*. 

##### **What is the primary input for this pipeline?**
##### This application takes over from the path for the folder that holds the BED files for the genomic regions defining a sample. There could be multiple samples, although, the following conditions must hold. First, the BED files must be without a header. Second, the data must be in the basic tab-separated, BED format with *Chromosome*, *Start*, and *End* values; and third, the files must be saved with the **GEO sample name** as the primary name (preferably) and **.bed** extension. 