
<br>
<ol>
<li> **What is GR-shiny?** </li>
<p> GR-Shiny is a shiny app that allows users to benchmark various genomic-regions enrichment tools against their data. This comparison can be parametrized via metrics as sensitivity, specificity, precision, and prioritization. The following graphic represents the overall workflow and the functions in use. 
<br>
![](gr-shiny.png)
<br>
</p>
<br>
<br>

<li> **How is GR-shiny built?** </li>
<p> As a reference, kindly follow this [workflow](https://nbviewer.jupyter.org/github/shauryajauhari/GSABenchmarkTestAnalysis/blob/master/testProtocol.ipynb), or visit the [repository](https://github.com/shauryajauhari/GSABenchmarkTestAnalysis) for function definitions. P.S. the function definitions for the notebook have been tweaked a little to adapt to the Shiny App. 
</p>
<br>
<br>

<li> **What is the scope of the current set of methods available in this application?** </li>
<p> Due to their availability as functions in R, this application is limited to three tools currently for the end user, viz. [ChIPEnrich](http://chip-enrich.med.umich.edu/), [BroadEnrich](http://broad-enrich.med.umich.edu/), and [Seq2Pathway](https://www.bioconductor.org/packages/release/bioc/html/seq2pathway.html). However, in our analysis we have also compiled results explicitly from [GREAT](http://bejerano.stanford.edu/great/public/html/) and [Enrichr](http://amp.pharm.mssm.edu/Enrichr/) as well, as the two are prominent methods in use. These provide for contemporary options for testing enrichment in genomic regions. In the future, there is a distinct poosibility of including novel tools as part of an extended anslysis.
</p>
<br>
<br>

<li> **Are the gold standard datasets extendible?** </li>
<p> No. In it's current implementation, the gold standard datasets available for processing are *Colorectal Cancer*, *Gastric Cancer*, *Prostate Cancer*, and *Alzheimer's Disease*. 
</p>
<br>
<br>

<li> **What is the primary input for this pipeline?** </li>
<p> This application takes over from the path for the folder that holds the BED files for the genomic regions defining a sample. There could be multiple samples, although, the following conditions must hold. First, the BED files must be without a header. Second, the data must be in the basic tab-separated, BED format with *Chromosome*, *Start*, and *End* values; and third, the files must be saved with the **GEO sample name** as the primary name (preferably) and **.bed** extension. 
<br>
Initially when the app executes, the message box under 'Analyze Data' sidebar panel displays an error message attributing to no input from the user, as pictured below.
<br>
![](messageBoxError.png)
<br>
This is ephimeral. After the path entered by the user is recognized as valid, the samples are loaded and the error message washes off.
<br>
![](messageBoxErrorResolved.png)
</p>
<br>
<br>


<li> **How to track processing of the Shiny app?** </li>
<p> The featured *animated icons* represent the background working of the application. Exceptionally, they also signify that a wait is active for the user's input.
</p>
<br>
<br>


<li> **How can the results be downloaded?** </li>
<p> The results are saved as RDS objects inside the *./data/results/* folder of the working directory. The resultant plots (with apt nomenclature) are also saved as PNG graphics at the same location.
</p>
<br>
<br>


<li> **What if I have a method and I want to test it with the listed protocol?** </li>
<p> There could be two possible approaches to accomplish this.
  <ol>
    <li> **Add as a staple method (permanently)** </li>
    
    <p>
      If you want your method to be added, please send us the following information:
      * The list of packages that need to be loaded
      * The following variables: method-name, etc
      * The following functions: (a function that takes xxx input and produces yyy output)
      * etc"
      Then, we should just need to take this information and add it to the specific places in the code (such procedure should be docummented too).

    </p>
    
    <li> **Test a new method locally** </li>
    
    <p>
      If you want to test your method locally, without adding it permanently yet, you need the following procedure:
      * Download our jupyter/shiny app.
      * The list of packages that need to be loaded: YOU MUST SHOW THEM WHERE TO ADD THIS INFORMATION IN THE NOTEBOOK/APP
      * The following variables: method-name, etc: YOU MUST SHOW THEM WHERE TO ADD THIS INFORMATION IN THE NOTEBOOK/APP
      * The following functions: (a function that takes xxx input and produces yyy output): YOU MUST SHOW THEM WHERE TO ADD THIS INFORMATION
    
    </p>
  </ol>

</p>
<br>
</ol>