# Q&A for MSI-CODEC

Thanks @Luca Mariani from Copenhagen University Hospital for asking these questions: 

1. What are the output from the CODEC msi c++ code
   
   A: `.all` and `.msi`. `.all` contains all identified homopolymers and STR in the sample, where `.msi` only contains sites with deletion(s) (compared to reference) in the tumor sample.

2. What does each column represent in the output. For example for a row as
   
   ```chr1  1195500  T  18  15  18,18  0.3  0.5  0.000001```

   A: Column from left to right represents chromosome, position, reference allele/single repeat unit, length of the repeat unit (in this case a homopolymer of T with a length of 18),
   tumor sample fragment count (a single fragment/observation of T15), normal sample fragment count (2 observations of T18),
   frequency of T in the previous 10bp window of the reference genome and frequency of T in the next 10bp window of the reference genome.
   The last column is this population VAF, which was used to filter common alleles in the population. 
3. Where does the reference file `msi_profile.json` originate from and what does it represent?

   A: The msi_profile.json was generated from the CODEC NA12878 sample published in the Nat Genet paper. It represents the CODEC error profile, as shown in Figure 5 of that publication.  

4. How is the final MSI score being calculated

   A: This calculation of the MSI scores is in the `msi.R` file. If you look at the last chunk of the R code (after `###LOAD ALL DATA!!!!!!!!!!!!` line),
   there is an example of how to run it on a test data set. I've also added a test data set in this directory.
   You can download it and replace the input paths in the R script with your local paths. 
