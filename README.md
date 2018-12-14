## Summary
The purpose of this project is to locate mCHH islands within several grass genomes. 

## Background
  In plants, cytosines can be methylated in three sequence contexts: CG, CHG, and CHH where H = any base other than G. One main purpose of methylation is to silence transposable elements (TEs), which are mostly deleterious. Plants methylate cytosines within TEs in all three of the sequence contexts, and this appears to keep them in a heterochromatic state where they aren’t expressed and can’t proliferate. Strangely, genes are also often highly methylated with the opposite outcome. Highly methylated genes are generally moderately expressed across a broad range of tissues. However, while methylated TEs have cytosines methylated in all three contexts, genes are generally only highly methylated in the CG context. 
  Looking across the genome, one can imagine long sections of TEs tightly wrapped in heterochromatin with CG, CHG, and CHH cytosines methylated, and genic sections in euchromatin with only CG cytosines methylated. An additional epigenetic feature in this genomic landscape was recently discovered in maize, where genes usually have small regions (~100 bp) dense with methylated CHH cytosines before the transcription start site and after the terminator. Presence of these regions near genes, called mCHH islands, has been associated with more highly expressed genes in maize (Li et al. 2015). mCHH islands have also been found in many other species as well, but it was unclear whether the association between island presence and expression was universal (Niederhuth et al. 2016).
  A post doc in my lab, Dr. Danelle Seymour, is currently writing a paper examining gene body methylation across 8 grass species. So, I’m using her bisulfite sequencing and RNA-seq data to investigate mCHH islands and try to associate them with expression. At the moment, I'm still trying to figure out the best way to define these islands (both in terms of what proportion of CHH sites need to be methylated and how to associate these regions with genes) but for this project I have used the definition which has been used previously in mCHH island literature. In this rigid (and arbitrary) definition, mCHH islands are 100 bp regions with >25% mCHH.
  