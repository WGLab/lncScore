# Introduction

lncScore is a python package for the identification of lncRNA from the assembled novel transcripts, and it also can be used to calculate the coding potential.

# Abstract

RNA-Seq based transcriptome assembly has been widely used in the identification of novel lncRNAs. However, it remains a challenge to identify lncRNAs from the assembled transcripts, particularly the partial-length ones, since partial-length protein-coding transcripts are more likely to be classified as lncRNAs due to their incomplete CDS. Furthermore, potential sequencing or assembly error that gain or abolish stop codons also complicates ORF-based prediction of lncRNAs. Here, we present a novel alignment-free tool, lncScore, which uses a logistic regression model with 11 carefully selected features derived from the open reading frame, exon, and the maximum coding subsequences. Compared to other state-of-the-art alignment-free tools (e.g. CPAT, CNCI, and PLEK), lncScore outperforms them on accurately identifying lncRNAs and mRNAs, especially on partial-length transcripts from the human and mouse datasets. Furthermore, lncScore also performed well on transcripts from five other species (Zebrafish, Fly, C. elegans, Rat, and Sheep), using models trained on human and mouse datasets. To speed up the prediction, multithreading is implemented within lncScore, and it only took 2 minute to classify 64,756 transcripts and 54 seconds to train a new model with 21,000 transcripts with 12 threads, which is much faster than other tools.

# Documentation

Documentation for the software is available at http://lncscore.openbioinformatics.org.

# Author

lncScore is developed by Jian Zhao (zhao_doctor@hotmail.com). For questions and comments, please contact Jian or submit an issue on github.

# Reference

- Zhao J, Song X\*, Wang K\*. lncScore: alignment-free identification of lncRNA from assembled novel transcripts. *Submitted*, 2016

# License

[WGLab MIT License](http://wglab.mit-license.org).
