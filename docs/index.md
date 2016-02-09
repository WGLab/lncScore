# lncScore

RNA-Seq based transcriptome assembly has been widely used in the identification of novel lncRNAs. However, it remains a challenge to identify lncRNAs from the assembled transcripts, particularly the partial-length ones, since partial-length protein-coding transcripts are more likely to be classified as lncRNAs due to their incomplete CDS. Furthermore, potential sequencing or assembly error that gain or abolish stop codons also complicates ORF-based prediction of lncRNAs. Here, we present a novel alignment-free tool, lncScore, which uses a logistic regression model with 11 carefully selected features derived from the open reading frame, exon, and the maximum coding subsequences. Compared to other state-of-the-art alignment-free tools (e.g. CPAT, CNCI, and PLEK), lncScore outperforms them on accurately identifying lncRNAs and mRNAs, especially on partial-length transcripts from the human and mouse datasets. Furthermore, lncScore also performed well on transcripts from five other species (Zebrafish, Fly, C. elegans, Rat, and Sheep), using models trained on human and mouse datasets. To speed up the prediction, multithreading is implemented within lncScore, and it only took 2 minute to classify 64,756 transcripts and 54 seconds to train a new model with 21,000 transcripts with 12 threads, which is much faster than other tools.

Please click the menu items to navigate through this website. Check [here](misc/whatsnew.md) to see what is new.

---

* ![new](img/new.png) 2016Feb09: Documentation for lncScore is added, and can now be accessed from http://lncscore.openbioinformatics.org.

* ![new](img/new.png) 2016Jan31: Initial version of lncScore (v1.0.0) is released to the public.

---

## Reference

- Zhao J, Song X\*, Wang K\*. lncScore: alignment-free identification of lncRNA from assembled novel transcripts. *Submitted*, 2016



