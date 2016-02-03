lncScore is a python package for the identification of lncRNA from the assembled novel transcripts, and it also can be used to calculate the coding potential.
Version: 1.0.0. Last Modified: 2016-1-31
Authors: Jian Zhao (zhaojian@nuaa.edu.cn or zhao_doctor@hotmail.com)

Software Prerequisites:
      The following two software should be installed in your cluster or computer before running the lncScore.py
?         Perl (>=5.10.1), https://www.perl.org/get.html.
?         Python (>= 2.7), https://www.python.org/downloads/.
                 The scikit-learn module, http://scikit-learn.org/stable/install.html.

Usage:
Usage: lncScore.py [options]
    Options:
      -h, --help            show this help message and exit
      -f input files, --file=input files
                            enter transcripts in .bed or .fasta format: if this is
                            a .bed format file, '-r' must be specified; if this is
                            a .fasta format file, ignore the '-r'.
      -g gtf file name, --gtf=gtf file name
                            please enter your gtf files
      -o output files, --out=output files
                            assign your output file
      -p prallel numbers, --parallel=prallel numbers
                            please enter your specified speed ratio
      -x hexamer matrix, --hex=hexamer matrix
                            Prebuilt hexamer frequency table (Human, Mouse, Fly,
                            Zebrafish, C. elegans, Sheep and Rat). Run
                            'make_hexamer_tab.py' to make this table out of your
                            own training dataset.
      -t training dataset, --train=training dataset
                            Please enter your specified training dataset
      -r reference genome files, --ref=reference genome files
                            Reference genome sequences in FASTA format. Ignore
                            this option if sequences file was provided to '-f'.
                            Reference genome file will be indexed automatically
                            (produce *.fa file along with the original *.bed file
                            within the same directory) if hasn't been done.

Note:
     1. If the input is a FASTA format file, then make sure that the description line
      (begins with a ¡°>¡± symbol) only contains the transcript ID, for example:
      ¡°ENST00000379410.7¡±, which should be found in the GTF format file.
     2. If the input file in .bed format, then an additional python package named
        'pysam' was required to be installed first.   

Example:
      python lncScore.py -f test/human_test.fa -g test/human_test.gtf -o result -p 1 -x dat/Human_Hexamer.tsv -t dat/Human_training.dat

Output file format:
      Columns (from left to right):
          1. Transcript id;
          2. Index;
          3. Coding score.
      For example:
          Transcript_id   Index   Coding_score
          ENST00000342066.7      coding  0.999999999941
          ENST00000327044.6      coding  1.0
          ENST00000338591.7      coding  0.999999999998
          ENST00000562538.1      noncoding       0.0950315486986
          ENST00000415166.1      noncoding       0.0264167916398
          ENST00000565955.1      noncoding       0.0149329291793
