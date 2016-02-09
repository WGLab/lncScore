# General usage:

```
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
```

Note:

1. If the input is a FASTA format file, then make sure that the description line
      (begins with a `>` symbol) only contains the transcript ID, for example:
      `ENST00000379410.7`, which should be found in the GTF format file.
2. If the input file in .bed format, then an additional python package named
        'pysam' is required to be installed first.   

Example:

      `python lncScore.py -f test/human_test.fa -g test/human_test.gtf -o result -p 1 -x dat/Human_Hexamer.tsv -t dat/Human_training.dat`

# Output file format:

Columns (from left to right) are 1. Transcript id; 2. Index; 3. Coding score. For example:

```
          Transcript_id   Index   Coding_score
          ENST00000342066.7      coding  0.999999999941
          ENST00000327044.6      coding  1.0
          ENST00000338591.7      coding  0.999999999998
          ENST00000562538.1      noncoding       0.0950315486986
          ENST00000415166.1      noncoding       0.0264167916398
          ENST00000565955.1      noncoding       0.0149329291793
```

<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-73542276-1', 'auto');
  ga('send', 'pageview');

</script>
