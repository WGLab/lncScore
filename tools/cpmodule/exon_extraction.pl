#!usr/bin/perl
use warnings;
use strict;

my $input = shift @ARGV;
my $origin = shift @ARGV;
my $output = shift @ARGV;
$output = ">" . $output;
my %id = ();
my %length = ();
my %exonlen = ();

open(ID,$input) || die "$!";
while(<ID>){
	if (/^>/){
		/^>(\S*)/;
		#print $1."\n";
	#	my @a = split(" ",$_);
#		print $1,"\n";
		if (!exists($id{$1})) {
			$id{$1} = 0;
			$exonlen{$1} = '';
		}else{
			print "error!!\n";
			exit;
		}
	}
}
close(ID);

open(ID,$origin) || die "$!";
open(IDNEW,$output) || die "$!";
while(<ID>) {
	if (/\sexon\s*.*transcript_id\s\"(\S*)\";/) {

		if (exists($id{$1})) {
            #print $1."exon \n";
#			$id{$1} += 1;
			my @a = split(" ",$_);
			my $len = $a[4] - $a[3] +1;
			my $s = $exonlen{$1};
			$s = $s . " " . $len;
			$exonlen{$1} = $s;
			#$exonlen{$1} .= chr($len);

		}

	}
}
close(ID);
# print IDNEW "ID\tExonnumber\n";
my $flag = 0;
open(ID,$input) || die "$!";
while(<ID>){
        if (/^>/){
                /^>(\S*)/;
		my $len = join(" ",$exonlen{$1});
		my @tmp = split(' ',$exonlen{$1});
		my $num = @tmp;
		if ($num){
			$flag = 1;
			print IDNEW ">";
			print IDNEW $1,"\t";
			print IDNEW $exonlen{$1},"\n";
#		print IDNEW $id{$1},"\n";
		}
	} else {
		if ($flag) {
			print IDNEW $_;
			$flag = 0;
		}
	}
}
close(ID);
close(IDNEW);
exit;
