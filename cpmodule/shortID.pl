#!/usr/bin/perl

use warnings;
use strict;

my $input = shift @ARGV;
my $output = '>' . shift @ARGV;

open(ID,$input) || die "$!";
open(IDNEW,$output) || die "$!";

while(<ID>){

	if(/^>/){
		my @a = split(" ",$_);
		my @b = split(/\|/,$a[0]);
		print IDNEW $b[0];
		print IDNEW "\n";
	}else {
		print IDNEW $_;
	} 
}



close(ID);
close(IDNEW);

exit;
