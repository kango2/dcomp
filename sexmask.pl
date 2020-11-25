#!/usr/bin/perl
use strict;
use warnings;

my %acc2chr = ();

open (R, "<$ARGV[0]") or die $!; ## _assembly_report.txt from RefSeq file has contigs/scaffolds listed against chromosomes in Assigned-Molecule column and RefSeq-Accn as sequence IDs in fasta file
while (<R>){
	chomp $_;
	next if (substr($_,0,1) eq "#");
	my @a = split ("\t", $_);
	$acc2chr{$a[6]} = $a[2];
}
close R;

open (F, "gunzip -c $ARGV[1] |") or die $!; ##input fasta for masking
open (O, "| gzip >$ARGV[2]") or die $!; ##output masked fasta
my $seqid = "";
my $seq = "";
while (my $line = <F>) {
	chomp $line;
	if ($line =~ />(\S+)/){
		$seqid = $1;
		print O "$line\n";
	}
	else {
		if ($acc2chr{$seqid} =~ /$ARGV[3]/){
			$line =~ s/\S/N/g;
			print O "$line\n";
		}
		else{
			print O "$line\n";
		}
	}
}
close F;
close O;

