#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util 'sum';

# returns the average for each bin that has >= 10 pairwise comparisons

my $minpairs = 10;
my $infile;
GetOptions(
	'minpairs=i' => \$minpairs,
	'infile=s'   => \$infile,
);

print "bin\tgenpergeo\tmeanlat\n";
open my $fh, '<', $infile or die $!;
my @header;
my $bin;
my ( @meanlat, @geodist, @gendist );
LINE: while(<$fh>) {
	chomp;
	if ( not @header ) {
		@header = split /\t/, $_;
		next LINE;
	}
	my @record = split /\t/, $_;	
	if ( not $bin ) {
		$bin = $record[0];
		@meanlat = ( $record[-1] );
		@geodist = ( $record[-2] );
		@gendist = ( $record[-3] );
		next LINE;
	}
	if ( $bin eq $record[0] ) {
		push @meanlat, $record[-1];
		push @geodist, $record[-2];
		push @gendist, $record[-3];
		next LINE;	
	}
	else {
		if ( @meanlat >= $minpairs ) {
			my $meanlat = sum(@meanlat)/scalar(@meanlat);
			my @genpergeo;
			for my $i ( 0 .. $#gendist ) {
				push @genpergeo, $gendist[$i] / $geodist[$i] if $geodist[$i];
			}
			my $genpergeo = sum(@genpergeo)/scalar(@genpergeo);
			print $bin, "\t", $genpergeo, "\t", $meanlat, "\n";
		}		
		$bin = $record[0];
		@meanlat = ( $record[-1] );
		@geodist = ( $record[-2] );
		@gendist = ( $record[-3] );
	}
	
}