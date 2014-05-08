#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util 'sum';
use Scalar::Util 'looks_like_number';

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
		warn "read header\n";
		next LINE;
	}
	my @record = split /\t/, $_;	
	if ( not $bin ) {
		$bin = $record[0];
		warn "started the first bin: $bin\n";
		next LINE if not looks_like_number $record[-2];
		@meanlat = ( $record[-1] );
		@geodist = ( $record[-2] );
		@gendist = ( $record[-3] );
		next LINE;
	}
	if ( $bin eq $record[0] ) {
		warn "extending bin $bin\n";
		next LINE if not looks_like_number $record[-2];		
		push @meanlat, $record[-1];
		push @geodist, $record[-2];
		push @gendist, $record[-3];
		next LINE;	
	}
	else {
		if ( @meanlat >= $minpairs ) {
			print_out();
		}
		else {
			warn "too few comparisons in $bin, skipping\n";
		}		
		$bin = $record[0];
		next LINE if not looks_like_number $record[-2];		
		@meanlat = ( $record[-1] );
		@geodist = ( $record[-2] );
		@gendist = ( $record[-3] );
		warn "started a new bin: $bin\n";
	}	
}

print_out();

sub print_out {
	my $meanlat = sum(@meanlat)/scalar(@meanlat);
	my @genpergeo;
	for my $i ( 0 .. $#gendist ) {
		push @genpergeo, $gendist[$i] / $geodist[$i] if $geodist[$i];
	}
	my $genpergeo = sum(@genpergeo)/scalar(@genpergeo);
	print $bin, "\t", $genpergeo, "\t", $meanlat, "\n";
}