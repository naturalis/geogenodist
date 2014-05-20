#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util 'sum';
use Scalar::Util 'looks_like_number';
use Bio::Phylo::Util::Logger ':levels';

# returns the average for each bin that has >= 10 pairwise comparisons

# process command line arguments
my $verbosity = WARN;
my $minpairs = 10;
my $infile;
my $absolute;
GetOptions(
	'minpairs=i' => \$minpairs,
	'infile=s'   => \$infile,
	'verbose+'   => \$verbosity,
	'absolute'   => \$absolute,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# print the output header
print "bin\tgenpergeo\tmeanlat\n";

# compose the parameters to open the file
my $fh;
if ( $infile =~ /\.gz$/ ) {
	open $fh, '-|', "gunzip -c $infile" or die $!;
}
else {
	open $fh, '<', $infile or die $!;
}

# start reading the file
my @header;
my $bin;
my ( @meanlat, @geodist, @gendist );
LINE: while(<$fh>) {
	chomp;
	if ( not @header ) {
		@header = split /\t/, $_;
		$log->info("read header");
		next LINE;
	}
	my @record = split /\t/, $_;	
	if ( not $bin ) {
		$bin = $record[0];
		$log->info("started the first bin: $bin");
		next LINE if not looks_like_number $record[-2];
		@meanlat = ( $record[-1] );
		@geodist = ( $record[-2] );
		@gendist = ( $record[-3] );
		next LINE;
	}
	if ( $bin eq $record[0] ) {
		$log->info("extending bin $bin");
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
			$log->warn("too few comparisons in $bin, skipping");
		}		
		$bin = $record[0];
		next LINE if not looks_like_number $record[-2];		
		@meanlat = ( $record[-1] );
		@geodist = ( $record[-2] );
		@gendist = ( $record[-3] );
		$log->info("started a new bin: $bin");
	}	
}

print_out();

sub print_out {
	my $meanlat = sum(@meanlat)/scalar(@meanlat);
	$meanlat = abs($meanlat) if $absolute;
	my @genpergeo;
	for my $i ( 0 .. $#gendist ) {
		push @genpergeo, $gendist[$i] / $geodist[$i] if $geodist[$i];
	}
	if ( scalar @genpergeo ) {
		my $genpergeo = sum(@genpergeo)/scalar(@genpergeo);
		print $bin, "\t", $genpergeo, "\t", $meanlat, "\n";
	}
}