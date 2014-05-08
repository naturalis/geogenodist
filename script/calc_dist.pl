#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Statistics::R;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Phylo::Matrices::Datum;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Phylo::Util::Logger ':levels';
use constant PI => 4 * atan2(1, 1);

# minimum bin size
my $binsize = 2;

# definition line format:
# >BOLD:AAA2367|KKCHE859-09|Lycosidae|Pardosa_lapponica|658|-94.186|58.78
my $tid = 0;
my $sid = 1;
my $fid = 2;
my $spp = 3;
my $len = 4;
my $lon = 5;
my $lat = 6;

# do the stupid muscle wrapper song and dance
@Bio::Tools::Run::Alignment::Muscle::MUSCLE_SWITCHES = grep { $_ ne 'profile' } @Bio::Tools::Run::Alignment::Muscle::MUSCLE_SWITCHES;

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new;
my $aln = Bio::Tools::Run::Alignment::Muscle->new( '-quiet' => 1 );
my $R = Statistics::R->new;

# return brief usage message
sub usage {
	"perl $0 -fasta <unaligned, sorted fasta> -b <bin size> [-verbose] > outfile.tsv";
}

# process command line arguments
sub get_args {
	my $verbosity = WARN;
	my $fasta;
	GetOptions(
		'verbose+'  => \$verbosity,
		'fasta=s'   => \$fasta,
		'binsize=i' => \$binsize,
	);	
	if ( not -e $fasta ) {
		die usage();
	}
	$log->VERBOSE( '-level' => $verbosity, '-class' => 'main' );
	return $fasta, $binsize;
}

# parse the BOLD definition lines
sub get_fields {
	my ( $defline, @indices ) = @_;
	my @parts = split /\|/, $defline;
	return @parts[@indices];
}

# move the cursor to the next bin
sub get_next_bin {
	my $seqio = shift;
	my @bin;
	if ( my $first = $seqio->next_seq ) {
		push @bin, $first;
	
		# get the bin id from the definition line
		my $id = $first->id();
		my $bin = get_fields($id,$tid);
	
		# advance the cursor until the bin id changes
		BIN: while( my $seq = $seqio->next_seq ) {
			my $binid = get_fields($seq->id(),$tid);
		
			# break out the loop if the bin id changes
			last BIN if $binid ne $bin;
		
			# grow the set
			push @bin, $seq;
		}
		if ( @bin >= $binsize ) {
			$log->info("going to align bin $bin");
			return $bin, $aln->align(\@bin);
		}
		else {
			$log->info("bin $bin has fewer than $binsize specimens, skipping");
			get_next_bin($seqio);
		}
	}
	return undef;
}

# calculate all pairwise distances within the alignment,
# returns a 2D array with id1, id2, distance
sub get_genetic_distances {
	my $alignment = shift;
	my @seq = $alignment->each_seq;
	my @result;
	
	# outer loop of pairwise comparisons
	for my $i ( 0 .. $#seq - 1 ) {
		my $s1 = Bio::Phylo::Matrices::Datum->new_from_bioperl( $seq[$i] );
		my $s1name = get_fields($s1->get_name,$sid);
		
		# inner loop of pairwise comparisons
		for my $j ( $i + 1 .. $#seq ) {
			my $s2 = Bio::Phylo::Matrices::Datum->new_from_bioperl( $seq[$j] );
			my $s2name = get_fields($s2->get_name,$sid);			
			my $dist = $s1->calc_distance($s2);	
			$log->info("calculating distance between $s1name and $s2name: $dist");			
			push @result, [ $s1name, $s2name, $dist ];
		}
	}
	return @result;
}

# calculate all great circle distances within the alignment,
# returns a 2D array with id1, id2, distance
sub get_geographical_distances {
	my $alignment = shift;
	my @seq = $alignment->each_seq;
	my @result;
	no warnings 'uninitialized';
	
	# outer loop of pairwise comparisons
	OUTER: for my $i ( 0 .. $#seq - 1 ) {
		my $s1 = Bio::Phylo::Matrices::Datum->new_from_bioperl( $seq[$i] );
		my $def1 = join ' ', $seq[$i]->id, $seq[$i]->desc;
		my $s1name = get_fields( $def1, $sid );
		my @s1xy = get_fields( $def1, $lon, $lat );
		
		# inner loop of pairwise comparisons
		INNER: for my $j ( $i + 1 .. $#seq ) {
			my $s2 = Bio::Phylo::Matrices::Datum->new_from_bioperl( $seq[$j] );
			my $def2 = join ' ', $seq[$j]->id, $seq[$j]->desc;
			my $s2name = get_fields( $def2, $sid );
			my @s2xy = get_fields( $def2, $lon, $lat );	
			$log->info("x1,y1=@s1xy x2,y2=@s2xy");								
			push @result, { 
				'seqid1'  => $s1name,
				's1lon'   => $s1xy[0],
				's1lat'   => $s1xy[1],
				'seqid2'  => $s2name, 
				's2lon'   => $s2xy[0],
				's2lat'   => $s2xy[1],				
				'geodist' => calc_dist(@s1xy,@s2xy),
				'meanlat' => ( ( $s1xy[1] + $s2xy[1] ) / 2 ),
			};
		}
	}
	return @result;

}

# calls R to calculate the great circle "Vincenty Ellipsoid" distance between two points. 
# requires the R package 'geosphere'. returns a single distance.
sub calc_dist {
	my @points = @_;
	my @r = @points;
	for my $i ( 0 .. $#r ) {
		if ( $r[$i] < -90 or $r[$i] > 90 ) {
			$log->warn("coordinate out of allowed range (-90 < value < 90): $r[$i]");
			return undef;
		}
	}
	my $command = <<"R";
library(geosphere);
xy1 <- rbind(c($r[0],$r[1]))
xy2 <- rbind(c($r[2],$r[3]))
result <- distm(xy1,xy2,fun=distHaversine)
R
	$log->debug("going to run R:\n$command");
	$R->run($command);
	my $result = $R->get('result');
	$log->info("distance is ".$result->[1]);
	return $result->[1] / 1000;
}

# analysis steps
# 1. for all sequences in a bin, do a multiple sequence alignment
# 2. for each sequence pair within a bin, calculate the (uncorrected?) distance
# 3. for each sequence pair within a bin, calculate the geographical distance
# 4. divide genetic distance by geographical distance, e.g. substitutions/km
# report: bin, seq id1, long1, lat1, seq id2, long2, lat2, genetic distance, geo distance, mean latitude
sub main {
	# process command line arguments
	my ( $fasta, $binsize ) = get_args();
	$log->info("going to process sequences from $fasta");
	$log->info("going to process all bins with more than $binsize specimens");
	
	# instantiate sequence reader
	my $seqio = Bio::SeqIO->new(
		'-format' => 'fasta',
		'-file'   => $fasta,
	);
	
	# print output header
	print join "\t", qw(bin sid1 long1 lat1 sid2 long2 lat2 gendist geodist meanlat);
	print "\n";
	
	# iterate over acceptable bins
	while( my ( $bin, $aln ) = get_next_bin($seqio) ) {
		last unless $bin and $aln;
		
		# calculate the genetic distances
		$log->info("going to calculate genetic distances in bin $bin");
		my @gendist = get_genetic_distances($aln);
		
		# calculate geographic distance
		$log->info("going to calculate geographic distances in bin $bin");
		my @geodist = get_geographical_distances($aln);
		
		# iterate over specimen pairs
		for my $i ( 0 .. $#gendist ) {
			my @record = (
				$bin,
				$gendist[$i]->[0],
				$geodist[$i]->{'s1lon'},
				$geodist[$i]->{'s1lat'},
				$gendist[$i]->[1],
				$geodist[$i]->{'s2lon'},
				$geodist[$i]->{'s2lat'},
				$gendist[$i]->[2],
				$geodist[$i]->{'geodist'},	
				$geodist[$i]->{'meanlat'},								
			);	
			print join "\t", @record;
			print "\n";	
		}	
	}
}

main();
