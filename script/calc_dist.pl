#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Phylo::Util::Logger ':levels';

# definition line format:
# >BOLD:AAA2367|KKCHE859-09|Lycosidae|Pardosa_lapponica|658|-94.186|58.78
# 1. taxon ID
# 2. sequence identifier
# 3. (optional) family
# 4. (optional) species
# 5. sequence length
# 6. longitude
# 7. latitude

# analysis steps
# 0. for all sequences in a bin, do a multiple sequence alignment
# 1. for each sequence pair within a bin, calculate the (uncorrected?) distance
# 2. for each sequence pair within a bin, calculate the geographical distance
# 3. divide genetic distance by geographical distance, e.g. substitutions/km


