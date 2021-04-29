#!/usr/bin/perl

#
# - iterate through numEss, SD, etc
# - get top X interactions (5000?)
# - calc LLS and record
#

$correlations_file = $ARGV[0];		# sorted corrs list
$output_file = $ARGV[1];
$MAXNUM = 200000;
if(defined $ARGV[2]){
	$MAXNUM = $ARGV[2];
}
#$sig_file = $ARGV[1];

=begin
open(REF, $gene_stats_file) || die "fail gene stats\n";		# gene-mean-std... file
while(<REF>) {
	chomp;
	($g, $m, $s, $obs, $ess) = split(/\t/);
	$stats{$g}{'mean'} = $m;
	$stats{$g}{'std'} = $s;
	$stats{$g}{'obs'} = $obs;
	$stats{$g}{'ess'} = $ess;
}
close(REF);
=cut

#open(REF, "/home/traver/Analysis/Reference/MolSigDB/c2.curated.kegg.v4.0.symbols.gmt") || die "Can't open specified reference file\n";
#open(REF, "kegg_pathway_Jul2018.rm_global.symbol.combined.under200") || die "Can't open specified reference file\n";
open(REF, "reference_gobp_experimental_rm_oxphos.gmt") || die "Can't open specified reference file\n";
#open(REF, $ARGV[1]) || die "can't open sig file\n";
# this is a Molecular Signatures Database file, one pathway/gene set per line
# ignore terms with > 50 genes (previous 200)
while(<REF>) {
        chomp;
        @data = split(/\t/);
        $name = $data[0];
        $numGenes = $#data -1;
        next if ($numGenes > 50);
        foreach $i (2..$#data) {
                $prot{$data[$i]} = 1;
        }
        foreach $i (2..$#data-1) {
                $g1 = $data[$i];
                foreach $j ($i+1..$#data) {
                        $g2 = $data[$j];
                        $int{$g1}{$g2} = 1;
                        $int{$g2}{$g1} = 1;
                }
        }
        $pathway{$name} = $numGenes;
}
close(REF);
@prots = keys %prot;
$num_prots_in_ref = $#prots + 1;
print "$num_prots_in_ref\n";
$num_ints = 0;
foreach $g1 ( keys %int ) {
  foreach $g2 ( keys %{ $int{$g1} } ) {
    $pos_ints++;
  }
}
$pos_ints = $pos_ints / 2;
$neg_ints = ($num_prots_in_ref * ($num_prots_in_ref-1))/2  - $pos_ints;
print "POS / NEG: $pos_ints // $neg_ints\n";
$binSize   = 1000;
#
# walk through through numEssObs cutoffs
# walk up sd cutoffs
# read top 25k pairs (that meet cutoff reqs)
# calc LLS on pairs
# write to file
#

my %generecall;

$fout_name = $output_file;
open(OUT, ">$fout_name");
$count = 0;
$sumPos = 0;
$sumNeg = 0;
print "Opening corrs file\n";
open(IN, $correlations_file) || die "fail opening $correlations_file\n";
while(<IN>) {
	chomp;
	($g1,$g2, $c, $p) = split(/\t/);

	$count++;
	$generecall{$g1} = 1;
	$generecall{$g2} = 1;

	#print "$count\n";
	if ( $prot{$g1} && $prot{$g2} ) {
		if ($int{$g1}{$g2} ) {
			$sumPos++;
			#print "$g1\t$g2\tPOS\n";
		} else {
			$sumNeg++;
			#print "$g1\t$g2\tNEG\n";
		}
	}
 	if ($count % $binSize == 0) {
		if ($sumPos && $sumNeg) {
			$obs = $sumPos/$sumNeg;
			$llr = (log( $obs ) - log( $pos_ints / $neg_ints)) / log(2);    # in ref set
		} else {
			$llr = 0;
		}
		$obs =  int( $obs * 1000 ) / 1000;
		$llr =  int( $llr * 1000 ) / 1000;
		$syscov = keys %generecall;
 		print OUT "$count\t$llr\t$obs\t$sumPos\t$sumNeg\t$syscov\n";
		#$sumPos = 0;  # set to zero for local, not cumulative, LLS calcs
		#$sumNeg = 0;
	}
	last if $count >= $MAXNUM;
}
if($count !=$MAXNUM && $count % $binSize != 0)  # added by Eiru
{
	if ($sumPos && $sumNeg) {
		$obs = $sumPos/$sumNeg;
		$llr = (log( $obs ) - log( $pos_ints / $neg_ints)) / log(2);    # in ref set
	} else {
		$llr = 0;
	}
	$obs =  int( $obs * 1000 ) / 1000;
	$llr =  int( $llr * 1000 ) / 1000;
	$syscov = keys %generecall;
 	print OUT "$count\t$llr\t$obs\t$sumPos\t$sumNeg\t$syscov\n";
}
close(IN);
close(OUT);
