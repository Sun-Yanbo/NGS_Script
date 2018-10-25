#!/usr/bin/perl
use File::Basename;

my $program = basename $0;
die "
Estimating the index of tissue specificity based on multi-tissue expression matrix
Usage: $program <expMatrixFile> <ReplicateMap> <ExpCutoff[1]> <TSICutoff[0.95]>

Note:
<ExpMatrix> FILE (sep=\\t; header=True; GeneID\\tTissue1\\tTissue2...):
GeneID	Sample-1	Sample-2 ...
geneX	1.234	2.345
...
	
<ReplicateMap> FILE:
Sample-1	TissueX
Sample-2	TissueX
Sample-3	TissueY

<ExpCutoff> INT: 1 (for FPKM)

<TSICutoff> FLOAT: 0.95 (for detection of tissue-specific genes)
" unless (@ARGV == 4);

my $expMatrixFile = shift;
my $ReplicateMap = shift;
my $ExpCutoff = shift;
my $TSICutoff = shift;

my (%hash_map,%hash_group);
open handle,"$ReplicateMap";
while(<handle>){
	chomp;
	my @cols = split(/\t/);
	if (@cols >= 2){
		$hash_map{$cols[0]} = $cols[1];
		$hash_group{$cols[1]}++;
	}
}
close handle;

open handle,"$expMatrixFile";
my $title = <handle>;chomp $title;
my @tissues = split(/\t/,$title);shift @tissues;

#transform the raw matrix to mean exp matrix:
open OUT,">$expMatrixFile\.mean";
select OUT; print "GeneID";
my @groups = sort keys %hash_group;
foreach (@groups){chomp; print "\t$_";} print "\n";

while(<handle>){
	chomp;
	my @cols = split(/\t/);
	my $GeneID = shift @cols;
	
	my (%hash_mean_exp,%hash_group_rep);
	foreach my $i (0..$#cols){
		my $group = $hash_map{$tissues[$i]};
		my $expr = $cols[$i];
		if($expr ne "NA" and $expr ne ""){
			$hash_mean_exp{$group} += $expr;
			$hash_group_rep{$group}++;
		}
	}
	
	print "$GeneID";
	foreach (@groups){
		if(exists $hash_mean_exp{$_}){
			my $meanExp = $hash_mean_exp{$_}/$hash_group_rep{$_};
			print "\t$meanExp";
		}else{
			print "\tNA";
		}
	}
	print "\n";
}
close OUT;

open OUT,">$expMatrixFile\.mean\.TSI"; select OUT;
open handle,"$expMatrixFile\.mean";
print "GeneID\tTSI\tMax.Exp\tMax.Exp.Tissue\n";
while(<handle>){
	chomp;
	my @cols = split(/\t/);
	my $GeneID = shift @cols;
	
	@cols_sort = sort {$b<=>$a} (@cols);
	my $exp_max = $cols_sort[0];
	next if ($exp_max == 0);
	
	my ($geneTSI,$tissueNum);
	for my $expr (@cols){
		if($expr ne "NA" and $expr ne ""){
			if ($expr >= $ExpCutoff){
				$geneTSI += (1-$expr/$exp_max);
			}else{
				$geneTSI ++;
			}
			$tissueNum++;
		}
	}
	$geneTSI /= ($tissueNum-1);
	
	if($geneTSI > 1){ #all expression values are less than 1
		print "$GeneID\tNA\t$exp_max\tNA\n";
	}elsif($geneTSI >=$TSICutoff){
		my $idx_max;
		for my $i (0..$#cols){
			if($cols[$i] == $exp_max){
				$idx_max = $i;last;
			}
		}
		print "$GeneID\t$geneTSI\t$exp_max\t$groups[$idx_max]\n";
	}else{
		print "$GeneID\t$geneTSI\t$exp_max\tNA\n";
	}
}
close handle;
close OUT;
