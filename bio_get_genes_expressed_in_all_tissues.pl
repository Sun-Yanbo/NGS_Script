#!/usr/bin/perl
use File::Basename;
my $Program = basename $0;
die "
Program for getting expressed genes with expression matrix of multi-tissues
$Program <ExpMatrix> <ReplicateMap> <ExpCutoff>

Note: 
<ExpMatrix> FILE:
GeneID	Sample-1	Sample-2 ...

<ReplicateMap> FILE:
Sample-1	TissueX
Sample-2	TissueX
Sample-3	TissueY
...

<ExpCutoff> INT: 1 (for FPKM)
" unless (@ARGV == 3);

my $ExpMatrix = shift;
my $ReplicateMap = shift;
my $ExpCutoff  = shift;

my (%hash_tissue,%hash_group);
open handle,"$ReplicateMap";
while(<handle>){
	chomp;
	my @cols = split(/\t/);
	next if(@cols < 2);
	$hash_tissue{$cols[0]} = $cols[1];
	$hash_group{$cols[1]}++;
}
close handle;

open handle,"$ExpMatrix";
$_ = <handle>; chomp;
my @tissues = split(/\t/);
shift @tissues;

open OUT,">$ExpMatrix\.AllTissueExp\.txt";
print OUT "GeneID";
for my $group (sort keys %hash_group){
	chomp $group;
	print OUT "\t$group";
}
print OUT "\n";

my $ExpGeneNum;
while(<handle>){
	chomp;
	next unless defined;
	my @cols = split(/\t/);
	my $GeneID = shift @cols;
	my (%hash_exp,%hash_rep,%hash_flag);
	for my $group (sort keys %hash_group){
		chomp $group;
		$hash_flag{$group}=0;
	}
	foreach my $i (0..$#cols){
		my $tissue = $tissues[$i];
		my $group = $hash_tissue{$tissue};
		next if($cols[$i] eq "NA" or $cols[$i] eq "");
		$hash_flag{$group} = 1 if (scalar ($cols[$i]) >= $ExpCutoff);
		$hash_exp{$group} += scalar($cols[$i]);
		$hash_rep{$group} ++;
	}
	print STDERR "[STDERR]: $GeneID\t",values %hash_flag,"\n";
	my $tmpFlag = 1;
	foreach (values %hash_flag){
		chomp;if ($_ == 0){$tmpFlag=0; last;}
	}
	if ($tmpFlag == 1){
		$ExpGeneNum++;
		print OUT "$GeneID";
		foreach my $tissue (sort keys %hash_flag){
			chomp $tissue;
			my $meanExp;
			if($hash_rep{$tissue} == 0){
				$meanExp = "NA";
			}else{
				$meanExp = $hash_exp{$tissue}/$hash_rep{$tissue};
			}
			print OUT "\t$meanExp";
		}
		print OUT "\n";
	}
}
print "The totle number of genes expressed in all tissues is: $ExpGeneNum\n";
print "The mean expression values can be found in output: $ExpMatrix\.AllTissueExp\.txt\n";
close handle;
close OUT;
