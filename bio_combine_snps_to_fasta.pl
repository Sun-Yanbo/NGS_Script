#!/usr/bin/perl
################################################
#Update: 11:36 2019/4/2 ÐÇÆÚ¶þ
#    01: improve depth estimate by omiting the Indel lines;
#    02: add an option '-d'(save degenerate bases);
#    03: add an option '-Q'(SNP calling QUAL);
#
#Update: 6:57 2018-6-20
#    01: try most avoid redundant bases;
#    02: change missing base to gap (not N)
#
#Update: 11:43 2016/6/13
#    01: print inf of depth into LOG file;
#    02: only handle SNP even when INDEL exist
################################################

use File::Basename;
use Getopt::Std;

my %opts;
my $update = "2019/04/02";
my $softName = basename $0;
getopts('hi:o:Q:q:m:M:n:l:rd',\%opts);
die("
Program for combining SNPs from filtered VCF to fasta.
Common usage:
    $softName -i <VCFFile> -o <OutFasta> [options]
Options:
    -i FILE   Input file name (in VCF[.gz] format)
    -o FILE   Output file name (in Fasta format)
    -Q INT    QUAL threshold [40]
    -q INT    genotype quality threshold (10) [No check]
    -m INT    the min. depth allowed for a site [auto:2.5%]
    -M INT    the Max. depth allowed for a site [auto:97.5%]
    -n INT    uncertained samples allowed at a site [0]
    -l INT    sequence length per line in output [60]
    -r        save the reference seq into the output
    -d        output degenerate bases when heterogeneity exist [ALT]
    -h        Print this help message
Version: $update.
Author: Yan-Bo Sun.
") if !$opts{i} || !$opts{o} || $opts{h};

########################################
print STDERR "\n===>Parameter parsing...\n";
my $VCFFile = $opts{i};print STDERR "    Input VCFFile:\t$VCFFile\n";
my $OutFasta = $opts{o};print STDERR "    Output Fasta:\t$OutFasta\n";
my $SNPQual = defined $opts{Q}?$opts{Q}:40;print STDERR "    SNP Calling Quality Cutoff:\t$SNPQual\n";
my $qualityCutoff = defined $opts{q}?$opts{q}:0;print STDERR "    Genotype Quality Cutoff:\t$qualityCutoff\n";
my $uncertainNumCutoff = defined $opts{n}?$opts{n}:0;print STDERR "    Max.NO. uncertained base:\t$uncertainNumCutoff\n";
my $seqLen4line = defined $opts{l}?$opts{l}:60;print STDERR "    Out Seq Length per Line:\t$seqLen4line\n";
my $saveRef = defined $opts{r}?'Yes':'No';print STDERR "    Whether save reference seq:\t$saveRef\n";
my $degenerate = defined $opts{d}?"Yes":"No"; print STDERR "    Whether save degenerate bases:\t$degenerate\n";
my ($minDepth,$maxDepth) = &estimateDepthCutoff();
my ($index4GQ,%hashName) = &estimateIndex4GQ();
my %hashBase = ('AG'=>'R','GA'=>'R','CT'=>'Y','TC'=>'Y','AC'=>'M','CA'=>'M',
				'TG'=>'K','GT'=>'K','GC'=>'S','CG'=>'S','AT'=>'W','TA'=>'W');

########################################
&main();
########################################
sub estimateDepthCutoff{
	print STDERR "\n===>Estimating depth cutoffs...\n";
	my $time1 = time();
	my ($depthL,$depthU);
	if(defined $opts{m} and defined $opts{M}){
		$depthL = $opts{m};
		$depthU = $opts{M};
	}else{
		if($VCFFile =~ /vcf\.gz/){
			open handIN,"zcat $VCFFile | " or die $!;
		}else{
			open handIN,"$VCFFile" or die $!;
		}
		my @depthList;
		my $lineNum;
		while(<handIN>){
			chomp;
			$lineNum++;
			next if(/^#/);
			my @cols = split(/\t/);
			my $infor = $cols[7];
			next if $cols[7] =~ /INDEL\;/; # pass if the line record is indel;
			print STDERR "\r$lineNum    $cols[0]\t$cols[1]\t$cols[2]    ";
			$infor =~ /DP=([0-9]+)/;
			push @depthList,$1;
		}
		close handIN;
		my @depthList_2 = sort {$a <=> $b} @depthList;
		my $indexL = int(0.025*scalar(@depthList_2));
		my $indexU = int(0.975*scalar(@depthList_2));
		($depthL,$depthU) = ($depthList_2[$indexL],$depthList_2[$indexU]);
		$depthL = $opts{m} if (defined $opts{m});
		$depthU = $opts{M} if (defined $opts{M});
	}
	my $time2 = time();my $timeused = $time2-$time1;
	print STDERR "    time used $timeused s\n";
	print STDERR "    minDepth: maxDepth: $depthL - $depthU\n";
	return ($depthL,$depthU);
}

########################################
sub estimateIndex4GQ{
	print STDERR "\n===>Determining the indexes of samples...";
	my $time1 = time();
	if($VCFFile =~ /.*\.gz/){
		open handIN,"zcat $VCFFile | " or die $!;
	}else{
		open handIN,"$VCFFile" or die $!;
	}
	my ($GQindex,%hash);
	my $check = -1;
	while(<handIN>){
		chomp;
		next if ($line =~ /^##/);
		my @cols = split(/\t/);
		if(/^#/){
			my $maxIndex = scalar(@cols)-1;
			foreach my $i (9..$maxIndex){$hash{$cols[$i]} = $i;}
		}else{
			my @formats = split(/\:/,$cols[8]);
			my $maxIndex = scalar(@formats)-1;
			foreach my $i (0..$maxIndex){if($formats[$i] eq "GQ"){$GQindex = $i;$check = 1;last;}}
			last;
		}
	}
	close handIN;
	if($check < 0 and $qualityCutoff > 0){
		die "\nError! GQ value not found!\n";
	}
	my $time2 = time();my $timeused = $time2-$time1;
	print STDERR "    \(GQ: $GQindex\) time used $timeused s\n";
	foreach my $sample (sort {$hash{$a}<=>$hash{$b}} keys %hash){
		print STDERR "    $sample: $hash{$sample}\n";
	}
	return ($GQindex,%hash);
}
########################################
sub getBase4sample_1{
	my $refBase = shift;
	my $altBase = shift;
	my $genotype = shift;
	my $genotypeQuality = shift;
	
	my $base2add;
	if($genotype eq '0/0'){
		$base2add = $genotypeQuality < $qualityCutoff?'N':$refBase;
	}elsif($genotype eq '1/1'){
		$base2add = $genotypeQuality < $qualityCutoff?'N':$altBase;
	}elsif($genotype eq '0/1' or $genotype eq '1/0'){
		if($degenerate eq "Yes"){
			my $tmp = $refBase.$altBase;
			$base2add = $genotypeQuality < $qualityCutoff?'N':$hashBase{$tmp};
		}else{
			$base2add = $genotypeQuality < $qualityCutoff?'N':$altBase;
		}
	}else{
		$base2add = '-';
	}
	return $base2add;
}
########################################
sub getBase4sample_2{
	my $refBase = shift;
	my $altBase1 = shift;
	my $altBase2 = shift;
	my $genotype = shift;
	my $genotypeQuality = shift;
	
	my $base2add;
	if($genotype eq '0/0'){
		$base2add = $genotypeQuality < $qualityCutoff?"N":$refBase;
	}elsif($genotype eq '1/1'){
		$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase1;
	}elsif($genotype eq '2/2'){
		$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase2;
	}elsif($genotype eq '0/1' or $genotype eq '1/0'){
		if($degenerate eq "Yes"){
			my $tmp = $refBase.$altBase1;
			$base2add = $genotypeQuality < $qualityCutoff?"N":$hashBase{$tmp};
		}else{
			$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase1;
		}
	}elsif($genotype eq '0/2' or $genotype eq '2/0'){
		if($degenerate eq "Yes"){
			my $tmp = $refBase.$altBase2;
			$base2add = $genotypeQuality < $qualityCutoff?"N":$hashBase{$tmp};
		}else{
			$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase2;
		}
	}elsif($genotype eq '1/2' or $genotype eq '2/1'){
		my $tmp = $altBase1.$altBase2;
		$base2add = $genotypeQuality < $qualityCutoff?"N":$hashBase{$tmp};
	}else{
		$base2add = '-';
	}
	return $base2add;
}
########################################
sub getBase4sample_3{
	my $refBase = shift;
	my $altBase1 = shift;
	my $altBase2 = shift;
	my $altBase3 = shift;
	my $genotype = shift;
	my $genotypeQuality = shift;
	
	my $base2add;
	if($genotype eq '0/0'){
		$base2add = $genotypeQuality < $qualityCutoff?"N":$refBase;
	}elsif($genotype eq '1/1'){
		$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase1;
	}elsif($genotype eq '2/2'){
		$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase2;
	}elsif($genotype eq '3/3'){
		$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase3;
	}elsif($genotype eq '0/1' or $genotype eq '1/0'){
		if($degenerate eq "Yes"){
			my $tmp = $refBase.$altBase1;
			$base2add = $genotypeQuality < $qualityCutoff?"N":$hashBase{$tmp};
		}else{
			$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase1;
		}
	}elsif($genotype eq '0/2' or $genotype eq '2/0'){
		if($degenerate == "Yes"){
			my $tmp = $refBase.$altBase2;
			$base2add = $genotypeQuality < $qualityCutoff?"N":$hashBase{$tmp};
		}else{
			$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase2;
		}
	}elsif($genotype eq '0/3' or $genotype eq '3/0'){
		if($degenerate eq "Yes"){
			my $tmp = $refBase.$altBase3;
			$base2add = $genotypeQuality < $qualityCutoff?"N":$hashBase{$tmp};
		}else{
			$base2add = $genotypeQuality < $qualityCutoff?"N":$altBase3;
		}
	}elsif($genotype eq '1/2' or $genotype eq '2/1'){
		my $tmp = $altBase1.$altBase2;$base2add = $genotypeQuality < $qualityCutoff?"N":$hashBase{$tmp};
	}elsif($genotype eq '1/3' or $genotype eq '3/1'){
		my $tmp = $altBase1.$altBase3;$base2add = $genotypeQuality < $qualityCutoff?"N":$hashBase{$tmp};
	}elsif($genotype eq '2/3' or $genotype eq '3/2'){
		my $tmp = $altBase2.$altBase3;$base2add = $genotypeQuality < $qualityCutoff?"N":$hashBase{$tmp};
	}else{
		$base2add = '-';
	}
	return $base2add;
}
########################################
sub printLOG{
	my $depth = shift;
	my $num4uncBase = shift;
	my $scaffold = shift;
	my $position = shift;
	
	my $whether2continue = 1;
	if($depth > $maxDepth){
		print LOG "$scaffold\t$position\tdepth \($depth\) too large!\n";
		$whether2continue = -1;
	}elsif($depth < $minDepth){
		print LOG "$scaffold\t$position\tdepth \($depth\) too small!\n";
		$whether2continue = -1;
	}elsif($num4uncBase > $uncertainNumCutoff){
		print LOG "$scaffold\t$position\ttoo many uncertained samples!\n";
		$whether2continue = -1;
	}else{
		print LOG "$scaffold\t$position\tpass\n";
	}
	return $whether2continue;
}
########################################
sub newLineSeq{
	my $rawSeq = shift;
	my $seqLen = length($rawSeq);
	my $outSeq;
	for(my $i=0;$i<$seqLen;$i+=$seqLen4line){
		$outSeq .= substr($rawSeq,$i,$seqLen4line);
		$outSeq .= "\n";
	}
	return $outSeq;
}
########################################
sub main{
	print STDERR "\n===>Extracting SNPs from each sample...\n";
	
	open SNP,">$OutFasta\.snp";
	print SNP "#Chrom\tPos\tRef";
	foreach my $sample (sort {$hashName{$a}<=>$hashName{$b}} keys %hashName){
		print SNP "\t$sample";
	}
	print SNP "\n";
	
	my $hashSeqs;
	if($VCFFile =~ /vcf\.gz/){
		open handIN,"zcat $VCFFile | " or die $!;
	}else{
		open handIN,"$VCFFile" or die $!;
	}
	open LOG,">$OutFasta\.inf";
	while(<handIN>){
		chomp;
		next if (/^#/);
		my @cols = split(/\t/);
		next if $cols[5] < $SNPQual;
		next if $cols[7] =~ /INDEL/; # pass if the line record is indel;
		
		my ($depth,$num4uncBase);
		$cols[7] =~ /DP=([0-9]+)/;$depth = $1;
		while(my($sample,$column) = each %hashName){
			my @formats = split(/\:/,$cols[$column]);
			my $curGQ = $formats[$index4GQ];
			$num4uncBase++ if ($curGQ < $qualityCutoff);
		}
		my $whether2continue = &printLOG($depth,$num4uncBase,$cols[0],$cols[1]);
		next if ($whether2continue < 0);
		
		my ($refBase,$altBase) = ($cols[3],$cols[4]);
		my @altBases = split(/\,/,$altBase);
		my $num4altBase = scalar(@altBases);
		$hashSeqs{'refSeq'}.="$refBase";
		
		$cols[0]=~/.*?(\d+).*/; #only report the scaffold number
		print SNP "$1\t$cols[1]\t$refBase";
		print STDERR "\r    $cols[0]    $cols[1]    ";
		foreach my $sample (sort {$hashName{$a}<=>$hashName{$b}} keys %hashName){
			my $base2add;
			my $column = $hashName{$sample};
			my @formats = split(/\:/,$cols[$column]);
			my ($genotype,$curGQ) = ($formats[0],$formats[$index4GQ]);
			if($num4altBase == 1){
				$base2add = &getBase4sample_1($refBase,$altBase,$genotype,$curGQ);
			}elsif($num4altBase == 2){
				my ($altBase1,$altBase2) = ($altBases[0],$altBases[1]);
				$base2add = &getBase4sample_2($refBase,$altBase1,$altBase2,$genotype,$curGQ);
			}elsif($num4altBase == 3){
				my ($altBase1,$altBase2,$altBase3) = ($altBases[0],$altBases[1],$altBases[2]);
				$base2add = &getBase4sample_3($refBase,$altBase1,$altBase2,$altBase3,$genotype,$curGQ);
			}else{
				$base2add = "N";
			}
			$hashSeqs{$sample}.="$base2add";
			print STDERR "$base2add    ";
			print SNP "\t$base2add";
		}
		print SNP "\n";
	}
	close LOG;
	close SNP;
	close handIN;
	
	print STDERR "\n\n===>Outputing the results...\n";
	`touch $OutFasta`;
	open handOUT,">$OutFasta" or die $!;
	while(my ($sample,$seq)=each %hashSeqs){
		my $seqLen = length($seq);
		next if ($sample eq 'refSeq' and $saveRef eq 'No');
		$sample =~ s/[\/\:]/./g;
		my $outSeq = &newLineSeq($seq);
		print handOUT ">$sample\n$outSeq";
		print STDERR "    >$sample: \t$seqLen bp\n";
	}
	close handOUT;
	print STDERR "\n===>Done!\n";
}
