#!/usr/bin/perl

use Getopt::Std;
use File::Basename;

my $program = basename $0;
my $version = "1.3";
my %opts; getopts('hi:o:n:',\%opts);

die("
Program for extracting the longest transcript for Trinity genes (also suitable for ENSEMBL genes)

Common usage:
		$program -i FILE [options]

Options:
    -i FILE   Input Trinity fasta file
    -o FILE   Output file name [auto: *.long.fas]
    -n        Output ID model [0]
              [0]: raw transcript ID;
              [1]: clean transcript ID;
              [2]: clean gene ID;
    -h        Print this help message

Example:
    $program -i Trinity.fasta
    
Version: $version (Author: Yan-Bo Sun).
") if !$opts{i} || $opts{h};

my $inFile = $opts{i};
my $outID = defined $opts{n}?$opts{n}:0;
my $outFile;
if (defined $opts{o}){
	$outFile = $opts{o};
}else{
	$inFile =~ /(.*)\..*$/;
	$outFile = $1.".longTrans.fas";
}

my (%hashID1,%hashID2,%hashSeq);
my ($geneID,$tranID,$rawID,$seq);

open handleIn,"$inFile";
my $line = <handleIn>; 
($geneID,$tranID,$rawID) = &getID($line);
while(my $line = <handleIn>){
	chomp $line;
	if($line =~ /\>/){
		if(exists $hashSeq{$geneID}){
			my $oldLen = length($hashSeq{$geneID});
			my $newLen = length($seq);
			if ($newLen > $oldLen){
				$hashSeq{$geneID} = $seq;
				$hashID1{$geneID} = $tranID;
				$hashID2{$geneID} = $rawID;
			}
		}else{
			$hashSeq{$geneID} = $seq;
			$hashID1{$geneID} = $tranID;
			$hashID2{$geneID} = $rawID;
		}
		($geneID,$tranID,$rawID) = &getID($line);
		$seq = "";
	}else{
		$seq .= "$line\n";
	}
}
if(exists $hashSeq{$geneID}){
	my $oldLen = length($hashSeq{$geneID});
	my $newLen = length($seq);
	if ($newLen > $oldLen){
		$hashSeq{$geneID} = $seq;
		$hashID1{$geneID} = $tranID;
		$hashID2{$geneID} = $rawID;
	}
}else{
	$hashSeq{$geneID} = $seq;
	$hashID1{$geneID} = $tranID;
	$hashID2{$geneID} = $rawID;
}
close handleIn;

my $finalGeneNumber;
open handleOut,">$outFile";
if($outID == 0){
	while(my ($geneID,$rawID) = each %hashID2){
		print handleOut ">$rawID\n$hashSeq{$geneID}";
		$finalGeneNumber++;
	}
}elsif($outID == 1){
	while(my ($geneID,$tranID) = each %hashID1){
		print handleOut ">$tranID\n$hashSeq{$geneID}";
		$finalGeneNumber++;
	}
}else{
	while(my ($geneID,$tranID) = each %hashID1){
		print handleOut ">$geneID\n$hashSeq{$geneID}";
		$finalGeneNumber++;
	}
}
close handleOut;
print STDERR "=" x 50,"\n";
print STDERR "[InFile]: $inFile\n";
print STDERR "[OutFile]: $outFile\n";
print STDERR "[SaveID]: $outID\n";
print STDERR "[Processed]: $finalGeneNumber genes\n";
###############################################################
sub getID{
	my ($geneID,$tranID,$rawID);
	my $input = shift;
	chomp $input;
	$input =~ />(.*)$/;	$rawID = $1;
	if($input =~ /gene\:/){# biomart output
		$input =~ />([A-Z0-9]+).*gene\:([A-Z0-9]+).*$/;
		($tranID,$geneID) = ($1,$2);
	}else{
		$input =~ /\>(\S+)/;
		$tranID = $1;
		if($tranID =~ /\|/){# code for dealing with ensembl biomart id
			my @cols = split(/\|/,$tranID);
			$geneID = $cols[0];
			$tranID = $cols[1];
		}else{
			my @cols = split(/_/,$tranID);
			if(scalar @cols == 5){
				$geneID = "$cols[0]\_$cols[1]\_$cols[2]\_$cols[3]"; #with new trinity transcript
			}else{
				$geneID = "$cols[0]\_$cols[1]";  #with old trinity transcript
			}
		}
	}
	return ($geneID,$tranID,$rawID);
}
__END__
>comp0_c0_seq1 len=425 path=[5:0-424]  old Trinity
>TRINITY_DN86_c0_g1_i1 len=544 path=[522:0-269 1067:270-293 793:294-543] [-1, 522, 1067, 793, -2] new Trinity
>ENSTTRG00000000019|ENSTTRT00000000019
>ENST00000390578.1 cds chromosome:GRCh38:14:105897957:105897987:-1 gene:ENSG00000211918.1 gene_biotype:IG_D_gene transcript_biotype:IG_D_gene gene_symbol:IGHD2-15
