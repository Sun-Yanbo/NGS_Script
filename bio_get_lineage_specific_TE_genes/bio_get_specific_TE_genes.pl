#!/usr/bin/perl
use File::Find;
use File::Basename;
use Cwd qw(abs_path);
use FindBin qw($Bin);

my $Program = basename $0;
die "
A pipeline for identifying lineage-specific TE-like genes.
Version: 1.0.0      (Modified  at  2017-7-18)
Author:  Yan-Bo Sun (sunyanbo\@mail.kiz.ac.cn)

         Usage: $Program <config_file>

Software required: RepeatMasker (with HMMER installation)
Detailed information available in specTEG Configure file.

" if (@ARGV != 1);

################################################################################
##             0.parse the config file and check required programs
################################################################################
print STDERR "===>[0]::Parsing Configure File...\n";
my $key_geneFasFolder = "geneFasFolder";
my $key_geneType = "geneType";
my $key_genomeFolder = "genomeFolder";
my $key_querySpecies = "querySpecies";
my $key_repeatMaskerPath = "repeatMaskerPath";
#my $key_blastBinaryFolder = "blastBinaryFolder";
my $key_forceRunR = "forceRunR";
my $key_engine = "engine";
my $key_search_parallel = "search_parallel";
my $key_searchType = "searchType";
my $key_bitScore = "bitScore";
my $key_seqDiverge = "seqDiverge";
my $key_TE_library = "TE_library";
my $key_forceRunH = "forceRunH";
my $key_genetic_code = "genetic_code";
my $key_blast_parallel = "blast_parallel";
my $key_Evalue = "Evalue";
my $key_reference = "reference";
my $key_inflation = "inflation";
my $key_minCoverage = "minCoverage";
#my $key_flankDistance = "flankDistance";
my $key_minTELenth = "minTELenth";
my $key_outputFolder = "outputFolder";
my $key_outputFile = "outputFile";

my %hash_config;

#parse config_file
open IN,"$ARGV[0]";
while(<IN>){
	chomp;
	$_ =~ s/\s//g; #remove spaces in this line;
	next if ($_ eq "" or $_ =~ /^\#/);
	my $line;
	if(/\#/){
		$_ =~ /^(.*?)\#/;
		$line = $1;
	}else{
		$line = $_;
	}
	my @cols = split(/=/,$line);#print "@cols\n";
	if ((scalar @cols == 2) and ($cols[0] ne "") and ($cols[1] ne "")){
		$hash_config{$cols[0]} = $cols[1];
	}else{
		die "    [CONFIG]::Error! No value assigned for $cols[0] in Configure file\n";
	}
}
close IN;

#check repeatMaskerPath
my $repeatMaskerPath = $hash_config{$key_repeatMaskerPath};
if(-e $repeatMaskerPath){print STDERR "    [CONFIG]::RepeatMasker found at $repeatMaskerPath\n";}
else{die   "    [CONFIG]::Error! RepeatMasker not found at $repeatMaskerPath!\n";}

if(defined $hash_config{$key_TE_library}){
	my $repeatMaskerLibPath = $hash_config{$key_TE_library};
	if(-e $repeatMaskerLibPath){print STDERR "    [CONFIG]::RepeatMaskerLib found at $repeatMaskerLibPath\n";}
	else{print STDERR   "    [CONFIG]::Warnning! RepeatMaskerLib not found at $repeatMaskerPath!\n";}
}

system "chmod +x $Bin/tools/*.pl";
system "chmod +x $Bin/tools/bin/diamond-0.8.25/diamond";
system "chmod +x $Bin/tools/bin/mcl-14-137/src/shmcl/mcl";
system "chmod +x $Bin/tools/bin/hmmer-3.1b2/binaries/hmmscan";
system "chmod +x $Bin/tools/bin/hmmer-3.1b2/binaries/hmmpress";

#check genetic_code: used for orf prediction and translation
my $genetic_code = $hash_config{$key_genetic_code};
if ($genetic_code == 0){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'L','CTC' => 'L','CTA' => 'L','CTG' => 'L','ATT' => 'I','ATC' => 'I','ATA' => 'I','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => '*','TAG' => '*','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'K','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => '*','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => 'R','AGG' => 'R','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}elsif($genetic_code == 1){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'L','CTC' => 'L','CTA' => 'L','CTG' => 'L','ATT' => 'I','ATC' => 'I','ATA' => 'M','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => '*','TAG' => '*','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'K','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => 'W','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => '*','AGG' => '*','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}elsif($genetic_code == 2){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'T','CTC' => 'T','CTA' => 'T','CTG' => 'T','ATT' => 'I','ATC' => 'I','ATA' => 'M','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => '*','TAG' => '*','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'K','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => 'W','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => 'R','AGG' => 'R','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}elsif($genetic_code == 3){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'L','CTC' => 'L','CTA' => 'L','CTG' => 'L','ATT' => 'I','ATC' => 'I','ATA' => 'I','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => '*','TAG' => '*','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'K','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => 'W','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => 'R','AGG' => 'R','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}elsif($genetic_code == 4){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'L','CTC' => 'L','CTA' => 'L','CTG' => 'L','ATT' => 'I','ATC' => 'I','ATA' => 'M','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => '*','TAG' => '*','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'K','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => 'W','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => 'S','AGG' => 'S','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}elsif($genetic_code == 5){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'L','CTC' => 'L','CTA' => 'L','CTG' => 'L','ATT' => 'I','ATC' => 'I','ATA' => 'I','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => 'Q','TAG' => 'Q','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'K','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => '*','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => 'R','AGG' => 'R','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}elsif($genetic_code == 8){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'L','CTC' => 'L','CTA' => 'L','CTG' => 'L','ATT' => 'I','ATC' => 'I','ATA' => 'I','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => '*','TAG' => '*','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'N','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => 'W','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => 'S','AGG' => 'S','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}elsif($genetic_code == 9){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'L','CTC' => 'L','CTA' => 'L','CTG' => 'L','ATT' => 'I','ATC' => 'I','ATA' => 'I','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => '*','TAG' => '*','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'K','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => 'C','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => 'R','AGG' => 'R','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}elsif($genetic_code == 10){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'L','CTC' => 'L','CTA' => 'L','CTG' => 'L','ATT' => 'I','ATC' => 'I','ATA' => 'I','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => '*','TAG' => '*','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'K','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => '*','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => 'R','AGG' => 'R','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}elsif($genetic_code == 11){
	%hash_gen_code = (
	'TTT' => 'F','TTC' => 'F','TTA' => 'L','TTG' => 'L','CTT' => 'L','CTC' => 'L','CTA' => 'L','CTG' => 'S','ATT' => 'I','ATC' => 'I','ATA' => 'I','ATG' => 'M','GTT' => 'V',
	'GTC' => 'V','GTA' => 'V','GTG' => 'V','TCT' => 'S','TCC' => 'S','TCA' => 'S','TCG' => 'S','CCT' => 'P','CCC' => 'P','CCA' => 'P','CCG' => 'P','ACT' => 'T','ACC' => 'T',
	'ACA' => 'T','ACG' => 'T','GCT' => 'A','GCC' => 'A','GCA' => 'A','GCG' => 'A','TAT' => 'Y','TAC' => 'Y','TAA' => '*','TAG' => '*','CAT' => 'H','CAC' => 'H','CAA' => 'Q',
	'CAG' => 'Q','AAT' => 'N','AAC' => 'N','AAA' => 'K','AAG' => 'K','GAT' => 'D','GAC' => 'D','GAA' => 'E','GAG' => 'E','TGT' => 'C','TGC' => 'C','TGA' => '*','TGG' => 'W',
	'CGT' => 'R','CGC' => 'R','CGA' => 'R','CGG' => 'R','AGT' => 'S','AGC' => 'S','AGA' => 'R','AGG' => 'R','GGT' => 'G','GGC' => 'G','GGA' => 'G','GGG' => 'G','---' => '-');
}
################################################################################
##             0.DNA2AA
################################################################################
my $timeStart = time();
print STDERR "===>[0]::CDS translation...\n";
my ($out_main_specTEG,$out_PepFolder);
if (exists $hash_config{$key_geneFasFolder}){
	print STDERR "    [CONFIG]::Gene sequence folder = $hash_config{$key_geneFasFolder}\n";
	if(exists $hash_config{$key_outputFolder} and $hash_config{$key_outputFolder} ne "auto"){
		$out_main_specTEG = $hash_config{$key_outputFolder};
	}else{
		$out_main_specTEG = "$hash_config{$key_geneFasFolder}\_specTEG_main_outs_under_$ARGV[0]";
	}
	mkdir "$out_main_specTEG";
	$out_PepFolder = "pep_inputs";mkdir "$out_main_specTEG/$out_PepFolder";
	my @gene_files = glob "$hash_config{$key_geneFasFolder}/*.fa*";
	
	foreach (@gene_files){
		my $file_tmp = basename $_;
		my $file_inter = &sub_get_filename($file_tmp);
		unless(-e "$out_main_specTEG/$out_PepFolder/$file_inter\.fna"){
			system("cp $hash_config{$key_geneFasFolder}/$file_tmp $out_main_specTEG/$out_PepFolder/$file_inter\.fna");
		}
	}
	
	if ($hash_config{$key_geneType} == 2){ #2=cDNA seq;
		foreach (@gene_files){
			my $file_tmp = basename $_;
			my $file_inter = &sub_get_filename($file_tmp);
			unless(-e "$out_main_specTEG/$out_PepFolder/$file_inter\.faa"){
				&sub_getorf_file("$hash_config{$key_geneFasFolder}/$file_tmp","$out_main_specTEG/$out_PepFolder/$file_inter\.faa");
			}
		}
	}elsif($hash_config{$key_geneType} == 1){# 1=CDS;
		foreach (@gene_files){
			my $file_tmp = basename $_;
			my $file_inter = &sub_get_filename($file_tmp);
			unless(-e "$out_main_specTEG/$out_PepFolder/$file_inter\.faa"){
				&sub_translate_file("$hash_config{$key_geneFasFolder}/$file_tmp","$out_main_specTEG/$out_PepFolder/$file_inter\.faa");
			}
		}
	}else{
		die "    [CONFIG]::Error! No correct value assigned for 'geneType' in Configure file!\n";
	}
}else{
	die "    [CONFIG]::Error! No correct value assigned for 'geneFasFolder' in Configure file!\n";
}

`cp $ARGV[0] $out_main_specTEG/`;
chdir "$out_main_specTEG";
my $timeEnd = time();
my $timeUsed = $timeEnd - $timeStart;
my $timeCur = &getCurTime();
my $timeInfo = "[INF]::Done::Step-0::DNA translation    TimeUsed::$timeUsed Sec    $timeCur\n";
&printInf2file($timeInfo,"runtime.log");

################################################################################
##             1.Homologous genes prediction by get_homologues
################################################################################
$timeStart = time();
print STDERR "===>[1]::Homologues Identificaiton...\n";
print STDERR "    [CLUSTER]::clustering pep sequences...\n";
if($hash_config{$key_forceRunH} == 1){
	system "rm -fr $out_PepFolder\_homologues";
	if($hash_config{$key_reference} eq "NA"){
		system "$Bin/tools/get_homologues.pl -d $out_PepFolder -n $hash_config{$key_blast_parallel} -t 0 -M -X -C $hash_config{$key_minCoverage} -E $hash_config{$key_Evalue} -F $hash_config{$key_inflation} 1> 01.get_homologues.log";
	}else{
		my $reference4homolog = &sub_get_filename($hash_config{$key_reference});
		system "$Bin/tools/get_homologues.pl -d $out_PepFolder -n $hash_config{$key_blast_parallel} -t 0 -M -X -C $hash_config{$key_minCoverage} -E $hash_config{$key_Evalue} -F $hash_config{$key_inflation} -r $reference4homolog 1> 01.get_homologues.log";
	}
	system "touch $out_PepFolder\_homologues/get_homologues.run.ok";
}else{
	unless(-e "$out_PepFolder\_homologues/get_homologues.run.ok"){
		if($hash_config{$key_reference} eq "NA"){
			system "$Bin/tools/get_homologues.pl -d $out_PepFolder -n $hash_config{$key_blast_parallel} -t 0 -M -X -C $hash_config{$key_minCoverage} -E $hash_config{$key_Evalue} -F $hash_config{$key_inflation} 1> 01.get_homologues.log";
		}else{
			my $reference4homolog = &sub_get_filename($hash_config{$key_reference});
			system "$Bin/tools/get_homologues.pl -d $out_PepFolder -n $hash_config{$key_blast_parallel} -t 0 -M -X -C $hash_config{$key_minCoverage} -E $hash_config{$key_Evalue} -F $hash_config{$key_inflation} -r $reference4homolog 1> 01.get_homologues.log";
		}
		system "touch $out_PepFolder\_homologues/get_homologues.run.ok";
	}
}
$timeEnd = time();
$timeUsed = $timeEnd - $timeStart;
$timeCur = &getCurTime();
$timeInfo = "[INF]::Done::Step-1::Homologous Clustering    TimeUsed::$timeUsed Sec    $timeCur\n";
&printInf2file($timeInfo,"runtime.log");

################################################################################
##             2.TE Identificaiton by RepeatMasker
################################################################################
$timeStart = time();
print STDERR "===>[2]::TE Identificaiton by RepeatMasker...\n";
my $out_TEFolder = "TE_prediction"; mkdir "$out_TEFolder";
my @gene_files = glob "$out_PepFolder/*.fna";
foreach (@gene_files){
	chomp;
	my $file_tmp = basename $_;
	print STDERR "    [RepeatMasker]::$file_tmp ==> $file_tmp\.out\n";
	if($hash_config{$key_forceRunR} != 1 and (-e "$out_TEFolder/$file_tmp\.out")){
		next;
	}else{
		if($hash_config{$key_searchType} eq "slow"){
			if(defined $hash_config{$key_TE_library}){
				system "$hash_config{$key_repeatMaskerPath} -e $hash_config{$key_engine} -pa $hash_config{$key_search_parallel} -s -norna -nolow -div $hash_config{$key_seqDiverge} -lib $hash_config{$key_TE_library} -cutoff 250 -dir $out_TEFolder $_ 1> 02.RepeatMasker.log";
			}else{
				system "$hash_config{$key_repeatMaskerPath} -e $hash_config{$key_engine} -pa $hash_config{$key_search_parallel} -s -norna -nolow -div $hash_config{$key_seqDiverge} -dir $out_TEFolder $_ 1> 02.RepeatMasker.log";
			}
		}else{
			if(defined $hash_config{$key_TE_library}){
				system "$hash_config{$key_repeatMaskerPath} -e $hash_config{$key_engine} -pa $hash_config{$key_search_parallel} -norna -nolow -div $hash_config{$key_seqDiverge} -lib $hash_config{$key_TE_library} -cutoff 250 -dir $out_TEFolder $_ 1> 02.RepeatMasker.log";
			}else{
				system "$hash_config{$key_repeatMaskerPath} -e $hash_config{$key_engine} -pa $hash_config{$key_search_parallel} -norna -nolow -div $hash_config{$key_seqDiverge} -dir $out_TEFolder $_ 1> 02.RepeatMasker.log";
			}
		}
	}
}
$timeEnd = time();
$timeUsed = $timeEnd - $timeStart;
$timeCur = &getCurTime();
$timeInfo = "[INF]::Done::Step-2::TE-like Sequence Identificaiton    TimeUsed::$timeUsed Sec    $timeCur\n";
&printInf2file($timeInfo,"runtime.log");

################################################################################
##             3.Identificaiton of linege-specific TE domestication
################################################################################
$timeStart = time();
print STDERR "===>[3]::Parse Outputs of TE and Homologous...\n";
my @hash_TEout_list; #my $out_TEFolder = "TE_prediction";
foreach (@gene_files){
	chomp;
	my $file_tmp = basename $_;
	my $te_out_file = "$out_TEFolder/$file_tmp.out";
	print STDERR "    [READOUT]::reading RepeatMasker outfile $te_out_file\n";
	my %tmp_hash;
	open IN,$te_out_file;
	while(<IN>){
		chomp;
		my @cols = split(/\s+/);
		my $insertLen = scalar($cols[7] - $cols[6]);
		if ($insertLen >= $hash_config{$key_minTELenth} and $cols[1] >= $hash_config{$key_bitScore}){
			$tmp_hash{$cols[5]}.="$_";
		}
		
	}
	close IN;
	push @hash_TEout_list,\%tmp_hash;
}

#construct query flags which will be compared with final_real flags to determine whether lineage-specific;
print STDERR "    [VECTOR]::constructing compared vector  ";
@gene_files = glob "$out_PepFolder/*.fna";
my (@query_flags,@real_flags);
foreach (@gene_files){
	chomp;
	$target_sp = basename $_;

	my $tmp_flag = 0;
	push (@real_flags,$tmp_flag);
	
	my @query_sps = split(",",$hash_config{$key_querySpecies});
	foreach my $query_sp (@query_sps){
		chomp $query_sp;
		$query_sp = &sub_get_filename($query_sp);
		if("$query_sp\.fna" eq $target_sp){
			$tmp_flag = 1;last;
		}
	}
	push (@query_flags,$tmp_flag);
}
print STDERR "\[",join(" ",@query_flags),"\]\n";

my $out_PepFolder_homologues = "$out_PepFolder\_homologues";
my $out_PepFolder_homologues_list;
find(\&sub_find_cluster_list,$out_PepFolder_homologues);
my $out_PepFolder_homologues_fold = &sub_get_filename_2($out_PepFolder_homologues_list);
print STDERR "    [READCLUSTER]::List = $out_PepFolder_homologues_list\n";
print STDERR "    [READCLUSTER]::Folder = $out_PepFolder_homologues_fold\n";

open IN,"$out_PepFolder_homologues_list";
open OUT,">$hash_config{$key_outputFile}";

#print out title:
print OUT "ClusterID";
foreach (@gene_files){
	chomp;
	print OUT "\t",&sub_get_filename(basename $_);
}
print OUT "\tIfSpecific\tTEInfo\n";

#parse the first cluster to determine the real_flags
my ($file2read,@id_list);
my $line = <IN>; chomp $line;
$line =~ /cluster.*?dnafile: (.*)$/;
$file2read = "$out_PepFolder_homologues_fold/$1"; #print STDERR "\n$file2read\n";
while(<IN>){
	chomp;
	$line = $_;
	if($line =~ /cluster.*?dnafile: (.*)$/){
		#print STDERR "\n@id_list\n";
		my $parse_out = &sub_read_cluster($file2read);
		print OUT "$file2read\t",join("\t",@real_flags),"\t$parse_out\n"; #print out to $hash_config{$key_outputFile};
		@id_list = ();
		map($_=0,@real_flags); #change the memeber back to 0;
		$file2read = "$out_PepFolder_homologues_fold/$1";
	}else{
		$line =~ /: (.*)/;
		push @id_list,&sub_get_filename($1);
	}
}
my $parse_out = &sub_read_cluster($file2read);
print OUT "$file2read\t",join("\t",@real_flags),"\t$parse_out\n"; #print out to $hash_config{$key_outputFile};
close IN;
close OUT;

print STDERR "\n===>[4]::DONE!";
print STDERR "\n    [SAVE]::OutFile = ",abs_path("$hash_config{$key_outputFile}"),"\n";

$timeEnd = time();
$timeUsed = $timeEnd - $timeStart;
$timeCur = &getCurTime();
$timeInfo = "[INF]::Done::Step-3::Lineage-Specific Identificaiton    TimeUsed::$timeUsed Sec    $timeCur\n";
&printInf2file($timeInfo,"runtime.log");

################################################################################
##             5.Sub_functions
################################################################################
sub sub_find_cluster_list{
    if(/.*cluster_list$/){
        $out_PepFolder_homologues_list = $File::Find::name;
    }
}

sub sub_read_cluster{
	my $file_tmp = shift; #print "@id_list\n"; #print STDERR "\n$file_tmp\n"; 
	print STDERR "\r    [READCLUSTER]::Parsing ", basename($file_tmp),"    ";
	my $file_num;
	my $file_te_inf;
	open handle,"$file_tmp";
	while(<handle>){
		chomp;
		if(/>(\S+)/){
			$file_num++;
			my $gene_id = $1;
			my $file_name = $id_list[$file_num-1]; #print STDERR "    GeneID = $gene_id\tFileName = $file_name\n";
			my $file_index = -1;
			foreach my $i (0..scalar(@gene_files)-1){
				chomp $gene_files[$i]; #print STDERR &sub_get_filename(basename $gene_files[$i]),"\t$file_name\n";
				if(&sub_get_filename(basename $gene_files[$i]) eq $file_name){
					$file_index = $i;last;
				}
			}
			#print STDERR "$file_name\t$file_index\t$gene_id\t$hash_TEout_list[$file_index]->{$gene_id}\n";
			if (exists $hash_TEout_list[$file_index]->{$gene_id}){
				$real_flags[$file_index]=1;
				$file_te_inf .= "\t$hash_TEout_list[$file_index]->{$gene_id}";
			}
			
		}
	}
	close handle;
	my $tmp_query = join("\t",@query_flags);
	my $tmp_target = join("\t",@real_flags);
	if($tmp_query eq $tmp_target){
		print STDERR "Yes    \n";
		return "Yes$file_te_inf";
	}else{
		return "No";
	}
}

sub sub_getorf_file{ #get the longest ORF per cDNA as well as the pep seq
	my $file_cDNA = shift;
	my $file_pep = shift; print STDERR "    [GETORF]::$file_cDNA ==> $file_pep\n";
	my $seq;
	open hanlde,"$file_cDNA";
	open OUT,">$file_pep";
	while(<hanlde>){
		chomp;
		if(/\>/){
			if ($seq ne ""){
				print OUT &sub_getorf_long($seq),"\n";
				$seq = "";
			}
			print OUT "$_\n";
		}else{
			$seq .= $_;
		}
	}
	print OUT &sub_getorf_long($seq),"\n";
	close OUT;
	close hanlde;
}

sub sub_getorf_long{
	my $seq = shift;
	my @pepList;
	foreach my $i (1..3){
		my $dna = substr ($seq,$i-1);
		my $protein = &sub_translate_frame($dna);
		my @tmp = split (/\*/, $protein);
		foreach (@tmp){chomp;push (@pepList,$_) if (length($_) > 50);}
	}
	$seq_revCmp = &sub_revComp($seq);
	foreach my $i (1..3){
		my $dna = substr ($seq_revCmp,$i-1);
		my $protein = &sub_translate_frame($dna);
		my @tmp = split (/\*/, $protein);
		foreach (@tmp){chomp;push (@pepList,$_) if (length($_) > 50);}
	}
	
	my $maxPEP = $pepList[0];
	my $maxLength = length $pepList[0];
	foreach my $i (0..$#pepList){
		my $curLength = length($pepList[$i]);
		if($curLength > $maxLength){
			$maxLength = $curLength;
			$maxPEP = $pepList[$i];
		}
	}
	return $maxPEP;
}

sub sub_revComp{
	my $seq = shift; chomp $seq;
	my $seq_revCmp;
	foreach my $i (0..length($seq)-1){
		my $base = uc substr ($seq, $i, 1);
		my %hash_revCmp = (
		'A'=>'T','T'=>'A','G'=>'C','C'=>'G','R'=>'Y','Y'=>'R','W'=>'S','S'=>'W',
		'K'=>'M','M'=>'K','B'=>'V','V'=>'B','H'=>'D','D'=>'H','N'=>'N','-'=>'-');
		$seq_revCmp = exists $hash_revCmp{$base}?$hash_revCmp{$base}.$seq_revCmp:" ".$seq_revCmp;
	}
	return $seq_revCmp;
}

sub sub_translate_frame{
	my $seq = shift; chomp $seq;
	my $seq_len = length($seq);
	my $protein;
	for(my $i=0; $i <= ($seq_len - 3); $i+=3){
		my $codon = uc(substr($seq, $i, 3));
		$protein .= $hash_gen_code{$codon};
	}
	return $protein;
}

sub sub_get_filename{
	my $filename = shift;
	chomp $filename;
	$filename =~ /^(\S+?)\./;
	return $1;
}

sub sub_get_filename_2{
	my $filename = shift;
	chomp $filename;
	$filename =~ /^(\S+)\./;
	return $1;
}

sub sub_translate_file{
	my $file_cds = shift; #in
	my $file_pep = shift; #out
	
	print STDERR "    [CDS2PEP]::$file_cds\n";
	my ($seq_cds,$seq_pep);
	open IN,"$file_cds";
	open OUT,">$file_pep";
	while(<IN>){chomp; last if (/\>/);}
	print OUT "$_\n";
	while(<IN>){
		chomp;
		if(/>/){
			my $seq_len = length $seq_cds;
			for(my $i=0;$i<=($seq_len-3);$i+=3){
				my $codon = uc substr($seq_cds,$i,3);
				$seq_pep .= $hash_gen_code{$codon};
			}
			print OUT "$seq_pep\n";
			$seq_cds = "";$seq_pep = "";$seq_len = 0;
			print OUT "$_\n";
		}else{
			$seq_cds .= $_;
		}
	}
	for(my $i=0;$i<=($seq_len-3);$i+=3){
		my $codon = uc substr($seq_cds,$i,3);
		$seq_pep .= $hash_gen_code{$codon};
	}
	print OUT "$seq_pep\n";
	close OUT;
	close IN;
}

sub getCurTime{
	my($sec,$min,$hour,$day,$mon,$year)=localtime;
	$year-=100;
	$mon++; $mon=$mon<10?"0$mon":$mon;
	$day=$day<10?"0$day":$day;
	$hour=$hour<10?"0$hour":$hour;
	$min=$min<10?"0$min":$min;
	$sec=$sec<10?"0$sec":$sec;
	return "$hour:$min:$sec $mon/$day/$year";
}

sub printInf2file{
	my $information = shift;
	my $fileName = shift;
	open handle,">>$fileName";
	print handle "$information\n";
	close handle;
}