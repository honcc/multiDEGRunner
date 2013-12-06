#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Math::Random;
use Math::Round;
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs abs2rel);
use List::Util qw (sum shuffle min max);
use Scalar::Util qw(looks_like_number);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This perl script to read the path of all bam files of multiple sample with multiple replicates and pipe the bams to cuffDiff and then convert the count output of cuffDiff to the input of DESeq and then run DESeq;
#	The script will launch cuffDiff in parallel and will start running DESeq when the corresponding cuffDiff output is finished;
#	
#	Options	
#		--sampleInfoPath=			path of the bam files, HTSeq Count output, as well as the sample names;
#		--refName=					the name of the reference sample; all comparison in Cuffdiff and DESeq will be paired with refSample as the reference sample; use 'all' to run all; multiple input alllowed;
#		--geneNameLenPath=			path of a txt file contains the genename and gene accession number;
#		--countType=				string [express_eff_count]; could be twoColumn, express_eff_count, express_est_count
#		--outDir=					output Dir;
#
#	History:
#		v0.1 
#			-built on v0.6 of multiCuffDiffDESeqRunner, but replacing CuffDiff with DESeq2;
#			-removed different choices of DESeq1 modes in options;
#
#		v0.2
#			[Thu 12 Sep 2013 10:49:34 CEST] cleaned;
#			[07/11/2013 15:31] will print senseAntisenseSameRow result table by recognizing the geneID;
#
#		v0.3
#			[08/11/2013 14:35] will be compatible with samples without replicates;
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-11-08 21:07]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/differentialExpression/multiDEGRunner/v0.3/multiDEGRunner_v0.3.pl --sampleInfoPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/differentialExpression/multiDEGRunner/v0.2/PF_4timePt_1rep_withAntisense.txt --refName=all --geneNameLenPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/plasmoDB/PF3D7/GFFSeqExtractorWithAntisense/log/mRNA.info.tsv --countType=twoColumn --outDir=/Volumes/F_Analysis/NGS/results/plasmodium/differentialExpression/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/differentialExpression/multiDEGRunner/v0.3/multiDEGRunner_v0.3.pl
#	--sampleInfoPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/differentialExpression/multiDEGRunner/v0.2/PF_4timePt_1rep_withAntisense.txt
#	--refName=all
#	--geneNameLenPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/plasmoDB/PF3D7/GFFSeqExtractorWithAntisense/log/mRNA.info.tsv
#	--countType=twoColumn
#	--outDir=/Volumes/F_Analysis/NGS/results/plasmodium/differentialExpression/BCV0.7/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $ARGVStr = join "\n", (&currentTime(), abs_path($0), @ARGV);#->355
my $globalScriptDirPath = dirname(rel2abs($0));
open DEBUGLOG, ">", "$globalScriptDirPath/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1918, readParameters|2187
#	secondaryDependOnSub: currentTime|355
#
#<section ID="startingTasks" num="0">
#---1 read the parameters
&printCMDLogOrFinishMessage("CMDLog");#->1918
my ($sampleInfoPath, $refNameAry_ref, $geneNameLenPath, $countType, $outDir) = &readParameters();#->2187
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my $DESeq1FDR = 0.05;
my $DESeq2FDR = 0.05;
my $EdgeRFDR = 0.05;
my $BCV = 0.7; #----dispersion value for EdgeR when there's no replicate, between suggest to try from 0.2 to 0.5
my $FDRTag = "D1_$DESeq1FDR\_D2_$DESeq2FDR\_Er_$EdgeRFDR";
my $minCountPerReplicate = -1;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_processInputData
#	primaryDependOnSub: checkSampleFiles|294, getAllCountInfo|1056, getGeneNameAndLength|1183
#	secondaryDependOnSub: reportStatus|2223, userDecisionToProceed|2467
#
#<section ID="processInputData" num="2">
#----read the BAM info
my ($refNameHsh_ref, $sampleInfoHsh_ref) = &checkSampleFiles($sampleInfoPath, $refNameAry_ref);#->294
my ($geneNameLenHsh_ref) = &getGeneNameAndLength($geneNameLenPath);#->1183
my ($allSampleCountInfoHsh_ref) = &getAllCountInfo($sampleInfoHsh_ref, $countType, $geneNameLenHsh_ref, $minCountPerReplicate);#->1056
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="3">
my $dirToPrintSubHTML = "$outDir/html";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_runDEAnalyses
#	primaryDependOnSub: runAllComparisons|2244
#	secondaryDependOnSub: generateComparisonCombination|469, generateDESeq1RScript|516, generateDESeq2RScript|698, generateEdgeRRScript|776, launchAndMonitorAllCMD|1219, overlapGeneExpDiffResults|1398, printAllComparisonFoldChange|1770, printCombineResults|1951
#
#<section ID="runDEAnalyses" num="4">
my ($allComparisonFoldChangeHsh_ref, $allResultFilePathHsh_ref, $allComparisonNameHsh_ref, $overallCountFilePathHsh_ref) = &runAllComparisons($sampleInfoPath, $outDir, $allSampleCountInfoHsh_ref, $DESeq1FDR, $DESeq2FDR, $EdgeRFDR, $FDRTag, $refNameHsh_ref, $geneNameLenHsh_ref, $sampleInfoHsh_ref, $countType, $BCV);#->2244
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_summarizeResults
#	primaryDependOnSub: printAllComparisonFoldChange|1770, printAllSampleCountInfo|1843
#	secondaryDependOnSub: >none
#
#<section ID="summarizeResults" num="5">
#---print All Comparison Fold Change
&printAllComparisonFoldChange($allComparisonFoldChangeHsh_ref, $geneNameLenHsh_ref, "$outDir/ovrlp/", 'allComparisons', $allComparisonNameHsh_ref, $FDRTag);#->1770
#----print allSampleCountResults
my ($allSampleCountPathHsh_ref, $sampleTotalFragNumHsh_ref, $sampleColumnAry_ref) = &printAllSampleCountInfo($allSampleCountInfoHsh_ref, $outDir);#->1843
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_PCAAndHeatmap
#	primaryDependOnSub: generateAllHeatmap|373, generateAllPCA|394
#	secondaryDependOnSub: PCAWithFactoMiner|220, generateAllSampleDistHeatmap|416
#
#<section ID="PCAAndHeatmap" num="6">
&generateAllHeatmap($outDir, $allSampleCountPathHsh_ref, $sampleColumnAry_ref, $allResultFilePathHsh_ref);#->373
&generateAllPCA($outDir, $allSampleCountInfoHsh_ref, $allSampleCountPathHsh_ref, $allResultFilePathHsh_ref);#->394
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_outputHTML
#	primaryDependOnSub: generateMasterHTML|891
#	secondaryDependOnSub: >none
#
#<section ID="outputHTML" num="7">
&generateMasterHTML($allResultFilePathHsh_ref, $dirToPrintSubHTML, $overallCountFilePathHsh_ref, $outDir);#->891
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 8_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1918
#	secondaryDependOnSub: currentTime|355
#
#<section ID="finishingTasks" num="8">
&printCMDLogOrFinishMessage("finishMessage");#->1918
close DEBUGLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	HTML [n=1]:
#		generateMasterHTML
#
#	general [n=7]:
#		checkSampleFiles, currentTime, printCMDLogOrFinishMessage
#		printCountTSVInput, readParameters, reportStatus
#		userDecisionToProceed
#
#	getTextInfo [n=1]:
#		getGeneNameAndLength
#
#	plotInR [n=2]:
#		pairwiseVennDiagram, tripleVennDiagram
#
#	specific [n=16]:
#		PCAWithFactoMiner, generateAllHeatmap, generateAllPCA
#		generateAllSampleDistHeatmap, generateComparisonCombination, generateDESeq1RScript
#		generateDESeq2RScript, generateEdgeRRScript, getAllCountInfo
#		getCountInput, launchAndMonitorAllCMD, overlapGeneExpDiffResults
#		printAllComparisonFoldChange, printAllSampleCountInfo, printCombineResults
#		runAllComparisons
#
#	taskManage [n=1]:
#		runAndCheckSerialTask
#
#====================================================================================================================================================#

sub PCAWithFactoMiner {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: generateAllPCA|394
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_PCAAndHeatmap|154
#	input: $allResultFilePathHsh_ref, $allSampleCountInfoHsh_ref, $allSampleCountPathHsh_ref, $factoMinerPath, $log2OrLinear, $valueType
#	output: none
#	toCall: &PCAWithFactoMiner($factoMinerPath, $allSampleCountInfoHsh_ref, $allSampleCountPathHsh_ref, $valueType, $log2OrLinear, $allResultFilePathHsh_ref);
#	calledInLine: 411
#....................................................................................................................................................#
	
	my ($factoMinerPath, $allSampleCountInfoHsh_ref, $allSampleCountPathHsh_ref, $valueType, $log2OrLinear, $allResultFilePathHsh_ref) = @_;

	my $outPrefix = join ".", ('allSamples', $valueType, $log2OrLinear);

	print "Running PCA using $log2OrLinear $valueType.\n";
	system "mkdir -p -m 777 $factoMinerPath/input/";
	system "mkdir -p -m 777 $factoMinerPath/output/";

	my @replicateStringAry;

	foreach my $geneID (sort {$a cmp $b} keys %{$allSampleCountInfoHsh_ref}) {
		foreach my $condition (sort {$a cmp $b} keys %{$allSampleCountInfoHsh_ref->{$geneID}}) {
			my $repNum = 0;
			foreach my $replicate (sort {$a <=> $b} keys %{$allSampleCountInfoHsh_ref->{$geneID}{$condition}}) {
				$repNum++;
			}
			push @replicateStringAry, "rep(\"$condition\", $repNum)";
		}
		last;
	}
	
	my $replicateString = join ",", @replicateStringAry;
	
	my $PC1vsPC2Path = "$factoMinerPath/output/$outPrefix.pca.PC1vsPC2.pdf";
	my $PC2vsPC3Path = "$factoMinerPath/output/$outPrefix.pca.PC2vsPC3.pdf";
	my $PC3vsPC4Path = "$factoMinerPath/output/$outPrefix.pca.PC3vsPC4.pdf";
	my $PC4vsPC5Path = "$factoMinerPath/output/$outPrefix.pca.PC4vsPC5.pdf";

	$allResultFilePathHsh_ref->{'PCA'}{$valueType}{$log2OrLinear}{'PC1vsPC2Path'} = $PC1vsPC2Path;
	$allResultFilePathHsh_ref->{'PCA'}{$valueType}{$log2OrLinear}{'PC2vsPC3Path'} = $PC2vsPC3Path;
	$allResultFilePathHsh_ref->{'PCA'}{$valueType}{$log2OrLinear}{'PC3vsPC4Path'} = $PC3vsPC4Path;
	$allResultFilePathHsh_ref->{'PCA'}{$valueType}{$log2OrLinear}{'PC4vsPC5Path'} = $PC4vsPC5Path;
	
	my $factoMinerScriptPath = "$factoMinerPath/input/$outPrefix.PCA.R";
	open (FACTOMINER, ">$factoMinerScriptPath");
	my $tsvDataPath = $allSampleCountPathHsh_ref->{$valueType}{$log2OrLinear};
	print FACTOMINER "library(FactoMineR);\n";
	print FACTOMINER "expData = read.table(\"$tsvDataPath\",header=T,sep=\"\\t\",dec=\".\",row.names=1);\n";
	print FACTOMINER "expData = as.data.frame(t(expData));\n";
	print FACTOMINER "expData = cbind.data.frame(as.factor(c($replicateString)), expData);\n";
	print FACTOMINER "colnames(expData)[1] = \"e\";\n";
	print FACTOMINER "res.pca = PCA(expData,quali.sup=1, graph=F);\n";
	print FACTOMINER "pdf(\"$PC1vsPC2Path\");\n";
	print FACTOMINER "plot(res.pca,habillage=1, axes=1:2);\n";
	print FACTOMINER "dev.off();\n";
	print FACTOMINER "pdf(\"$PC2vsPC3Path\");\n";
	print FACTOMINER "plot(res.pca,habillage=1, axes=2:3);\n";
	print FACTOMINER "dev.off();\n";
	print FACTOMINER "pdf(\"$PC3vsPC4Path\");\n";
	print FACTOMINER "plot(res.pca,habillage=1, axes=3:4);\n";
	print FACTOMINER "dev.off();\n";
	print FACTOMINER "pdf(\"$PC4vsPC5Path\");\n";
	print FACTOMINER "plot(res.pca,habillage=1, axes=4:5);\n";
	print FACTOMINER "dev.off();\n";
	print FACTOMINER "sink(file = \"$factoMinerPath/output/$outPrefix.pca.info.txt\");\n";
	print FACTOMINER "head(res.pca,habillage=1);\n";
	print FACTOMINER "sink();\n";
	close FACTOMINER;
	my $RCmd = "R --slave --no-restore --vanilla --file=$factoMinerScriptPath";
	system("$RCmd >$factoMinerScriptPath.error.log.txt 2>&1;");
}
sub checkSampleFiles {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: userDecisionToProceed|2467
#	appearInSub: >none
#	primaryAppearInSection: 2_processInputData|108
#	secondaryAppearInSection: >none
#	input: $refNameAry_ref, $sampleInfoPath
#	output: \%refNameHsh, \%sampleInfoHsh
#	toCall: my (\%refNameHsh, \%sampleInfoHsh) = &checkSampleFiles($sampleInfoPath, $refNameAry_ref);
#	calledInLine: 114, 306
#....................................................................................................................................................#

	#--- ($refNameHsh_ref, $sampleInfoHsh_ref) = &checkSampleFiles($sampleInfoPath, $refNameHsh_ref);#->294
	
	my ($sampleInfoPath, $refNameAry_ref) = @_;
	
	my %sampleInfoHsh = ();
	my %repNumCountHsh = ();
	my %refNameHsh = ();
	
	#---check bam file and collect reference sample name
	open (SAMINFO, "$sampleInfoPath");
	while (my $theLine = <SAMINFO>) {
		chomp $theLine;
		next if length $theLine < 2; #----empty line
		next if $theLine =~ m/^#/;
		my ($sampleName, $countFilePath)= split /\s+/, $theLine;
		$repNumCountHsh{$sampleName}++;
		if ((grep {$_ eq 'all'} @{$refNameAry_ref}) or (grep {$_ eq $sampleName} @{$refNameAry_ref})) {
			$refNameHsh{$sampleName}++;
		}
		
		#---replicate number is arranged by order of appearance, first is zero
		my $repNum = $repNumCountHsh{$sampleName} - 1;
		${${$sampleInfoHsh{$sampleName}}{$repNum}}{'countFilePath'} = $countFilePath;
		open (TEST, $countFilePath) || die "Can't open $countFilePath\n"; close TEST;
		print $sampleName." ".$repNum." expressCount checked\n";
	}
	close SAMINFO;
	
	my $refSampleToProcess = keys %refNameHsh;
	print "\n$refSampleToProcess reference sample will be ran,\n";
	foreach my $refName (sort {$a cmp $b} keys %refNameHsh) {
		print "\n============================================================\n";
		print "For sample $refName, the following comparison will be run,\n";
		foreach my $sampleName (sort {$a cmp $b} keys %sampleInfoHsh) {
			next if $sampleName eq $refName;
			my $refSample = $refName;
			my $qrySample = $sampleName;
			my $comparisonName = "$refSample.VS.$qrySample";
			my $refNum = keys %{$sampleInfoHsh{$refSample}};
			my $qryNum = keys %{$sampleInfoHsh{$qrySample}};
			print "Paired $comparisonName in a $refNum vs $qryNum replicate setting.\n";
		}
	}
	
	&userDecisionToProceed();#->2467
	
	return \%refNameHsh, \%sampleInfoHsh;
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|1918, reportStatus|2223
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|81, 8_finishingTasks|175
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 70, 1938, 1941, 1946, 2239
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub generateAllHeatmap {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: generateAllSampleDistHeatmap|416
#	appearInSub: >none
#	primaryAppearInSection: 6_PCAAndHeatmap|154
#	secondaryAppearInSection: >none
#	input: $allResultFilePathHsh_ref, $allSampleCountPathHsh_ref, $outDir, $sampleColumnAry_ref
#	output: none
#	toCall: &generateAllHeatmap($outDir, $allSampleCountPathHsh_ref, $sampleColumnAry_ref, $allResultFilePathHsh_ref);
#	calledInLine: 159
#....................................................................................................................................................#
	my ($outDir, $allSampleCountPathHsh_ref, $sampleColumnAry_ref, $allResultFilePathHsh_ref, ) = @_;
	
	#----generate all sample heatmaps
	my $distHeatMapPath = "$outDir/allSample/distHeatmap/";
	#foreach (qw/inputCount fpm fpkm/) {
	foreach (qw/inputCount/) {
		&generateAllSampleDistHeatmap($allSampleCountPathHsh_ref, $distHeatMapPath, $sampleColumnAry_ref, $_, $allResultFilePathHsh_ref, $outDir);#->416
	}
}
sub generateAllPCA {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: PCAWithFactoMiner|220
#	appearInSub: >none
#	primaryAppearInSection: 6_PCAAndHeatmap|154
#	secondaryAppearInSection: >none
#	input: $allResultFilePathHsh_ref, $allSampleCountInfoHsh_ref, $allSampleCountPathHsh_ref, $outDir
#	output: none
#	toCall: &generateAllPCA($outDir, $allSampleCountInfoHsh_ref, $allSampleCountPathHsh_ref, $allResultFilePathHsh_ref);
#	calledInLine: 160
#....................................................................................................................................................#

	my ($outDir, $allSampleCountInfoHsh_ref, $allSampleCountPathHsh_ref, $allResultFilePathHsh_ref) = @_;
	
	#----generate PCA with factominer
	#foreach (qw/vsd_inputCount vsd_fpm vsd_fpkm/) {
	foreach (qw/vsd_inputCount/) {
		&PCAWithFactoMiner("$outDir/allSample/PCA/", $allSampleCountInfoHsh_ref, $allSampleCountPathHsh_ref, $_, 'log2', $allResultFilePathHsh_ref);#->220
	}

}
sub generateAllSampleDistHeatmap {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: generateAllHeatmap|373
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_PCAAndHeatmap|154
#	input: $allResultFilePathHsh_ref, $allSampleCountPathHsh_ref, $countType, $distHeatMapPath, $outDir, $sampleColumnAry_ref
#	output: none
#	toCall: &generateAllSampleDistHeatmap($allSampleCountPathHsh_ref, $distHeatMapPath, $sampleColumnAry_ref, $countType, $allResultFilePathHsh_ref, $outDir);
#	calledInLine: 390
#....................................................................................................................................................#
	
	my ($allSampleCountPathHsh_ref, $distHeatMapPath, $sampleColumnAry_ref, $countType, $allResultFilePathHsh_ref, $outDir) = @_;
	

	system "mkdir -p -m 777 $distHeatMapPath/input";
	system "mkdir -p -m 777 $distHeatMapPath/output";
	
	print "Clustering all samples based on variance stablized $countType data from DESeq.\n";
	
	my $sampleColumnStr = join "\",\"", @{$sampleColumnAry_ref};
	$sampleColumnStr = "\"".$sampleColumnStr."\"";

	my $heatmapPath = "$distHeatMapPath/output/allHeatMap.$countType.heatmap.pdf";
	$allResultFilePathHsh_ref->{'distanceHeatmap'}{$countType} = $heatmapPath;

	my $allSampleDistHeatmapRScriptPath = "$outDir/allSample/distHeatmap/input/$countType.distHeatmapWithVsd.R";
	my $varianceStablizedCountDataPath = "$outDir/allSample/distHeatmap/output/$countType.varianceStablizedCountData.log.xls";
	open (R, ">$allSampleDistHeatmapRScriptPath");

	my $tsvDataPath = $allSampleCountPathHsh_ref->{$countType}{'linear'};
	
	print R "countTable <- read.table(\"$tsvDataPath\", header=TRUE, row.names=1);\n";
	print R "roundCountTable <- round(countTable, digits = 0);\n"; #----round to integer
	print R "conds <- factor( c($sampleColumnStr));\n";
	print R "library(\"DESeq\");\n";
	print R "cds <- newCountDataSet( roundCountTable, conds );\n";
	print R "cds <- estimateSizeFactors( cds );\n";
	print R "cdsBlind <- estimateDispersions( cds, method=\"blind\");\n";
	print R "vsd <- getVarianceStabilizedData( cdsBlind )\n";
	print R "write.table(vsd, file = \"$varianceStablizedCountDataPath\", sep = \"\\t\", col.names=NA, quote=FALSE);\n";
	print R "pdf(\"$heatmapPath\")\n";
	print R "dists <- dist( t( vsd ) );\n";
	print R "heatmap( as.matrix( dists ), symm=TRUE, scale=\"none\", margins=c(10,10), col = colorRampPalette(c(\"darkblue\",\"white\"))(100), labRow = paste( pData(cdsBlind)\$condition));\n";
	print R "invisible(dev.off());\n";
	close R;

	my $RCmd = "R --slave --no-restore --vanilla  --file=$allSampleDistHeatmapRScriptPath >$allSampleDistHeatmapRScriptPath.$countType.error.log.txt 2>&1";
	system("$RCmd;");
	
	$allSampleCountPathHsh_ref->{'vsd_'.$countType}{'log2'} = $varianceStablizedCountDataPath;
}
sub generateComparisonCombination {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: runAllComparisons|2244
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_runDEAnalyses|131
#	input: $allComparisonNameHsh_ref, $outDir, $refName, $sampleInfoHsh_ref
#	output: $comparisonInfoHsh_ref
#	toCall: my ($comparisonInfoHsh_ref) = &generateComparisonCombination($sampleInfoHsh_ref, $refName, $allComparisonNameHsh_ref, $outDir);
#	calledInLine: 2274
#....................................................................................................................................................#

	my ($sampleInfoHsh_ref, $refName, $allComparisonNameHsh_ref, $outDir) = @_;
	
	my $comparisonInfoHsh_ref;
	
	foreach my $sampleName (sort {$a cmp $b} keys %{$sampleInfoHsh_ref}) {
		next if $sampleName eq $refName;
		my $refSample = $refName;
		my $qrySample = $sampleName;
		my $comparisonName = "$refSample.VS.$qrySample";
		$comparisonInfoHsh_ref->{$comparisonName}{'refSample'} = $refSample;
		$comparisonInfoHsh_ref->{$comparisonName}{'qrySample'} = $qrySample;
		$allComparisonNameHsh_ref->{$comparisonName}++;
		
		foreach my $repNum (sort {$a <=> $b} keys %{$sampleInfoHsh_ref->{$refSample}}) {
			push @{$comparisonInfoHsh_ref->{$comparisonName}{'refExpressCount'}}, $sampleInfoHsh_ref->{$refSample}{$repNum}{'countFilePath'};
		}

		foreach my $repNum (sort {$a <=> $b} keys %{$sampleInfoHsh_ref->{$qrySample}}) {
			push @{$comparisonInfoHsh_ref->{$comparisonName}{'qryExpressCount'}}, $sampleInfoHsh_ref->{$qrySample}{$repNum}{'countFilePath'};
		}
		
		my $refRepNum = @{$comparisonInfoHsh_ref->{$comparisonName}{'refExpressCount'}};
		my $qryRepNum = @{$comparisonInfoHsh_ref->{$comparisonName}{'qryExpressCount'}};
		print "Paired $comparisonName in a $refRepNum vs $qryRepNum replicate setting.\n";

		my $rootOutDir = "$outDir/$refSample/$comparisonName/";
		$comparisonInfoHsh_ref->{$comparisonName}{'rootOutDir'} = $rootOutDir;
		system "mkdir -p -m 777 $rootOutDir"
		
	}

	return $comparisonInfoHsh_ref;
	
}
sub generateDESeq1RScript {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: runAllComparisons|2244
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_runDEAnalyses|131
#	input: $DESeq1FDR, $DESeq1FitType, $DESeq1Method, $DESeq1SharingMode, $comparisonInfoHsh_ref, $outDir
#	output: none
#	toCall: &generateDESeq1RScript($comparisonInfoHsh_ref, $DESeq1FDR, $DESeq1FitType, $DESeq1SharingMode, $DESeq1Method, $outDir);
#	calledInLine: 2280
#....................................................................................................................................................#

	my ($comparisonInfoHsh_ref, $DESeq1FDR, $DESeq1FitType, $DESeq1SharingMode, $DESeq1Method, $outDir) = @_;
	
	my $DESeqResultPathHsh_ref;
	
	foreach my $comparisonName (sort {$a cmp $b} keys %{$comparisonInfoHsh_ref}) {

		my $refSample = $comparisonInfoHsh_ref->{$comparisonName}{'refSample'};
		my $qrySample = $comparisonInfoHsh_ref->{$comparisonName}{'qrySample'};
		my $refSampleNum = @{$comparisonInfoHsh_ref->{$comparisonName}{'refExpressCount'}};
		my $qrySampleNum = @{$comparisonInfoHsh_ref->{$comparisonName}{'qryExpressCount'}};

		if ($qrySampleNum == 1 or $refSampleNum == 1) {
			$DESeq1FitType = 'local';
			$DESeq1SharingMode = 'fit-only';
			$DESeq1Method = 'blind';
		}

		my $DESeq1InDir = "$outDir/$refSample/$comparisonName/DESeq1/input";
		my $DESeq1OutDir = "$outDir/$refSample/$comparisonName/DESeq1/output";

		system "mkdir -p -m 777 $DESeq1OutDir";
		system "mkdir -p -m 777 $DESeq1InDir";

		my $DESeq1InputTSVPath = "$DESeq1InDir/readCount.xls";
		my $DESeq1SummaryPath = "$DESeq1OutDir/summary.log.txt";
		my $DESeq1ScrLog = "$outDir/$refSample/$comparisonName/DESeq1.ScreenLog.txt";
		my $DESeq1RScriptPath = "$DESeq1InDir/$comparisonName.DESeq1.R";
		my $DESeq1CMD = "R --slave --no-restore --vanilla --file=$DESeq1RScriptPath";
		$DESeq1CMD =~ s/\s+/ /g;
		my $DESeq1DEGResultPath = "$DESeq1OutDir/DESeq1.DEG.xls";
		my $DESeq1ModLog2FCCSVPath = "$DESeq1OutDir/moderated.log.foldChange.Results.csv";

		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq1InDir'} = $DESeq1InDir;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq1OutDir'} = $DESeq1OutDir;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq1InputTSVPath'} = $DESeq1InputTSVPath;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq1SummaryPath'} = $DESeq1SummaryPath;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq1ScrLog'} = $DESeq1ScrLog;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq1RScriptPath'} = $DESeq1RScriptPath;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq1CMD'} = $DESeq1CMD;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq1DEGResultPath'} = $DESeq1DEGResultPath;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq1ModLog2FCCSVPath'} = $DESeq1ModLog2FCCSVPath;

		my @sampleColumnAry = ();
		foreach (@{$comparisonInfoHsh_ref->{$comparisonName}{'refExpressCount'}}) {
			push @sampleColumnAry, "\"$refSample\"";
		}
		foreach (@{$comparisonInfoHsh_ref->{$comparisonName}{'qryExpressCount'}}) {
			push @sampleColumnAry, "\"$qrySample\"";
		}
		
		my $sampleColumnStr = join ",", @sampleColumnAry;
		
		open R, ">$DESeq1RScriptPath";

		print R "cat(\"Starting DESeq analysis.\\n\")\n";
		print R "cat(\"Reading count table.\\n\")\n";
		print R "countTable <- read.table(\"$DESeq1InputTSVPath\", header=TRUE, row.names=1)\n";
		print R "conds <- factor( c($sampleColumnStr))\n";
		print R "library(\"DESeq\")\n";
		print R "cat(\"Normalizing count data and estimating variance in mode.\\n\")\n";
		print R "cds <- newCountDataSet( countTable, conds )\n";
		print R "cds <- estimateSizeFactors( cds )\n";
		print R "cds <- estimateDispersions( cds, fitType=\"$DESeq1FitType\", sharingMode=\"$DESeq1SharingMode\", method=\"$DESeq1Method\" )\n";
		print R "cat(\"Plotting variance data.\\n\")\n";
		
		if ($DESeq1Method eq 'per-condition') {
			print R "plotDispEsts <- function( cds )";
			print R "{ plot(rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(cds, \"$refSample\")\$perGeneDispEsts, pch = '.', log=\"xy\")\n";
			print R "xg <- 10^seq( -.5, 5, length.out=300 )\n";
			print R "lines( xg, fitInfo(cds, \"$refSample\")\$dispFun( xg ), col=\"red\") }\n";
			print R "pdf(\"$DESeq1OutDir/DispEsts.$refSample.pdf\")\n";
			print R "plotDispEsts( cds )\n";
			print R "invisible(dev.off())\n";
		
			print R "plotDispEsts <- function( cds )";
			print R "{ plot(rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(cds, \"$qrySample\")\$perGeneDispEsts, pch = '.', log=\"xy\")\n";
			print R "xg <- 10^seq( -.5, 5, length.out=300 )\n";
			print R "lines( xg, fitInfo(cds, \"$qrySample\")\$dispFun( xg ), col=\"red\") }\n";
			print R "pdf(\"$DESeq1OutDir/DispEsts.$qrySample.pdf\")\n";
			print R "plotDispEsts( cds )\n";
			print R "invisible(dev.off())\n";
			
		} else {
			print R "plotDispEsts <- function( cds )";
			print R "{ plot(rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(cds)\$perGeneDispEsts, pch = '.', log=\"xy\")\n";
			print R "xg <- 10^seq( -.5, 5, length.out=300 )\n";
			print R "lines( xg, fitInfo(cds)\$dispFun( xg ), col=\"red\") }\n";
			print R "pdf(\"$DESeq1OutDir/dispEsts.pdf\")\n";
			print R "plotDispEsts( cds )\n";
			print R "invisible(dev.off())\n";
		}
		
		print R "cat(\"Performing nbinorm test.\\n\")\n";
		
		print R "res <- nbinomTest( cds, \"$refSample\", \"$qrySample\" )\n";
		print R "resSig <- res[ res\$padj < $DESeq1FDR, ]\n";
		
		print R "pdf(\"$DESeq1OutDir/pValueHistogram.pdf\")\n";
		print R "hist(res\$pval, breaks=100, col=\"skyblue\", border=\"slateblue\", main=\"\")\n";
		print R "invisible(dev.off())\n";
		
		print R "plotDE <- function(res) {plot(res\$baseMean, res\$log2FoldChange, log=\"x\", pch=20, cex=.3, col = ifelse( res\$padj < $DESeq1FDR, \"red\", \"black\" ))}\n";
		print R "pdf(\"$DESeq1OutDir/DE.pdf\")\n";
		print R "plotDE(res)\n";
		print R "invisible(dev.off())\n";

		print R "write.table( res, col.names=NA, sep =\"\t\", file=\"$DESeq1DEGResultPath\")\n";
		
		if ($qrySampleNum > 1 and $refSampleNum > 1) {#---has replicate
			print R "pdf(\"$DESeq1OutDir/$refSample.FCVsEXP.pdf\")\n";
			print R "ncu <- counts( cds, normalized=TRUE )[ , conditions(cds)==\"$refSample\"]\n";
			print R "plot( rowMeans(ncu), log2( ncu[,2] / ncu[,1] ), pch=\".\", log=\"x\")\n";
			print R "invisible(dev.off())\n";
			print R "pdf(\"$DESeq1OutDir/$qrySample.FCVsEXP.pdf\")\n";
			print R "ncu <- counts( cds, normalized=TRUE )[ , conditions(cds)==\"$qrySample\"]\n";
			print R "plot( rowMeans(ncu), log2( ncu[,2] / ncu[,1] ), pch=\".\", log=\"x\")\n";
			print R "invisible(dev.off())\n";

			print R "cat(\"Performing variance stabilizing transformation.\\n\")\n";
			print R "cdsBlind <- estimateDispersions( cds, method=\"blind\")\n";
			print R "vsd <- getVarianceStabilizedData( cdsBlind )\n";
			print R "library(\"vsn\")\n";
			print R "pdf(\"$DESeq1OutDir/meadSDAfterVarStabilize.pdf\")\n";
			print R "meanSdPlot(vsd)\n";
			print R "invisible(dev.off())\n";
			
			print R "cat(\"Moderating the fold change.\\n\")\n";
			print R "mod_lfc <- (rowMeans( vsd[, conditions(cds)==\"$qrySample\", drop=FALSE] ) - rowMeans( vsd[, conditions(cds)==\"$refSample\", drop=FALSE] ))\n";
			print R "write.csv( mod_lfc, file=\"$DESeq1ModLog2FCCSVPath\")\n";

			print R "lfc <- res\$log2FoldChange\n";
			print R "finite <- is.finite(lfc)\n";
			print R "largeNumber <- 7\n";
			print R "lfc <- ifelse(finite, lfc, sign(lfc) * largeNumber)\n";
			print R "logdecade <- 1 + round( log10( 1+rowMeans(counts(cdsBlind, normalized=TRUE))))\n";
			print R "colors <- colorRampPalette( c( \"gray\", \"blue\" ) )(6)[logdecade]\n";
			print R "pdf(\"$DESeq1OutDir/moderatedFoldChange.pdf\")\n";
			print R "plot( lfc, mod_lfc, pch=20, cex=.4, asp=1, col = ifelse( finite, colors, \"purple\"))\n";
			print R "abline( a=0, b=1, col=\"#40C04040\")\n";
			
			print R "cat(\"Plotting the expression distance heatmap of the samples.\\n\")\n";
			print R "dists <- dist( t( vsd ) )\n";
			print R "pdf(\"$DESeq1OutDir/sampleDist.heatmap.pdf\")\n";
			print R "heatmap( as.matrix( dists ), symm=TRUE, scale=\"none\", margins=c(10,10), col = colorRampPalette(c(\"darkblue\",\"white\"))(100), labRow = paste( pData(cdsBlind)\$condition))\n";
			print R "invisible(dev.off())\n";
			
			print R "cat(\"Finsihed DESeq analysis.\\n\")\n";
		
		} else {#---no replicate
		
			print R "cat(\"Performing variance stabilizing transformation.\\n\")\n";
			print R "vsd <- getVarianceStabilizedData( cds )\n";
			print R "library(\"vsn\")\n";
			print R "pdf(\"$DESeq1OutDir/meadSDAfterVarStabilize.pdf\")\n";
			print R "meanSdPlot(vsd)\n";
			print R "invisible(dev.off())\n";

			print R "cat(\"Moderating the fold change.\\n\")\n";
			print R "mod_lfc <- (rowMeans( vsd[, conditions(cds)==\"$qrySample\", drop=FALSE] ) - rowMeans( vsd[, conditions(cds)==\"$refSample\", drop=FALSE] ))\n";
			print R "write.csv( mod_lfc, file=\"$DESeq1ModLog2FCCSVPath\")\n";
			
			print R "cat(\"Finsihed DESeq analysis.\\n\")\n";

		}
		
		close R;
	}

}
sub generateDESeq2RScript {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: runAllComparisons|2244
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_runDEAnalyses|131
#	input: $DESeq2FDR, $DESeq2FitType, $comparisonInfoHsh_ref, $outDir
#	output: none
#	toCall: &generateDESeq2RScript($comparisonInfoHsh_ref, $DESeq2FDR, $DESeq2FitType, $outDir);
#	calledInLine: 2284
#....................................................................................................................................................#

	my ($comparisonInfoHsh_ref, $DESeq2FDR, $DESeq2FitType, $outDir) = @_;
	
	foreach my $comparisonName (sort {$a cmp $b} keys %{$comparisonInfoHsh_ref}) {

		my $refSample = $comparisonInfoHsh_ref->{$comparisonName}{'refSample'};
		my $qrySample = $comparisonInfoHsh_ref->{$comparisonName}{'qrySample'};
		my $refSampleNum = @{$comparisonInfoHsh_ref->{$comparisonName}{'refExpressCount'}};
		my $qrySampleNum = @{$comparisonInfoHsh_ref->{$comparisonName}{'qryExpressCount'}};

		my $DESeq2InDir = "$outDir/$refSample/$comparisonName/DESeq2/input";
		my $DESeq2OutDir = "$outDir/$refSample/$comparisonName/DESeq2/output";

		system "mkdir -p -m 777 $DESeq2OutDir";
		system "mkdir -p -m 777 $DESeq2InDir";

		my $DESeq2InputTSVPath = "$DESeq2InDir/readCount.xls";
		my $DESeq2ScrLog = "$outDir/$refSample/$comparisonName/DESeq2.ScreenLog.txt";
		my $DESeq2RScriptPath = "$DESeq2InDir/$comparisonName.DESeq2.R";
		my $DESeq2CMD = "R --slave --no-restore --vanilla --file=$DESeq2RScriptPath";
		$DESeq2CMD =~ s/\s+/ /g;
		my $DESeq2DEGResultPath = "$DESeq2OutDir/DESeq2.DEG.xls";

		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq2InDir'} = $DESeq2InDir;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq2OutDir'} = $DESeq2OutDir;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq2InputTSVPath'} = $DESeq2InputTSVPath;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq2ScrLog'} = $DESeq2ScrLog;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq2RScriptPath'} = $DESeq2RScriptPath;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq2CMD'} = $DESeq2CMD;
		$comparisonInfoHsh_ref->{$comparisonName}{'DESeq2DEGResultPath'} = $DESeq2DEGResultPath;

		my @sampleColumnAry = ();
		foreach (@{$comparisonInfoHsh_ref->{$comparisonName}{'refExpressCount'}}) {
			push @sampleColumnAry, "\"ref_$refSample\"";
		}
		foreach (@{$comparisonInfoHsh_ref->{$comparisonName}{'qryExpressCount'}}) {
			push @sampleColumnAry, "\"qry_$qrySample\"";
		}
		
		my $sampleColumnStr = join ",", @sampleColumnAry;
		
		open R, ">$DESeq2RScriptPath";
		print R "library(\"DESeq2\")\n";
		print R "cat(\"Starting DESeq2 analysis.\\n\")\n";
		print R "cat(\"Reading count table.\\n\")\n";
		print R "countData <- read.table(\"$DESeq2InputTSVPath\", header=TRUE, row.names=1)\n";
		print R "colData <- DataFrame(condition=factor(c($sampleColumnStr)))\n";
		print R "cat(\"Normalizing count data and estimating variance.\\n\")\n";
		print R "dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)\n";
		print R "dds <- estimateSizeFactors(dds)\n";
		print R "dds <- estimateDispersions(dds, fitType=\"$DESeq2FitType\", maxit=100)\n";
		print R "cat(\"Performing nbinorm test.\\n\")\n";
		print R "dds <- nbinomWaldTest(dds)\n";
		print R "results <- results(dds)\n";
		print R "pdf(\"$DESeq2OutDir/plotMA.pdf\")\n";
		print R "plotMA(dds, pvalCutoff=$DESeq2FDR)\n";
		print R "invisible(dev.off())\n";
		print R "pdf(\"$DESeq2OutDir/pValueHistogram.pdf\")\n";
		print R "hist(results\$pval, breaks=100, col=\"skyblue\", border=\"slateblue\", main=\"\")\n";
		print R "invisible(dev.off())\n";
		print R "write.table(as.data.frame(results), col.names=NA, sep =\"\t\", file=\"$DESeq2DEGResultPath\")\n";
		close R;
	}

}
sub generateEdgeRRScript {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: runAllComparisons|2244
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_runDEAnalyses|131
#	input: $BCV, $comparisonInfoHsh_ref, $outDir
#	output: none
#	toCall: &generateEdgeRRScript($comparisonInfoHsh_ref, $BCV, $outDir);
#	calledInLine: 2287
#....................................................................................................................................................#
	
	my ($comparisonInfoHsh_ref, $BCV, $outDir) = @_;

	foreach my $comparisonName (sort {$a cmp $b} keys %{$comparisonInfoHsh_ref}) {
		my $refSample = $comparisonInfoHsh_ref->{$comparisonName}{'refSample'};
		my $qrySample = $comparisonInfoHsh_ref->{$comparisonName}{'qrySample'};
		
		my $refSampleNum = @{$comparisonInfoHsh_ref->{$comparisonName}{'refExpressCount'}};
		my $qrySampleNum = @{$comparisonInfoHsh_ref->{$comparisonName}{'qryExpressCount'}};

		my $noReplicate = 'no';
		if ($qrySampleNum == 1 or $refSampleNum == 1) {
			$noReplicate = 'yes';
		}
		
		my $EdgeRInDir = "$outDir/$refSample/$comparisonName/EdgeR/input";
		my $EdgeROutDir = "$outDir/$refSample/$comparisonName/EdgeR/output";
		system "mkdir -p -m 777 $EdgeRInDir";
		system "mkdir -p -m 777 $EdgeROutDir";
		my $EdgeRInputTSVPath = "$EdgeRInDir/readCountFromHTSeqCount.xls";
		my $EdgeRSummaryPath = "$EdgeROutDir/summary.log.txt";
		$comparisonInfoHsh_ref->{$comparisonName}{'EdgeRInputTSVPath'} = $EdgeRInputTSVPath;
		$comparisonInfoHsh_ref->{$comparisonName}{'EdgeROutDir'} = $EdgeROutDir;
		$comparisonInfoHsh_ref->{$comparisonName}{'EdgeRScrLog'} = "$outDir/$refSample/$comparisonName/EdgeR.ScreenLog.txt";
		open (SUMMARY, ">$EdgeRSummaryPath");
		print SUMMARY "FDR10\tFDR5\tFDR1\n";
		close SUMMARY;
		
		my @sampleColumnAry = ();
		foreach (@{$comparisonInfoHsh_ref->{$comparisonName}{'refExpressCount'}}) {
			push @sampleColumnAry, "\"$refSample\"";
		}
		foreach (@{$comparisonInfoHsh_ref->{$comparisonName}{'qryExpressCount'}}) {
			push @sampleColumnAry, "\"$qrySample\"";
		}
	
		my $sampleColumnStr = join ",", @sampleColumnAry;
		
		my $EdgeRRScriptPath = "$EdgeRInDir/$comparisonName.EdgeR.R";
		$comparisonInfoHsh_ref->{$comparisonName}{'EdgeRRScriptPath'} = $EdgeRRScriptPath;
	
		my $EdgeRCMD = "R --slave --no-restore --vanilla --file=$EdgeRRScriptPath";
		$EdgeRCMD =~ s/\s+/ /g;
		$comparisonInfoHsh_ref->{$comparisonName}{'EdgeRCMD'} = $EdgeRCMD;
	
		open R, ">$EdgeRRScriptPath";
		my $EdgeRGenesExactTestTSVPath = "$EdgeROutDir/EdgeR.exactTest.DEG.xls";
		my $EdgeRExpressionValuesTSVPath = "$EdgeROutDir/EdgeR.allExpressionValues.xls";
		my $volcanoPDFPath = "$EdgeROutDir/plotSmear.volcano.pdf";
		my $DEGNumPath = "$EdgeROutDir/DEGNum.txt";
		my $plotBCVPDFPath = "$EdgeROutDir/plotBCV.dispersions.pdf";
		my $PCAPDFPath = "$EdgeROutDir/plotMDS.PCA.pdf";
		my $sizeFactorPath = "$EdgeROutDir/normalizeFactor.txt";
		
		$comparisonInfoHsh_ref->{$comparisonName}{'EdgeRGenesExactTestTSVPath'} = $EdgeRGenesExactTestTSVPath;
		$comparisonInfoHsh_ref->{$comparisonName}{'EdgeRExpressionValuesTSVPath'} = $EdgeRExpressionValuesTSVPath;
	
		system ("mkdir -p -m 777 $EdgeROutDir");#---create the DIR first
		print R "library(edgeR)\n";
		print R "cat(\"Reading count table.\\n\")\n";
		print R "counttable <- read.table(\"$EdgeRInputTSVPath\", header=TRUE, row.names=1)\n";
		print R "sampleGroup <- factor(c($sampleColumnStr))\n";
		print R "DGEListObj <- DGEList(counts=counttable, group=sampleGroup)\n";
		print R "cat(\"Normalizing count table.\\n\")\n";
		print R "DGEListObj\$samples\$lib.size <- colSums(DGEListObj\$counts)\n";
		print R "DGEListObj <- calcNormFactors(DGEListObj)\n";
		print R "normalizeFactor<-DGEListObj\$samples\n";
		print R "write.table(normalizeFactor, col.names=NA, sep =\"\\t\", file=\"$sizeFactorPath\")\n";

		if ($noReplicate eq 'no') {
			print R "pdf(\"$PCAPDFPath\")\n";
			print R "plotMDS(DGEListObj)\n";
			print R "invisible(dev.off())\n";
			print R "cat(\"Estimating dispersion.\\n\")\n";
			print R "DGEListObj <- estimateCommonDisp(DGEListObj, verbose=TRUE)\n";
			print R "DGEListObj <- estimateTagwiseDisp(DGEListObj, verbose=TRUE)\n";
			print R "pdf(\"$plotBCVPDFPath\")\n";
			print R "plotBCV(DGEListObj)\n";
			print R "invisible(dev.off())\n";
			print R "cat(\"performing exact test.\\n\")\n";
			print R "exactTestObj <- exactTest(DGEListObj, pair=c(\"$refSample\",\"$qrySample\"))\n";
		} else {
			print R "cat(\"performing exact test.\\n\")\n";
			print R "exactTestObj <- exactTest(DGEListObj, pair=c(\"$refSample\",\"$qrySample\") ,dispersion=$BCV^2)\n";
		}

		print R "allGenes <-topTags(exactTestObj, n=9999999999, adjust.method = \"BH\")\n";
		print R "write.table(allGenes, col.names=NA, sep =\"\\t\", file=\"$EdgeRGenesExactTestTSVPath\")\n";
		print R "allExpressionValues <- cpm(DGEListObj)\n";
		print R "cat(\"Outputing DEG data.\\n\")\n";
		print R "write.table(allExpressionValues, col.names=NA, sep =\"\t\", file=\"$EdgeRExpressionValuesTSVPath\")\n";
		print R "allDEG <- decideTestsDGE(exactTestObj, adjust.method = \"BH\", p.value = 0.05)\n";
		print R "DEGNum <- summary(allDEG)\n";
		print R "write.table(DEGNum, col.names=NA, sep =\"\\t\", file=\"$DEGNumPath\")\n";
		print R "detags <- rownames(DGEListObj)[as.logical(allDEG)]\n";
		print R "pdf(\"$volcanoPDFPath\")\n";
		print R "plotSmear(exactTestObj, de.tags=detags)\n";
		print R "abline(h=c(-1, 1), col=\"blue\")\n";
		print R "invisible(dev.off())\n";
		print R "cat(\"Finished EdgeR analysis.\\n\")\n";
		close R;
	}
}
sub generateMasterHTML {
#....................................................................................................................................................#
#	subroutineCategory: HTML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 7_outputHTML|165
#	secondaryAppearInSection: >none
#	input: $allResultFilePathHsh_ref, $dirToPrintSubHTML, $outDir, $overallCountFilePathHsh_ref
#	output: none
#	toCall: &generateMasterHTML($allResultFilePathHsh_ref, $dirToPrintSubHTML, $overallCountFilePathHsh_ref, $outDir);
#	calledInLine: 170
#....................................................................................................................................................#

	my ($allResultFilePathHsh_ref, $dirToPrintSubHTML, $overallCountFilePathHsh_ref, $outDir) = @_;

	#$allResultFilePathHsh_ref->{"overlapBetweenMethodsPDF"}{$upDnAll."RegulatedGenes"}{"cutoff.RPKM.$minRPKM.FC$minLinearFC"}{$comparisonName} = $mergErD1D2PdfPath;
	#$allResultFilePathHsh_ref->{'fullResultTable'}{$FDRTag}{$comparisonName} = $combinedResultInfoPath;
	#$allResultFilePathHsh_ref->{'PCA'}{$valueType}{$log2OrLinear}{'PC1vsPC2Path'} = $PC1vsPC2Path;
	#$allResultFilePathHsh_ref->{'distanceHeatmap'}{$countType} = $heatmapPath;
	
	system ("mkdir -pm 777 $dirToPrintSubHTML");
	my $overlapBetweenMethodsPDFHTML = "$dirToPrintSubHTML/overlapBetweenMethodsPDF.html";
	my $resultTableHTML = "$dirToPrintSubHTML/resultTable.html";
	my $PCAHTML = "$dirToPrintSubHTML/PCA.html";
	my $distanceHeatmapHTML = "$dirToPrintSubHTML/distanceHeatmap.html";
	my $overallCountHTML = "$dirToPrintSubHTML/overallCount.html";
	my $masterIndexHTML = "$outDir/master_index.html";

	#---get the relative paths
	my $rel_overlapBetweenMethodsPDFHTML = File::Spec->abs2rel($overlapBetweenMethodsPDFHTML, $outDir);
	my $rel_resultTableHTML = File::Spec->abs2rel($resultTableHTML, $outDir);
	my $rel_PCAHTML = File::Spec->abs2rel($PCAHTML, $outDir);
	my $rel_overallCountHTML = File::Spec->abs2rel($overallCountHTML, $outDir);
	my $rel_distanceHeatmapHTML = File::Spec->abs2rel($distanceHeatmapHTML, $outDir);
	
	open (HTML, ">", $overlapBetweenMethodsPDFHTML);
	print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print HTML join "", ('<html>', "\n");
	print HTML join "", ("<head><title>Overlapping of differential expression results between different methods</title></head>\n");
	print HTML join "", ('<body>', "\n");
	print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
	print HTML join "", ('</div><div>', "\n");
	
	foreach my $upDnAll (sort keys %{$allResultFilePathHsh_ref->{"overlapBetweenMethodsPDF"}}) {
		print HTML join "", ("<h2>$upDnAll</h2>\n");
		print HTML join "", ('<ul>', "\n");
		foreach my $cutoff (sort keys %{$allResultFilePathHsh_ref->{"overlapBetweenMethodsPDF"}{$upDnAll}}) {
			print HTML join "", ("<h4>$cutoff</h4>\n");
			print HTML join "", ('<ul>', "\n");
			foreach my $comparisonName (sort keys %{$allResultFilePathHsh_ref->{"overlapBetweenMethodsPDF"}{$upDnAll}{$cutoff}}) {
				my $pdfLink = $allResultFilePathHsh_ref->{"overlapBetweenMethodsPDF"}{$upDnAll}{$cutoff}{$comparisonName};
				my $rel_path = File::Spec->abs2rel($pdfLink, $dirToPrintSubHTML);
				print HTML "<li><a href=\'$rel_path\'>$comparisonName</a>\n";
			}
			print HTML join "", ('</ul>', "\n");
		}
		print HTML join "", ('</ul>', "\n");

	}
	close HTML;
	
	open (HTML, ">", $resultTableHTML);
	print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print HTML join "", ('<html>', "\n");
	print HTML join "", ("<head><title>Combined differential expression result in tables</title></head>\n");
	print HTML join "", ('<body>', "\n");
	print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
	print HTML join "", ('</div><div>', "\n");
	
	foreach my $fullOrAbstract (qw/full abstract/) {
		foreach my $regularOrSenseAntisenseSameRow (qw/regular senseAntisenseSameRow/) {
			foreach my $FDRTag (sort keys %{$allResultFilePathHsh_ref->{$fullOrAbstract.$regularOrSenseAntisenseSameRow.'ResultTable'}}) {
				print HTML join "", ("<h2>$fullOrAbstract $regularOrSenseAntisenseSameRow table for FDR at $FDRTag</h2>\n");
				print HTML join "", ('<ul>', "\n");
				foreach my $comparisonName (sort keys %{$allResultFilePathHsh_ref->{$fullOrAbstract.$regularOrSenseAntisenseSameRow.'ResultTable'}{$FDRTag}}) {
					my $xlsLink = $allResultFilePathHsh_ref->{$fullOrAbstract.$regularOrSenseAntisenseSameRow.'ResultTable'}{$FDRTag}{$comparisonName};
					my $rel_path = File::Spec->abs2rel($xlsLink, $dirToPrintSubHTML);
					print HTML "<li><a href=\'$rel_path\'>$comparisonName</a>\n";
				}
				print HTML join "", ('</ul>', "\n");
			}
		}
	}
	close HTML;

	open (HTML, ">", $distanceHeatmapHTML);
	print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print HTML join "", ('<html>', "\n");
	print HTML join "", ("<head><title>Distance Heatmap</title></head>\n");
	print HTML join "", ('<body>', "\n");
	print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
	print HTML join "", ('</div><div>', "\n");
	
	foreach my $countType (sort keys %{$allResultFilePathHsh_ref->{"distanceHeatmap"}}) {
		print HTML join "", ("<h4>heatmap based on $countType</h4>\n");
		print HTML join "", ('<ul>', "\n");
		my $pdfLink = $allResultFilePathHsh_ref->{"distanceHeatmap"}{$countType};
		my $rel_path = File::Spec->abs2rel($pdfLink, $dirToPrintSubHTML);
		print HTML "<li><a href=\'$rel_path\'>heatmap</a>\n";
		print HTML join "", ('</ul>', "\n");
	}
	close HTML;

	open (HTML, ">", $PCAHTML);
	print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print HTML join "", ('<html>', "\n");
	print HTML join "", ("<head><title>Principle Component Analyses</title></head>\n");
	print HTML join "", ('<body>', "\n");
	print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
	print HTML join "", ('</div><div>', "\n");
	
	foreach my $valueType (sort keys %{$allResultFilePathHsh_ref->{"PCA"}}) {
		print HTML join "", ("<h2>value $valueType</h2>\n");
		print HTML join "", ('<ul>', "\n");
		foreach my $log2OrLinear (sort keys %{$allResultFilePathHsh_ref->{"PCA"}{$valueType}}) {
			print HTML join "", ("<h4>scale $log2OrLinear</h4>\n");
			print HTML join "", ('<ul>', "\n");
			foreach my $PCCombination (sort keys %{$allResultFilePathHsh_ref->{"PCA"}{$valueType}{$log2OrLinear}}) {
				my $pdfLink = $allResultFilePathHsh_ref->{"PCA"}{$valueType}{$log2OrLinear}{$PCCombination};
				my $rel_path = File::Spec->abs2rel($pdfLink, $dirToPrintSubHTML);
				print HTML "<li><a href=\'$rel_path\'>$PCCombination</a>\n";
			}
			print HTML join "", ('</ul>', "\n");
		}
		print HTML join "", ('</ul>', "\n");
	}
	close HTML;

	open (HTML, ">", $overallCountHTML);
	print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print HTML join "", ('<html>', "\n");
	print HTML join "", ("<head><title>Overall count of differentially expressed genes in different reference samples</title></head>\n");
	print HTML join "", ('<body>', "\n");
	print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
	print HTML join "", ('</div><div>', "\n");
	
	foreach my $FDRTag (keys %{$overallCountFilePathHsh_ref}) {
		print HTML join "", ("<h2>FDR at $FDRTag</h2>\n");
		print HTML join "", ('<ul>', "\n");
		foreach my $refName (sort keys %{$overallCountFilePathHsh_ref->{$FDRTag}}) {
			my $xlsLink = $overallCountFilePathHsh_ref->{$FDRTag}{$refName};
			my $rel_path = File::Spec->abs2rel($xlsLink, $dirToPrintSubHTML);
			print HTML "<li><a href=\'$rel_path\'>$refName vs the others</a>\n";
		}
		print HTML join "", ('</ul>', "\n");
	}
	close HTML;
	
	open (HTML, ">", $masterIndexHTML);
	print HTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print HTML join "", ('<html>', "\n");
	print HTML join "", ("<head><title>Summary of differential expression analyses</title></head>\n");
	print HTML join "", ('<body>', "\n");
	print HTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
	print HTML join "", ('</div><div>', "\n");
	print HTML join "", ("<h4>Click the following links to access the relevant information</h4>\n");
	print HTML join "", ('<ul>', "\n");
	
	print HTML "<li><a href=\'$rel_overlapBetweenMethodsPDFHTML\'>Overlapping between method DeSeq1, DeSeq2 and EdgeR</a>\n";
	print HTML "<li><a href=\'$rel_resultTableHTML\'>Tables of the combined results for each comparison</a>\n";
	print HTML "<li><a href=\'$rel_PCAHTML\'>Graphs of Principle component analyses</a>\n";
	print HTML "<li><a href=\'$rel_distanceHeatmapHTML\'>Inter-sample distance heatmap</a>\n";
	print HTML "<li><a href=\'$rel_overallCountHTML\'>Overall count of differentially expressed genes in each reference sample</a>\n";
	close HTML;
}
sub getAllCountInfo {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: reportStatus|2223
#	appearInSub: >none
#	primaryAppearInSection: 2_processInputData|108
#	secondaryAppearInSection: >none
#	input: $countType, $geneNameLenHsh_ref, $minCountPerReplicate, $sampleInfoHsh_ref
#	output: $allSampleCountInfoHsh_ref
#	toCall: my ($allSampleCountInfoHsh_ref) = &getAllCountInfo($sampleInfoHsh_ref, $countType, $geneNameLenHsh_ref, $minCountPerReplicate);
#	calledInLine: 116
#....................................................................................................................................................#
	
	my ($sampleInfoHsh_ref, $countType, $geneNameLenHsh_ref, $minCountPerReplicate) = @_;

	my $allSampleCountInfoHsh_ref;
	
	foreach my $sampleName (keys %{$sampleInfoHsh_ref}) {
		foreach my $repNum (keys %{$sampleInfoHsh_ref->{$sampleName}}) {
			
			&reportStatus("Reading counts from $sampleName replicate $repNum", 10, "\n");#->2223
			
			open (COUNTFILE, "<", $sampleInfoHsh_ref->{$sampleName}{$repNum}{'countFilePath'});
			if ($countType eq 'express_est_count' or $countType eq 'express_eff_count' or $countType eq 'express_uniq_count') {
				my $header = <COUNTFILE>;
				my @headerSplitAry = split /\t+/, $header;
				die "Bad express count file format. quitting" if @headerSplitAry != 14;
				while (<COUNTFILE>) {
					chomp $_;
					my %tmpHsh = ();
					my ($bundle_id, $geneID, $length, $eff_length, $tot_counts, $uniq_counts , $est_counts, $eff_counts, $ambig_distr_alpha, $ambig_distr_beta, $fpkm, $fpkm_conf_low, $fpkm_conf_high, $solvable) = split /\t+/, $_;
					if ($geneNameLenHsh_ref->{$geneID}) {
						$tmpHsh{'express_est_count'} = $est_counts;
						$tmpHsh{'express_eff_count'} = $eff_counts;
						$tmpHsh{'express_uniq_count'} = $uniq_counts;
						$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'inputCount'} = sprintf "%.0f", $tmpHsh{$countType};
					}
				}
			} elsif ($countType eq 'twoColumn') {
				while (<COUNTFILE>) {
					chomp $_;
					next if length $_ < 3;
					my ($geneID, $count) = split /\t+/, $_;
					if ($geneNameLenHsh_ref->{$geneID}) {
						$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'inputCount'} = sprintf "%.0f", $count;
					}
				}
			} else {
				die "invalid countType = $countType\n";
			}
			close COUNTFILE;
			
			#---add zero to allSampleCountInfoHsh_ref if geneID is absent
			foreach my $geneID (keys %{$geneNameLenHsh_ref}) {
				$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'inputCount'} = 0 if not $allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'inputCount'};
			}
			
		}
	}
	
	#-----remove genes less than certain count
	foreach my $geneID (keys %{$allSampleCountInfoHsh_ref}) {
		my $countInGene = 0;
		my $totalReNum = 0;
		foreach my $sampleName (keys %{$allSampleCountInfoHsh_ref->{$geneID}}) {
			foreach my $repNum (keys %{$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}}) {
				$countInGene += $allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'inputCount'};
				$totalReNum++;
			}
		}
		if ($countInGene/$totalReNum < $minCountPerReplicate) {
			delete $allSampleCountInfoHsh_ref->{$geneID};
			delete $geneNameLenHsh_ref->{$geneID};
		}
	}

	#---get the total number of counts
	my $sampleTotalFragNumHsh_ref = {};
	foreach my $geneID (keys %{$allSampleCountInfoHsh_ref}) {
		foreach my $sampleName (keys %{$allSampleCountInfoHsh_ref->{$geneID}}) {
			foreach my $repNum (keys %{$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}}) {
				$sampleTotalFragNumHsh_ref->{$sampleName}{$repNum} = 0 if not $sampleTotalFragNumHsh_ref->{$sampleName}{$repNum};
				$sampleTotalFragNumHsh_ref->{$sampleName}{$repNum} += $allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'inputCount'};
			}
		}
	}
	
	#---calculate the fragment per million
	foreach my $geneID (keys %{$allSampleCountInfoHsh_ref}) {
		foreach my $sampleName (keys %{$allSampleCountInfoHsh_ref->{$geneID}}) {
			foreach my $repNum (keys %{$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}}) {
				$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'fpm'} = $allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'inputCount'}/($sampleTotalFragNumHsh_ref->{$sampleName}{$repNum}/1000000);
				$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'fpkm'} = ($allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'inputCount'}/($sampleTotalFragNumHsh_ref->{$sampleName}{$repNum}/1000000))/($geneNameLenHsh_ref->{$geneID}{'len'}/1000);
			}
		}
	}

	return ($allSampleCountInfoHsh_ref);

}
sub getCountInput {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: launchAndMonitorAllCMD|1219
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $allSampleCountInfoHsh_ref, $qrySample, $refSample
#	output: $countInfoHsh_ref
#	toCall: my ($countInfoHsh_ref) = &getCountInput($refSample, $qrySample, $allSampleCountInfoHsh_ref);
#	calledInLine: 1348
#....................................................................................................................................................#

	my ($refSample, $qrySample, $allSampleCountInfoHsh_ref) = @_;
	
	my $countInfoHsh_ref = {};
	
	foreach my $geneID (keys %{$allSampleCountInfoHsh_ref}) {
		foreach my $sampleName (($refSample, $qrySample)) {
			foreach my $repNum (keys %{$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}}) {
				$countInfoHsh_ref->{$geneID}{$sampleName}{$repNum} = $allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{'inputCount'};
			}
		}
	}
	return $countInfoHsh_ref;
	
}
sub getGeneNameAndLength {
#....................................................................................................................................................#
#	subroutineCategory: getTextInfo
#	dependOnSub: reportStatus|2223
#	appearInSub: >none
#	primaryAppearInSection: 2_processInputData|108
#	secondaryAppearInSection: >none
#	input: $geneNameLenPath
#	output: $geneNameLenHsh_ref
#	toCall: my ($geneNameLenHsh_ref) = &getGeneNameAndLength($geneNameLenPath);
#	calledInLine: 115
#....................................................................................................................................................#

	my ($geneNameLenPath) = @_;
	
	my $geneNameLenHsh_ref;
	
	open (GENENAME, "<", $geneNameLenPath);
	<GENENAME>;
	while (<GENENAME>) {
		chomp $_;
		next if length $_ < 3;
		my ($geneID, $geneName, $strnd, $cntg, $locationTag, $geneRng, $exonRng, $CDSRng, $UTR3Rng, $UTR5Rng) = split /\t/;
		$geneNameLenHsh_ref->{$geneID}{'name'} = $geneName;
		$geneNameLenHsh_ref->{$geneID}{'len'} = $exonRng;
		$geneNameLenHsh_ref->{$geneID}{'locationTag'} = $locationTag;
	}
	close GENENAME;
	
	my $geneNum = keys %{$geneNameLenHsh_ref};
	print "\n";
	&reportStatus("Totally $geneNum gene names read", 10, "\n");#->2223

	return $geneNameLenHsh_ref;

}
sub launchAndMonitorAllCMD {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: getCountInput|1156, printCountTSVInput|2122
#	appearInSub: runAllComparisons|2244
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_runDEAnalyses|131
#	input: $allSampleCountInfoHsh_ref, $comparisonInfoHsh_ref, $countType, $outDir
#	output: none
#	toCall: &launchAndMonitorAllCMD($comparisonInfoHsh_ref, $outDir, $allSampleCountInfoHsh_ref, $countType);
#	calledInLine: 2290
#....................................................................................................................................................#

	my ($comparisonInfoHsh_ref, $outDir, $allSampleCountInfoHsh_ref, $countType) = @_;
	
	my $processInfoHsh_ref;#--- queued, running or finished
	
	open (CMDLOG, ">$outDir/launchCmd.log.txt");
	
	#---define the processInfoHsh for checking and running
	foreach my $comparisonName (sort {$a cmp $b} keys %{$comparisonInfoHsh_ref}) {

		my $DESeq1CMD = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq1CMD'};
		my $DESeq2CMD = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq2CMD'};
		my $EdgeRCMD = $comparisonInfoHsh_ref->{$comparisonName}{'EdgeRCMD'};

		my $DESeq1InputTSVPath = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq1InputTSVPath'};
		my $DESeq2InputTSVPath = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq2InputTSVPath'};
		my $EdgeRInputTSVPath = $comparisonInfoHsh_ref->{$comparisonName}{'EdgeRInputTSVPath'};

		my $EdgeRGenesExactTestTSVPath = $comparisonInfoHsh_ref->{$comparisonName}{'EdgeRGenesExactTestTSVPath'};
		my $EdgeRExpressionValuesTSVPath = $comparisonInfoHsh_ref->{$comparisonName}{'EdgeRExpressionValuesTSVPath'};

		my $DESeq1DEGResultPath = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq1DEGResultPath'};
		my $DESeq1ModLog2FCCSVPath = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq1ModLog2FCCSVPath'};
		my $DESeq2DEGResultPath = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq2DEGResultPath'};

		my $DESeq1ScrLog = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq1ScrLog'};
		my $DESeq2ScrLog = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq2ScrLog'};
		my $EdgeRScrLog = $comparisonInfoHsh_ref->{$comparisonName}{'EdgeRScrLog'};
		

		@{$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'inFile'}} = ($DESeq1InputTSVPath);
		@{$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'outFile'}} = ($DESeq1DEGResultPath,$DESeq1ModLog2FCCSVPath);
		$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'inTSV'} = $DESeq1InputTSVPath;
		$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'status'} = "queued";
		$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'psGrepStr'} = "ps xu| grep -e \'$DESeq1CMD\' | grep -v grep";;
		$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'mCheckStr'} = "R --slave --no-restore --vanilla";
		$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'allinFileExist'} = "no";
		$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'alloutFileExist'} = "no";
		$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'runCMD'} = $DESeq1CMD;
		$processInfoHsh_ref->{$comparisonName}{"DESeq1CMD"}{'scrLog'} = $DESeq1ScrLog;

		@{$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'inFile'}} = ($DESeq2InputTSVPath);
		@{$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'outFile'}} = ($DESeq2DEGResultPath);
		$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'inTSV'} = $DESeq2InputTSVPath;
		$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'status'} = "queued";
		$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'psGrepStr'} = "ps xu| grep -e \'$DESeq2CMD\' | grep -v grep";;
		$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'mCheckStr'} = "R --slave --no-restore --vanilla";
		$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'allinFileExist'} = "no";
		$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'alloutFileExist'} = "no";
		$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'runCMD'} = $DESeq2CMD;
		$processInfoHsh_ref->{$comparisonName}{"DESeq2CMD"}{'scrLog'} = $DESeq2ScrLog;

		@{$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'inFile'}} = ($EdgeRInputTSVPath);
		@{$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'outFile'}} = ($EdgeRGenesExactTestTSVPath, $EdgeRExpressionValuesTSVPath);
		$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'inTSV'} = $EdgeRInputTSVPath;
		$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'status'} = "queued";
		$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'psGrepStr'} = "ps xu| grep -e \'$EdgeRCMD\' | grep -v grep";;
		$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'mCheckStr'} = "R --slave --no-restore --vanilla";
		$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'allinFileExist'} = "no";
		$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'alloutFileExist'} = "no";
		$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'runCMD'} = $EdgeRCMD;
		$processInfoHsh_ref->{$comparisonName}{"EdgeRCMD"}{'scrLog'} = $EdgeRScrLog;

		print CMDLOG $comparisonName."\t".'DESeq1CMD'."\t".$DESeq1CMD."\n";
		print CMDLOG $comparisonName."\t".'DESeq2CMD'."\t".$DESeq2CMD."\n";
		print CMDLOG $comparisonName."\t".'EdgeRCMD'."\t".$EdgeRCMD."\n";
	}
	close CMDLOG;
	
	my $allCMDDone = 'no';

	while ($allCMDDone ne 'yes') {

		foreach my $comparisonName (sort {$a cmp $b} keys %{$comparisonInfoHsh_ref}) {

			my $refSample = $comparisonInfoHsh_ref->{$comparisonName}{'refSample'};
			my $qrySample = $comparisonInfoHsh_ref->{$comparisonName}{'qrySample'};

			foreach my $CMD (keys %{$processInfoHsh_ref->{$comparisonName}}) {
				my $status = $processInfoHsh_ref->{$comparisonName}{$CMD}{'status'};
				
				#---check if CMD is running
				if ($status eq 'running') {
					my $psGrepStr = $processInfoHsh_ref->{$comparisonName}{$CMD}{"psGrepStr"};
					my $mCheckStr = $processInfoHsh_ref->{$comparisonName}{$CMD}{"mCheckStr"};
					my $theProcesses = `$psGrepStr`; #---grep the process
					if ($theProcesses !~ m/$mCheckStr/) {#---the process is not running
	
						#---check if all outFile are there
						$processInfoHsh_ref->{$comparisonName}{$CMD}{'allOutFileExist'} = "yes";
						foreach my $outFile (@{$processInfoHsh_ref->{$comparisonName}{$CMD}{'outFile'}}) {
							$processInfoHsh_ref->{$comparisonName}{$CMD}{'allOutFileExist'} = "no" if (not (-s $outFile));
						}
	
						if ($processInfoHsh_ref->{$comparisonName}{$CMD}{'allOutFileExist'} eq "yes") {#---all outputs are here: i.e. process finished
							$processInfoHsh_ref->{$comparisonName}{$CMD}{'status'} = "finished";

						} else {
							die "$CMD finished without generating all outputs\n";
						}
					}
				
				#---check if CMD is queued
				} elsif ($status eq 'queued') {
					
					#---check if all outFile are there
					$processInfoHsh_ref->{$comparisonName}{$CMD}{'allOutFileExist'} = "yes";
					foreach my $outFile (@{$processInfoHsh_ref->{$comparisonName}{$CMD}{'outFile'}}) {
						$processInfoHsh_ref->{$comparisonName}{$CMD}{'allOutFileExist'} = "no" if (not (-s $outFile));
					}
					
					if ($processInfoHsh_ref->{$comparisonName}{$CMD}{'allOutFileExist'} eq "yes") {#---all outputs are here: i.e. process finished, so skipp
						$processInfoHsh_ref->{$comparisonName}{$CMD}{'status'} = "skipped";
					
					} else {
					
						my $TSVToPrintPath = $processInfoHsh_ref->{$comparisonName}{$CMD}{'inTSV'};
						my $meanCountCutoff = my $minCountInOneSample = 0;
						my ($countInfoHsh_ref) = &getCountInput($refSample, $qrySample, $allSampleCountInfoHsh_ref);#->1156
						&printCountTSVInput($countInfoHsh_ref, $meanCountCutoff, $minCountInOneSample, $TSVToPrintPath, $refSample, $qrySample);#->2122

						#---check if all inFile are there
						$processInfoHsh_ref->{$comparisonName}{$CMD}{'allInFileExist'} = "yes";
						foreach my $inFile (@{$processInfoHsh_ref->{$comparisonName}{$CMD}{'inFile'}}) {
							$processInfoHsh_ref->{$comparisonName}{$CMD}{'allInFileExist'} = "no" if (not (-s $inFile));
						}
						
						
						if ($processInfoHsh_ref->{$comparisonName}{$CMD}{'allInFileExist'} eq "yes") {
							my $runCMD = $processInfoHsh_ref->{$comparisonName}{$CMD}{'runCMD'};
							my $scrLog = $processInfoHsh_ref->{$comparisonName}{$CMD}{'scrLog'};
							system ("rm -f $scrLog");
							system(qq|$runCMD 1>>$scrLog 2>>$scrLog &|);
							$processInfoHsh_ref->{$comparisonName}{$CMD}{'status'} = "running";
						}
					}

				} elsif ($status eq 'finished') {
				
				} elsif ($status eq 'aborted') {
				
				} elsif ($status eq 'skipped') {
				
				}
			}
		}

		$allCMDDone = 'yes';
		system ("clear");
		foreach my $comparisonName (sort {$a cmp $b} keys %{$comparisonInfoHsh_ref}) {
			foreach my $CMD (sort {$a cmp $b}  keys %{$processInfoHsh_ref->{$comparisonName}}) {
				my $status = $processInfoHsh_ref->{$comparisonName}{$CMD}{'status'};
				
				$allCMDDone = 'no' if (($status eq 'queued') or ($status eq 'running'));
				
				my $scrLog = $processInfoHsh_ref->{$comparisonName}{$CMD}{'scrLog'};
				my $lastLineScrLog = "no scrLog";
				$lastLineScrLog = `tail -1 $scrLog` if -s $scrLog;
				chomp $lastLineScrLog;
				print $comparisonName."\t".$CMD."\t".$status."\t".$lastLineScrLog."\n";
			}
		}
		sleep 5;
	}
	
	system "clear";
}
sub overlapGeneExpDiffResults {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: pairwiseVennDiagram|1696, runAndCheckSerialTask|2316, tripleVennDiagram|2338
#	appearInSub: runAllComparisons|2244
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_runDEAnalyses|131
#	input: $DESeq1FDR, $DESeq2FDR, $EdgeRFDR, $allComparisonFoldChangeHsh_ref, $allResultFilePathHsh_ref, $allSampleCountInfoHsh_ref, $combinedCutoffOvrlpHsh_ref, $combinedResultCountHsh_ref, $combinedResultInfoHsh_ref, $comparisonInfoHsh_ref, $geneNameLenHsh_ref, $minLinearFC, $minRPKM, $upDnAll
#	output: $combinedResultCountHsh_ref
#	toCall: my ($combinedResultCountHsh_ref) = &overlapGeneExpDiffResults($comparisonInfoHsh_ref, $DESeq1FDR, $DESeq2FDR, $EdgeRFDR, $minLinearFC, $upDnAll, $combinedResultInfoHsh_ref, $combinedResultCountHsh_ref, $combinedCutoffOvrlpHsh_ref, $allComparisonFoldChangeHsh_ref, $minRPKM, $geneNameLenHsh_ref, $allSampleCountInfoHsh_ref, $allResultFilePathHsh_ref);
#	calledInLine: 2302
#....................................................................................................................................................#

	my ($comparisonInfoHsh_ref, $DESeq1FDR, $DESeq2FDR, $EdgeRFDR, $minLinearFC, $upDnAll, $combinedResultInfoHsh_ref, $combinedResultCountHsh_ref, $combinedCutoffOvrlpHsh_ref, $allComparisonFoldChangeHsh_ref, $minRPKM, $geneNameLenHsh_ref, $allSampleCountInfoHsh_ref, $allResultFilePathHsh_ref) = @_;

	foreach my $comparisonName (sort {$a cmp $b} keys %{$comparisonInfoHsh_ref}) {
		my @allPdfToMerge;
		my $refSample = $comparisonInfoHsh_ref->{$comparisonName}{'refSample'};
		my $qrySample = $comparisonInfoHsh_ref->{$comparisonName}{'qrySample'};

		print "overlapping results of $comparisonName at $upDnAll FC $minLinearFC RPKM $minRPKM.\n";
		
		my $genePassRPKMHsh_ref = {};
		my $genePassDESeq1FCHsh_ref = {};
		my $genePresenceHsh_ref = {};
		
		############Get geneIDs that passed minRPKM
		foreach my $valueType (qw/fpkm fpm inputCount/) {
			foreach my $geneID (keys %{$geneNameLenHsh_ref}) {
				foreach my $sampleName ($refSample, $qrySample) {
					my @valueAry;
					foreach my $repNum (keys %{$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}}) {
						push @valueAry, sprintf "%.4f", $allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{$valueType};
					}
					my $avgValue = sum(@valueAry)/@valueAry;
					$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"Avg_$valueType\_ref"} = $avgValue if $sampleName eq $refSample;
					$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"Avg_$valueType\_qry"} = $avgValue if $sampleName eq $qrySample;
					if ($valueType eq 'fpkm') {
						$genePassRPKMHsh_ref->{$geneID}++ if ($avgValue >= $minRPKM);
					}
				}
			}
		}
		
		#############---Read DESeq1 results
		my $DESeq1DEGResultPath = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq1DEGResultPath'};

		my $DESeq1SigNum = 0;
		my %DESeq1SigGeneHsh;
		
		open (DESEQ1DE, "$DESeq1DEGResultPath");
		<DESEQ1DE>;

		while (<DESEQ1DE>) {
			chomp;
			$_ =~ s/\"//g;
			my ($num, $geneID, $baseMean, $value1, $value2, $foldChange, $log2FC, $pValue, $qValue) = split /\t/;
			$foldChange = 0 if $foldChange eq 'NA';

			$log2FC = -99999 if $log2FC eq "-Inf";
			$log2FC = 99999 if $log2FC eq "Inf";

			$log2FC = 0 if $log2FC eq 'NA';
			$pValue = 1 if $pValue eq 'NA';
			$qValue = 1 if $qValue eq 'NA';
			
			my $DESeq1_linearFC;
			
			if ($log2FC == -99999 or $log2FC == 99999) {
				$DESeq1_linearFC = $log2FC;
			} else {
				$DESeq1_linearFC = exp(abs($log2FC) * log (2));
				$DESeq1_linearFC *= -1 if $log2FC < 0;
			}

			$genePresenceHsh_ref->{"DESeq1"}{$geneID}++;

			$genePassDESeq1FCHsh_ref->{$geneID}++ if (abs $DESeq1_linearFC >= $minLinearFC);

			if (($qValue <= $DESeq1FDR) and (exists $genePassDESeq1FCHsh_ref->{$geneID}) and (exists $genePassRPKMHsh_ref->{$geneID})){
				if (($upDnAll eq "all")
				or (($upDnAll eq "up") and ($log2FC > 0))
				or (($upDnAll eq "dn") and ($log2FC < 0))) {
					$DESeq1SigNum++;
					$DESeq1SigGeneHsh{$geneID} = $log2FC;
				}
			}

			#----take the chosen DESeq model data
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_Value1"} = $value1;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_Value2"} = $value2;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_Log2FC"} = $log2FC;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_pValue"} = $pValue;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_qValue"} = $qValue;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_linearFC"} = sprintf "%.3f", $DESeq1_linearFC;
		}
		close DESEQ1DE;

		#---get moderate FC from DESeq1
		my $DESeq1ModLog2FCCSVPath = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq1ModLog2FCCSVPath'};
		open (DESEQMFC, "$DESeq1ModLog2FCCSVPath");
		<DESEQMFC>;
		while (<DESEQMFC>) {
			chomp;
			$_ =~ s/\"//g;
			my ($geneID, $modLog2FC) = split /\,/;
			$modLog2FC = 0 if $modLog2FC eq 'NA';
			#----take the chosen DESeq model data
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_modLog2FC"} = $modLog2FC;
		}
		close DESEQMFC;

		#############---Read DESeq2 results
		my $DESeq2DEGResultPath = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq2DEGResultPath'};

		my $DESeq2SigNum = 0;
		my %DESeq2SigGeneHsh;
		
		open (DESEQ2DE, "$DESeq2DEGResultPath");
		<DESEQ2DE>;

		while (<DESEQ2DE>) {
			chomp;
			$_ =~ s/\"//g;
			my ($geneID, $baseMean, $log2FoldChange, $pvalue, $FDR) = split /\t/;
			$log2FoldChange = 0 if $log2FoldChange eq 'NA';
			$pvalue = 1 if $pvalue eq 'NA';
			$FDR = 1 if $FDR eq 'NA';

			my $DESeq2_linearFC = exp(abs($log2FoldChange) * log (2));
			$DESeq2_linearFC *= -1 if $log2FoldChange < 0;

			$genePresenceHsh_ref->{"DESeq2"}{$geneID}++;

			if (($FDR <= $DESeq2FDR) and (exists $genePassDESeq1FCHsh_ref->{$geneID}) and (exists $genePassRPKMHsh_ref->{$geneID})) {
				if (($upDnAll eq "all")
				or (($upDnAll eq "up") and ($combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_Log2FC"} > 0))
				or (($upDnAll eq "dn") and ($combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_Log2FC"} < 0))) {
					$DESeq2SigNum++;
					$DESeq2SigGeneHsh{$geneID} = $log2FoldChange;
				}
			}

			#----take the chosen DESeq model data
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq2_Log2FC"} = $log2FoldChange;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq2_pValue"} = $pvalue;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq2_qValue"} = $FDR;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq2_linearFC"} = sprintf "%.3f", $DESeq2_linearFC;
		}
		close DESEQ2DE;

		#############---Read EdgeR results
		my $EdgeRGenesExactTestTSVPath = $comparisonInfoHsh_ref->{$comparisonName}{'EdgeRGenesExactTestTSVPath'};
		my $EdgeRExpressionValuesTSVPath = $comparisonInfoHsh_ref->{$comparisonName}{'EdgeRExpressionValuesTSVPath'};
		open (EDGERET, "$EdgeRGenesExactTestTSVPath");
		my $EdgeRHeader = <EDGERET>;
		my $EdgeRSigNum = 0;
		my %EdgeRSigGeneHsh;
		while (<EDGERET>) {
			chomp;
			$_ =~ s/\"//g;
			my ($geneID, $log2FC, $logCPM, $pValue, $qValue) = split /\t/;

			$genePresenceHsh_ref->{"EdgeR"}{$geneID}++;

			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"EdgeR_logCPM"} = $logCPM;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"EdgeR_log2FC"} = $log2FC;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"EdgeR_pValue"} = $pValue;
			$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"EdgeR_qValue"} = $qValue;

			if (($qValue <= $EdgeRFDR) and (exists $genePassDESeq1FCHsh_ref->{$geneID}) and (exists $genePassRPKMHsh_ref->{$geneID})) {
				if (($upDnAll eq "all")
				or (($upDnAll eq "up") and ($log2FC > 0))
				or (($upDnAll eq "dn") and ($log2FC < 0))) {
					$EdgeRSigNum++;
					$EdgeRSigGeneHsh{$geneID} = $log2FC;
				}
			}
		}
		close EDGERET;
		
		#----fill the empty geneIDs (those were not tested in DESeq but in the gene list. it happens when cuffdiff minimum read cutoff, these genes will not appear in the DeSeq input);
		foreach my $geneID (keys %{$geneNameLenHsh_ref}) {
			if (not exists $genePresenceHsh_ref->{"EdgeR"}{$geneID}) {
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"EdgeR_logCPM"} = -1;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"EdgeR_log2FC"} = 0;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"EdgeR_pValue"} = 1;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"EdgeR_qValue"} = 1;
			}

			if (not exists $genePresenceHsh_ref->{"DESeq1"}{$geneID}) {
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_Value1"} = -1;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_Value2"} = -1;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_modLog2FC"} = 0;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_Log2FC"} = 0;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_linearFC"} = 0;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_pValue"} = 1;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_qValue"} = 1;
			}

			if (not exists $genePresenceHsh_ref->{"DESeq2"}{$geneID}) {
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq2_Log2FC"} = 0;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq2_pValue"} = 1;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq2_qValue"} = 1;
				$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq2_linearFC"} = 0;
			}
		}
		
		my $DESeq1OutDir = $comparisonInfoHsh_ref->{$comparisonName}{'DESeq1OutDir'};

		system "mkdir -p -m 777 $DESeq1OutDir/overlap/FC$minLinearFC/";
		my $D2D1PdfPathPrefix = "$DESeq1OutDir/overlap/FC$minLinearFC/$upDnAll.RPKM.$minRPKM.FC$minLinearFC.ds2.$DESeq2FDR.vs.ds1.$DESeq1FDR.$comparisonName";
		my ($CDOvrlpHsh_ref, $CDMergeHsh_ref, undef, undef, $CDMergeNum, $CDOvrlpNum, undef, undef, $D2D1PdfPath) = &pairwiseVennDiagram(\%DESeq2SigGeneHsh, \%DESeq1SigGeneHsh, 'DESeq2', "DESeq1", $D2D1PdfPathPrefix, $D2D1PdfPathPrefix.".log.txt", 'red', 'blue');#->1696
		push @allPdfToMerge, $D2D1PdfPath;

		my $ErD1PdfPathPrefix = "$DESeq1OutDir/overlap/FC$minLinearFC/$upDnAll.RPKM.$minRPKM.FC$minLinearFC.edg.$EdgeRFDR.vs.ds1.$DESeq1FDR.$comparisonName";
		my ($EDOvrlpHsh_ref, $EDMergeHsh_ref, undef, undef, $EDMergeNum, $EDOvrlpNum, undef, undef, $ErD1PdfPath) = &pairwiseVennDiagram(\%EdgeRSigGeneHsh, \%DESeq1SigGeneHsh, 'EdgeR', "DESeq1", $ErD1PdfPathPrefix, $ErD1PdfPathPrefix.".log.txt", 'green', 'blue');#->1696
		push @allPdfToMerge, $ErD1PdfPath;

		my $ErD2PdfPathPrefix = "$DESeq1OutDir/overlap/FC$minLinearFC/$upDnAll.RPKM.$minRPKM.FC$minLinearFC.edg.$EdgeRFDR.vs.ds2.$DESeq2FDR.$comparisonName";
		my ($ECOvrlpHsh_ref, $ECMergeHsh_ref, undef, undef, $ECMergeNum, $ECOvrlpNum, undef, undef, $ErD2PdfPath) = &pairwiseVennDiagram(\%EdgeRSigGeneHsh, \%DESeq2SigGeneHsh, 'EdgeR', "DESeq2", $ErD2PdfPathPrefix, $ErD2PdfPathPrefix.".log.txt", 'green', 'red');#->1696
		push @allPdfToMerge, $ErD2PdfPath;

		my $D2D1ErPdfPathPrefix = "$DESeq1OutDir/overlap/FC$minLinearFC/$upDnAll.RPKM.$minRPKM.FC$minLinearFC.ds2.$DESeq2FDR.vs.ds1.$DESeq1FDR.vs.edg.$EdgeRFDR.$comparisonName";
		my (undef, undef, undef, $CDEOvrlpHsh_ref, $CDEMergeHsh_ref, $D2D1ErPdfPath) = &tripleVennDiagram(\%DESeq2SigGeneHsh, \%DESeq1SigGeneHsh, \%EdgeRSigGeneHsh, 'DESeq2', "DESeq1", 'EdgeR', $D2D1ErPdfPathPrefix, $D2D1ErPdfPathPrefix.".log.txt", 'red', 'blue', 'green');#->2338
		push @allPdfToMerge, $D2D1ErPdfPath;
		
		my %combinationCountHsh;
		$combinationCountHsh{'ovrlp_1_any'} = 0;
		$combinationCountHsh{'ovrlp_2_any'} = 0;
		$combinationCountHsh{'ovrlp_3_all'} = 0;
		$combinationCountHsh{'ovrlp_orphan'} = 0;
		
		foreach my $geneID (keys %{$geneNameLenHsh_ref}) {
			
			my $significantNum = 0;
			
			$significantNum++ if exists $EdgeRSigGeneHsh{$geneID};
			$significantNum++ if exists $DESeq2SigGeneHsh{$geneID};
			$significantNum++ if exists $DESeq1SigGeneHsh{$geneID};
			
			$combinationCountHsh{'ovrlp_1_any'}++ if $significantNum >= 1;
			$combinationCountHsh{'ovrlp_2_any'}++ if $significantNum >= 2;
			$combinationCountHsh{'ovrlp_3_all'}++ if $significantNum == 3;
			$combinationCountHsh{'ovrlp_orphan'}++ if $significantNum == 1;

		}
		
		#----store boolean for all data
		if (($upDnAll eq 'all') and ($minLinearFC == 0)) {
			foreach my $geneID (keys %{$geneNameLenHsh_ref}) {
				
				my @significantAry;
				
				push @significantAry, 'Er' if exists $EdgeRSigGeneHsh{$geneID};
				push @significantAry, 'D2' if exists $DESeq2SigGeneHsh{$geneID};
				push @significantAry, 'D1' if exists $DESeq1SigGeneHsh{$geneID};
				
				$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{'ovrlp_1_any'} = 0;
				$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{'ovrlp_2_any'} = 0;
				$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{'ovrlp_3_all'} = 0;
				$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{'ovrlp_orphan'} = 0;

				$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{'ovrlp_1_any'} = join "", @significantAry if @significantAry >= 1;
				$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{'ovrlp_2_any'} = join "", @significantAry if @significantAry >= 2;;
				$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{'ovrlp_3_all'} = join "", @significantAry if @significantAry == 3;;
				$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{'ovrlp_orphan'} = join "", @significantAry if @significantAry == 1;;

				foreach my $FDRCutoffCombination ((qw/ovrlp_1_any ovrlp_2_any ovrlp_3_all ovrlp_orphan/)) {
					$allComparisonFoldChangeHsh_ref->{$FDRCutoffCombination}{$geneID}{$comparisonName} = 'na';
					$allComparisonFoldChangeHsh_ref->{$FDRCutoffCombination}{$geneID}{$comparisonName} = $combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{"DESeq1_linearFC"} if $combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{$FDRCutoffCombination} ne '0';
				}
			}
		}
		
		$combinedResultCountHsh_ref->{"RPKM.$minRPKM\tFC.$minLinearFC\tDESeq1\tFDR.$DESeq1FDR"}{"$upDnAll.$comparisonName"} = $DESeq1SigNum;
		$combinedResultCountHsh_ref->{"RPKM.$minRPKM\tFC.$minLinearFC\tDESeq2\tFDR.$DESeq2FDR"}{"$upDnAll.$comparisonName"} = $DESeq2SigNum;
		$combinedResultCountHsh_ref->{"RPKM.$minRPKM\tFC.$minLinearFC\tEdgeR\tFDR.$EdgeRFDR"}{"$upDnAll.$comparisonName"} = $EdgeRSigNum;
		$combinedResultCountHsh_ref->{"RPKM.$minRPKM\tFC.$minLinearFC\tovrlp_1_any\tD2_.$DESeq2FDR.D1_.$DESeq1FDR.Er_.$EdgeRFDR"}{"$upDnAll.$comparisonName"} = $combinationCountHsh{'ovrlp_1_any'};
		$combinedResultCountHsh_ref->{"RPKM.$minRPKM\tFC.$minLinearFC\tovrlp_2_any\tD2_.$DESeq2FDR.D1_.$DESeq1FDR.Er_.$EdgeRFDR"}{"$upDnAll.$comparisonName"} = $combinationCountHsh{'ovrlp_2_any'};
		$combinedResultCountHsh_ref->{"RPKM.$minRPKM\tFC.$minLinearFC\tovrlp_3_all\tD2_.$DESeq2FDR.D1_.$DESeq1FDR.Er_.$EdgeRFDR"}{"$upDnAll.$comparisonName"} = $combinationCountHsh{'ovrlp_3_all'};
		$combinedResultCountHsh_ref->{"RPKM.$minRPKM\tFC.$minLinearFC\tovrlp_orphan\tD2_.$DESeq2FDR.D1_.$DESeq1FDR.Er_.$EdgeRFDR"}{"$upDnAll.$comparisonName"} = $combinationCountHsh{'ovrlp_orphan'};

		#---merge all pdfs
		my $rootOutDir = $comparisonInfoHsh_ref->{$comparisonName}{'rootOutDir'};
		print "Merging all pdfs in $comparisonName\n";
		system "mkdir -p -m 777 $rootOutDir/overlap/pdf/FC$minLinearFC/";
		my $mergErD1D2PdfPath = "$rootOutDir/overlap/pdf/FC$minLinearFC/$upDnAll.RPKM.$minRPKM.FC$minLinearFC.D1_.$DESeq1FDR.D2_.$DESeq2FDR.Er_.$EdgeRFDR.$comparisonName.pdf";
		my $allPdfPath = join " ", @allPdfToMerge;
		my $mergePdfCMD = "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=\'$mergErD1D2PdfPath\' $allPdfPath";
		my $grepCMD = "ps -ef | grep $allPdfPath | grep -v grep";
		&runAndCheckSerialTask($grepCMD, $mergePdfCMD, $mergePdfCMD, "/dev/null");#->2316

		$allResultFilePathHsh_ref->{"overlapBetweenMethodsPDF"}{$upDnAll."RegulatedGenes"}{"cutoff.RPKM.$minRPKM.FC$minLinearFC"}{"$comparisonName Er[$EdgeRSigNum] D1[$DESeq1SigNum] D2[$DESeq2SigNum]"} = $mergErD1D2PdfPath;
	}
	
	return $combinedResultCountHsh_ref;
}
sub pairwiseVennDiagram {
#....................................................................................................................................................#
#	subroutineCategory: plotInR
#	dependOnSub: >none
#	appearInSub: overlapGeneExpDiffResults|1398
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $RScrLog, $color1, $color2, $listAHsh_ref, $listAName, $listBHsh_ref, $listBName, $pdfPathPrefix
#	output: $ABMergeNum, $ABOvrlpNum, $AOnlyNum, $BOnlyNum, $pdfPath, \%ABMergeHsh, \%ABOvrlpHsh, \%AOnlyHsh, \%BOnlyHsh
#	toCall: my (\%ABOvrlpHsh, \%ABMergeHsh, \%AOnlyHsh, \%BOnlyHsh, $ABOvrlpNum, $ABMergeNum, $AOnlyNum, $BOnlyNum, $pdfPath) = &pairwiseVennDiagram($listAHsh_ref, $listBHsh_ref, $listAName, $listBName, $pdfPathPrefix, $RScrLog, $color1, $color2);
#	calledInLine: 1609, 1613, 1617
#....................................................................................................................................................#

	my ($listAHsh_ref, $listBHsh_ref, $listAName, $listBName, $pdfPathPrefix, $RScrLog, $color1, $color2) = @_;
	
	my $listANum = keys %{$listAHsh_ref};
	my $listBNum = keys %{$listBHsh_ref};
	
	my %ABOvrlpHsh;
	my %ABMergeHsh;
	my %AOnlyHsh = %{$listAHsh_ref};
	my %BOnlyHsh = %{$listBHsh_ref};
	
	my $AOnlyNum = my $BOnlyNum = my $ABOvrlpNum = my $ABMergeNum = -1;
	my $pdfPath = "";

	if (($listANum > 0) and ($listBNum > 0)) {
		foreach my $id (keys %{$listAHsh_ref}) {
			if (exists $listBHsh_ref->{$id}) {
				$ABOvrlpHsh{$id}++;
				delete $BOnlyHsh{$id};
				delete $AOnlyHsh{$id};
			}
		}
		
		my ($name1, $name2) = ($listAName."[$listANum]", $listBName."[$listBNum]");
		
		if ($listBNum > $listANum) {
			($name1, $name2) = ($name2, $name1);
			($color1, $color2) = ($color2, $color1);
		}
		
		
		foreach my $id (keys %ABOvrlpHsh) {$ABMergeHsh{$id} = $ABOvrlpHsh{$id};}
		foreach my $id (keys %BOnlyHsh) {$ABMergeHsh{$id} = $BOnlyHsh{$id};}
		foreach my $id (keys %AOnlyHsh) {$ABMergeHsh{$id} = $AOnlyHsh{$id};}

		$AOnlyNum = keys %AOnlyHsh;
		$BOnlyNum = keys %BOnlyHsh;
		$ABOvrlpNum = keys %ABOvrlpHsh;
		$ABMergeNum = keys %ABMergeHsh;
		
		$pdfPath = "$pdfPathPrefix.A.$AOnlyNum.B.$BOnlyNum.Ovrlp.$ABOvrlpNum.Merge.$ABMergeNum.pdf";
		open (RSCRLOG, ">$RScrLog");
		print RSCRLOG "########## R CMD START #############\n";
		print RSCRLOG "library(VennDiagram)\n";
		print RSCRLOG "venn.plot <- draw.pairwise.venn($listANum, $listBNum, $ABOvrlpNum, category = c(\"$name1\", \"$name2\"), fill = c(\"$color1\", \"$color2\"), lty = \"blank\", cex = 2, cat.cex = 2, cat.pos = c(0, 180), cat.dist = 0.05, cat.just = list(c(0.5, 2), c(0.5, 0)), cat.col = c(\"$color1\", \"$color2\"), ext.pos = 30, ext.dist = -0.05, ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = \"dashed\" )\n";
		print RSCRLOG "pdf(\"$pdfPath\")\n";
		print RSCRLOG "grid.draw(venn.plot)\n";
		print RSCRLOG "dev.off()\n";
		print RSCRLOG "########## R CMD END #############\n";
		close RSCRLOG;

		open R, "|/usr/bin/R --vanilla --slave 1>>$RScrLog 2>>$RScrLog";
		print R "library(VennDiagram)\n";
		print R "venn.plot <- draw.pairwise.venn( area1 = $listANum, area2 = $listBNum, cross.area = $ABOvrlpNum, category = c(\"$name1\", \"$name2\"), fill = c(\"$color1\", \"$color2\"), lty = \"blank\", cex = 2, cat.cex = 2, cat.pos = c(0, 180), cat.dist = 0.05, cat.just = list(c(0.5, 2), c(0.5, 0)), cat.col = c(\"$color1\", \"$color2\"), ext.pos = 30, ext.dist = -0.05, ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = \"dashed\" )\n";
		print R "pdf(\"$pdfPath\")\n";
		print R "grid.draw(venn.plot)\n";
		print R "dev.off()\n";
		close R;
	}
	
	return (\%ABOvrlpHsh, \%ABMergeHsh, \%AOnlyHsh, \%BOnlyHsh, $ABOvrlpNum, $ABMergeNum, $AOnlyNum, $BOnlyNum, $pdfPath);
}
sub printAllComparisonFoldChange {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: runAllComparisons|2244
#	primaryAppearInSection: 5_summarizeResults|141
#	secondaryAppearInSection: 4_runDEAnalyses|131
#	input: $FDRTag, $allComparisonFoldChangeHsh_ref, $comparisonToPrintHsh_ref, $dirToPrint, $filePrefix, $geneNameLenHsh_ref
#	output: none
#	toCall: &printAllComparisonFoldChange($allComparisonFoldChangeHsh_ref, $geneNameLenHsh_ref, $dirToPrint, $filePrefix, $comparisonToPrintHsh_ref, $FDRTag);
#	calledInLine: 147, 2308
#....................................................................................................................................................#

	my ($allComparisonFoldChangeHsh_ref, $geneNameLenHsh_ref, $dirToPrint, $filePrefix, $comparisonToPrintHsh_ref, $FDRTag) = @_;
	
	system ("mkdir -pm 777 $dirToPrint");
	
	foreach my $FDRCutoffCombination (keys %{$allComparisonFoldChangeHsh_ref}) {
		open (COMBFC, ">$dirToPrint/$filePrefix.combined.FC.$FDRCutoffCombination.$FDRTag.xls");

		#----print header
		my $comparisonNum = 0;
		foreach my $geneID (sort keys %{$allComparisonFoldChangeHsh_ref->{$FDRCutoffCombination}}) {
			my @outputAry = ();
			push @outputAry, 'geneID';
			push @outputAry, 'name';
			foreach my $comparisonName (sort keys %{$allComparisonFoldChangeHsh_ref->{$FDRCutoffCombination}{$geneID}}) {
				next unless (exists $comparisonToPrintHsh_ref->{$comparisonName});
				push @outputAry, $comparisonName;
				$comparisonNum++;
			}
			push @outputAry, (qw/upNum dnNum bothNum overall avgFC/);
			print COMBFC join "", ((join "\t", @outputAry),"\n");
			last;
		}

		#----print content
		foreach my $geneID (sort keys %{$allComparisonFoldChangeHsh_ref->{$FDRCutoffCombination}}) {
			my @outputAry = ();
			push @outputAry, $geneID;
			push @outputAry, $geneNameLenHsh_ref->{$geneID}{'name'};
			my $upNum = my $dnNum = my $bothNum = 0;
			my @allFCAry;
			my $comparisonNum = 0;
			foreach my $comparisonName (sort keys %{$allComparisonFoldChangeHsh_ref->{$FDRCutoffCombination}{$geneID}}) {
				next unless (exists $comparisonToPrintHsh_ref->{$comparisonName});
				$comparisonNum++;
				my $DESeq1_linearFC = $allComparisonFoldChangeHsh_ref->{$FDRCutoffCombination}{$geneID}{$comparisonName};
				push @outputAry, $DESeq1_linearFC;
				if (looks_like_number($DESeq1_linearFC) and $comparisonName !~ m/pooledOthers/) {
					push @allFCAry, $DESeq1_linearFC;
					$bothNum++;
					$upNum++ if $DESeq1_linearFC > 0;
					$dnNum++ if $DESeq1_linearFC < 0;
				}
			}
			my $avgFC = 'na';
			$avgFC = sprintf "%.3f", sum(@allFCAry)/@allFCAry if @allFCAry > 0;
			my $overall;
			if ($upNum > $comparisonNum/2 and $dnNum == 0) {
				$overall = 'overallUp';
			} elsif ($dnNum > $comparisonNum/2 and $upNum == 0) {
				$overall = 'overallDown';
			} else {
				$overall = 'inconsistent';
			}
			push @outputAry, ($upNum, $dnNum, $bothNum, $overall, $avgFC);
			print COMBFC join "", ((join "\t", @outputAry),"\n");
		}
		
		close COMBFC;
	}
}
sub printAllSampleCountInfo {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 5_summarizeResults|141
#	secondaryAppearInSection: >none
#	input: $allSampleCountInfoHsh_ref, $outDir
#	output: $sampleTotalFragNumHsh_ref, \%allSampleCountPathHsh, \@sampleColumnAry
#	toCall: my (\%allSampleCountPathHsh, $sampleTotalFragNumHsh_ref, \@sampleColumnAry) = &printAllSampleCountInfo($allSampleCountInfoHsh_ref, $outDir);
#	calledInLine: 149, 1855
#....................................................................................................................................................#
	
	#---my ($allSampleCountPathHsh_ref, $sampleTotalFragNumHsh_ref, $sampleColumnAry_ref) = &printAllSampleCountInfo(\%allSampleCountInfoHsh, $outDir);#->1843

	my ($allSampleCountInfoHsh_ref, $outDir) = @_;
	
	my $sampleTotalFragNumHsh_ref = {};
	my %allSampleCountPathHsh;
	my @sampleColumnAry;
	
	system ('clear');
	print "Printing count data for all samples\n";
	
	#----print the values
	my $pushSampleName = 'yes';
	system "mkdir -p -m 777 $outDir/allSample/count/";
	
	foreach my $log2OrLinear ('log2', 'linear') {
		foreach my $valueType ('inputCount', 'fpm', 'fpkm') {
			my $outPath = "$outDir/allSample/count/allSampleCount.$valueType.$log2OrLinear.xls";
			${$allSampleCountPathHsh{$valueType}}{$log2OrLinear} = $outPath;
			open (OUTFILE, ">$outPath");
			foreach my $geneID (sort {$a cmp $b} keys %{$allSampleCountInfoHsh_ref}) {
				my @outputAry;
				push @outputAry, 'geneID';
				foreach my $sampleName (sort {$a cmp $b} keys %{$allSampleCountInfoHsh_ref->{$geneID}}) {
					foreach my $repNum (sort {$a <=> $b} keys %{$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}}) {
						push @outputAry, $sampleName."_".$repNum;
					}
				}
				print OUTFILE join "", ((join "\t", @outputAry), "\n");
				last;
			}
			
			foreach my $geneID (sort {$a cmp $b} keys %{$allSampleCountInfoHsh_ref}) {
				my @outputAry;
				push @outputAry, $geneID;
				foreach my $sampleName (sort {$a cmp $b} keys %{$allSampleCountInfoHsh_ref->{$geneID}}) {
					foreach my $repNum (sort {$a <=> $b} keys %{$allSampleCountInfoHsh_ref->{$geneID}{$sampleName}}) {
						
						push @sampleColumnAry, $sampleName if ($pushSampleName eq 'yes');
						
						my $value = $allSampleCountInfoHsh_ref->{$geneID}{$sampleName}{$repNum}{$valueType};

						if ($log2OrLinear eq 'log2') {
							if ($value > 0) {
								$value = log($value)/log(2);
							} else {
								$value = -9999;
							}
						} 
						push @outputAry, $value;
					}
				}
				print OUTFILE join "", ((join "\t", @outputAry), "\n");
				$pushSampleName = 'no' if @sampleColumnAry > 0;
			}
			close OUTFILE;
		}
	}
	
	return \%allSampleCountPathHsh, $sampleTotalFragNumHsh_ref, \@sampleColumnAry;
	
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|355
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|81, 8_finishingTasks|175
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 87, 180
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->355
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->355
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->355
		print "=========================================================================\n\n";
	}
}
sub printCombineResults {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: runAllComparisons|2244
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_runDEAnalyses|131
#	input: $FDRTag, $allResultFilePathHsh_ref, $allResultList_FH, $allSampleCountInfoHsh_ref, $combinedCutoffOvrlpHsh_ref, $combinedResultCountHsh_ref, $combinedResultInfoHsh_ref, $comparisonInfoHsh_ref, $dirToPrint, $geneNameLenHsh_ref, $overallCountFilePathHsh_ref, $refName
#	output: none
#	toCall: &printCombineResults($comparisonInfoHsh_ref, $combinedResultCountHsh_ref, $combinedResultInfoHsh_ref, $combinedCutoffOvrlpHsh_ref, $geneNameLenHsh_ref, $dirToPrint, $FDRTag, $allResultList_FH, $allResultFilePathHsh_ref, $allSampleCountInfoHsh_ref, $refName, $overallCountFilePathHsh_ref);
#	calledInLine: 2307
#....................................................................................................................................................#
	
	my ($comparisonInfoHsh_ref, $combinedResultCountHsh_ref, $combinedResultInfoHsh_ref, $combinedCutoffOvrlpHsh_ref, $geneNameLenHsh_ref, $dirToPrint, $FDRTag, $allResultList_FH, $allResultFilePathHsh_ref, $allSampleCountInfoHsh_ref, $refName, $overallCountFilePathHsh_ref) = @_;
	
	system ("mkdir -pm 777 $dirToPrint/overallCount/");
	my $overallCountPath = "$dirToPrint/overallCount/overallCount.$FDRTag.xls";
	open (COUNT, ">", $overallCountPath);
	
	$overallCountFilePathHsh_ref->{$FDRTag}{$refName} = $overallCountPath;
	
	foreach my $cutoffType (sort {$a cmp $b} keys %{$combinedResultCountHsh_ref}) {
		my @outputAry;
		push @outputAry, "minRPKM\tminFoldChange\tsoftware\tFDR";
		foreach my $countType (sort {$a cmp $b} keys %{$combinedResultCountHsh_ref->{$cutoffType}}) {
			push @outputAry, $countType;
		}
		print COUNT join "", ((join "\t", @outputAry), "\n");
		last;
	}

	foreach my $cutoffType (sort {$a cmp $b} keys %{$combinedResultCountHsh_ref}) {
		my @outputAry;
		push @outputAry, $cutoffType;
		foreach my $countType (sort {$a cmp $b} keys %{$combinedResultCountHsh_ref->{$cutoffType}}) {
			push @outputAry, $combinedResultCountHsh_ref->{$cutoffType}{$countType};
		}
		print COUNT join "", ((join "\t", @outputAry), "\n");
	}
	close COUNT;
	
	my %abstractItemToPrintHsh = (
		'Avg_fpkm_qry' => 1,
		'Avg_fpkm_ref' => 1,
		'Avg_inputCount_qry' => 1,
		'Avg_inputCount_ref' => 1,
		'DESeq1_linearFC' => 1,
		'DESeq1_Log2FC' => 1,
		'DESeq1_modLog2FC' => 1,
		'DESeq1_qValue' => 1,
		'DESeq2_qValue' => 1,
		'EdgeR_qValue' => 1,
		'ovrlp_2_any' => 1,
	);
	
	foreach my $fullOrAbstract (qw/full abstract/) {
		foreach my $comparisonName (sort {$a cmp $b} keys %{$combinedResultInfoHsh_ref}) {
			my $regularResultInfoPath = "$dirToPrint/$fullOrAbstract/$comparisonName.$fullOrAbstract.result.$FDRTag.regular.xls";
			my $saSameRowResultInfoPath = "$dirToPrint/$fullOrAbstract/$comparisonName.$fullOrAbstract.result.$FDRTag.senseAntisenseSameRow.xls";
			system ("mkdir -pm 777 $dirToPrint/$fullOrAbstract/");
		
			$allResultFilePathHsh_ref->{$fullOrAbstract.'regular'.'ResultTable'}{$FDRTag}{$comparisonName} = $regularResultInfoPath;
			$allResultFilePathHsh_ref->{$fullOrAbstract.'senseAntisenseSameRow'.'ResultTable'}{$FDRTag}{$comparisonName} = $saSameRowResultInfoPath;

			my $refSample = $comparisonInfoHsh_ref->{$comparisonName}{'refSample'};
			my $qrySample = $comparisonInfoHsh_ref->{$comparisonName}{'qrySample'};

			if ($fullOrAbstract eq 'full') {
				print {$allResultList_FH} join "", ((join "\t", ($comparisonName, $refSample, $qrySample, $regularResultInfoPath)), "\n");
			}
		
			open (SASAMEROW, ">", "$saSameRowResultInfoPath");
			open (REGULAR, ">", "$regularResultInfoPath");
			foreach my $geneID (sort {$a cmp $b} keys %{$allSampleCountInfoHsh_ref}) {
				my @regularOutputAry;
				my @saSameRowOutputAry;
				
				push @regularOutputAry, "geneID";
				push @regularOutputAry, "geneName";
				push @regularOutputAry, "locationTag";
				push @regularOutputAry, "geneLen";

				push @saSameRowOutputAry, "geneID";
				push @saSameRowOutputAry, "geneName";
				push @saSameRowOutputAry, "locationTag";
				push @saSameRowOutputAry, "geneLen";
				
				foreach my $infoType (sort {$a cmp $b} keys %{$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}}) {
					if (($fullOrAbstract eq 'full') or ($fullOrAbstract eq 'abstract' and $abstractItemToPrintHsh{$infoType})) {
						$infoType =~ s/_ref$/_$refSample/;
						$infoType =~ s/_qry$/_$qrySample/;
						push @regularOutputAry, $infoType;

						foreach my $dirtn (qw/a s/) {
							push @saSameRowOutputAry, $dirtn."_".$infoType;
						}
					}
				}
				foreach my $cutoffBooleanType (sort {$a cmp $b} keys %{$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}}) {
					if (($fullOrAbstract eq 'full') or ($fullOrAbstract eq 'abstract' and $abstractItemToPrintHsh{$cutoffBooleanType})) {
						push @regularOutputAry, $cutoffBooleanType;

						foreach my $dirtn (qw/a s/) {
							push @saSameRowOutputAry, $dirtn."_".$cutoffBooleanType;
						}
					}
				}
				print REGULAR join "", ((join "\t", @regularOutputAry), "\n");
				print SASAMEROW join "", ((join "\t", @saSameRowOutputAry), "\n");
				last;
			}

			my %printedGeneIDNoDirtnHsh = ();
			
			foreach my $geneID (sort {$a cmp $b} keys %{$allSampleCountInfoHsh_ref}) {
				
				#---[07/11/2013 15:39] see whether the geneID contains the a_ and s_ tag for sense anf antisense
				my $geneIDNoDirtn = $geneID;
				$geneIDNoDirtn =~ s/^a|s_//;
				my $printSaSameRow = 'no';
				if (exists $allSampleCountInfoHsh_ref->{"a_".$geneIDNoDirtn} and exists $allSampleCountInfoHsh_ref->{"a_".$geneIDNoDirtn} and not exists $printedGeneIDNoDirtnHsh{$geneIDNoDirtn}) {
					$printSaSameRow = 'yes';
					$printedGeneIDNoDirtnHsh{$geneIDNoDirtn}++;
				}
				my @regularOutputAry;
				push @regularOutputAry, $geneID;
				push @regularOutputAry, $geneNameLenHsh_ref->{$geneID}{'name'};
				push @regularOutputAry, $geneNameLenHsh_ref->{$geneID}{'locationTag'};
				push @regularOutputAry, $geneNameLenHsh_ref->{$geneID}{'len'};

				my @saSameRowOutputAry;
				if ($printSaSameRow eq 'yes') {
					push @saSameRowOutputAry, $geneIDNoDirtn;
					push @saSameRowOutputAry, $geneNameLenHsh_ref->{$geneID}{'name'};
					push @saSameRowOutputAry, $geneNameLenHsh_ref->{$geneID}{'locationTag'};
					push @saSameRowOutputAry, $geneNameLenHsh_ref->{$geneID}{'len'};
				}

				foreach my $infoType (sort {$a cmp $b} keys %{$combinedResultInfoHsh_ref->{$comparisonName}{$geneID}}) {
					if (($fullOrAbstract eq 'full') or ($fullOrAbstract eq 'abstract' and $abstractItemToPrintHsh{$infoType})) {
						push @regularOutputAry, $combinedResultInfoHsh_ref->{$comparisonName}{$geneID}{$infoType};

						if ($printSaSameRow eq 'yes') {
							foreach my $dirtn (qw/a s/) {
								push @saSameRowOutputAry, $combinedResultInfoHsh_ref->{$comparisonName}{$dirtn."_".$geneIDNoDirtn}{$infoType};
							}
						}
					}
				}
				foreach my $cutoffBooleanType (sort {$a cmp $b} keys %{$combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}}) {
					if (($fullOrAbstract eq 'full') or ($fullOrAbstract eq 'abstract' and $abstractItemToPrintHsh{$cutoffBooleanType})) {
						push @regularOutputAry, $combinedCutoffOvrlpHsh_ref->{$comparisonName}{$geneID}{$cutoffBooleanType};

						if ($printSaSameRow eq 'yes') {
							foreach my $dirtn (qw/a s/) {
								push @saSameRowOutputAry, $combinedCutoffOvrlpHsh_ref->{$comparisonName}{$dirtn."_".$geneIDNoDirtn}{$cutoffBooleanType};
							}
						}
					}
				}
				print REGULAR join "", ((join "\t", @regularOutputAry), "\n");

				if ($printSaSameRow eq 'yes') {
					print SASAMEROW join "", ((join "\t", @saSameRowOutputAry), "\n");
				}
			}
			close REGULAR;
			close SASAMEROW;
		}
	}
}
sub printCountTSVInput {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: launchAndMonitorAllCMD|1219
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $TSVToPrintPath, $countInfoHsh_ref, $meanCountCutoff, $minCountInOneSample, $qrySample, $refSample
#	output: none
#	toCall: &printCountTSVInput($countInfoHsh_ref, $meanCountCutoff, $minCountInOneSample, $TSVToPrintPath, $refSample, $qrySample);
#	calledInLine: 1349
#....................................................................................................................................................#

	my ($countInfoHsh_ref, $meanCountCutoff, $minCountInOneSample, $TSVToPrintPath, $refSample, $qrySample) = @_;
	
	#---filtering
	my %validIDHsh;
	foreach my $geneID (keys %{$countInfoHsh_ref}) {
		my $countSum = 0;
		my $sampleNum = 0;
		my $moreThanMinCount = 'no';
		foreach my $sampleName (keys %{$countInfoHsh_ref->{$geneID}}) {
			foreach my $repNum (keys %{$countInfoHsh_ref->{$geneID}{$sampleName}}) {
				$sampleNum++;
				$countSum += $countInfoHsh_ref->{$geneID}{$sampleName}{$repNum};
				$moreThanMinCount = 'yes' if $countInfoHsh_ref->{$geneID}{$sampleName}{$repNum} >= $minCountInOneSample;
			}
		}
		my $meanCount = $countSum/$sampleNum;
		if (($meanCount >= $meanCountCutoff) and ($moreThanMinCount eq 'yes')) {
			$validIDHsh{$geneID}++ 
		}
	}

	my $totalGeneNum = keys %{$countInfoHsh_ref};
	my $validGeneNum = keys %validIDHsh;
	
	#print "$validGeneNum out of $totalGeneNum passed the filter of meanCountCutoff $meanCountCutoff minCountInOneSample $minCountInOneSample\n";

	open (DESEQTSV, ">$TSVToPrintPath");
	my @headerAry;
	foreach my $geneID (sort {$a cmp $b} keys %{$countInfoHsh_ref}) {
		push @headerAry, 'geneID';
		foreach my $sampleName (($refSample, $qrySample)) {
			foreach my $repNum (sort {$a <=> $b} keys %{$countInfoHsh_ref->{$geneID}{$sampleName}}) {
				push @headerAry, $sampleName."_".$repNum;
			}
		}
		last;
	}
	print DESEQTSV join "", ((join "\t", @headerAry),"\n");
	
	foreach my $geneID (sort {$a cmp $b} keys %{$countInfoHsh_ref}) {
		next unless exists $validIDHsh{$geneID};
		my @outputAry;
		push @outputAry, $geneID;
		foreach my $sampleName (($refSample, $qrySample)) {
			foreach my $repNum (sort {$a <=> $b} keys %{$countInfoHsh_ref->{$geneID}{$sampleName}}) {
				push @outputAry, int ($countInfoHsh_ref->{$geneID}{$sampleName}{$repNum});
			}
		}
		print DESEQTSV join "", ((join "\t", @outputAry),"\n");
	}
	close DESEQTSV;
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|81
#	secondaryAppearInSection: >none
#	input: none
#	output: $countType, $geneNameLenPath, $outDir, $refNameAry_ref, $sampleInfoPath
#	toCall: my ($sampleInfoPath, $refNameAry_ref, $geneNameLenPath, $countType, $outDir) = &readParameters();
#	calledInLine: 88
#....................................................................................................................................................#
	
	my ($sampleInfoPath, $refNameAry_ref, $geneNameLenPath, $countType, $outDir);

	$outDir = "default"; #---same location as the sample info dir
	$countType = 'express_eff_count';
	
	GetOptions(	
		'sampleInfoPath=s' => \$sampleInfoPath, 
		'refName:s@' => \$refNameAry_ref,
		'geneNameLenPath=s' => \$geneNameLenPath, 
		'countType:s' => \$countType, 
		'outDir:s' => \$outDir,
	);
	
	my ($sampleInfoName,$sampleInfoDir,$sampleInfoSuffix) = fileparse($sampleInfoPath, qr/\.[^.]*/);
	
	chomp (my $pwd = `pwd`);
	$sampleInfoDir = $pwd if $sampleInfoDir eq "./";
	$outDir = $sampleInfoDir."/".$sampleInfoName if ($outDir eq "default");
	
	system ("mkdir -p -m 777 $outDir");
	
	return ($sampleInfoPath, $refNameAry_ref, $geneNameLenPath, $countType, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|355
#	appearInSub: getAllCountInfo|1056, getGeneNameAndLength|1183
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 2_processInputData|108
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 1075, 1213
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->355

	return ();
}
sub runAllComparisons {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: generateComparisonCombination|469, generateDESeq1RScript|516, generateDESeq2RScript|698, generateEdgeRRScript|776, launchAndMonitorAllCMD|1219, overlapGeneExpDiffResults|1398, printAllComparisonFoldChange|1770, printCombineResults|1951
#	appearInSub: >none
#	primaryAppearInSection: 4_runDEAnalyses|131
#	secondaryAppearInSection: >none
#	input: $BCV, $DESeq1FDR, $DESeq2FDR, $EdgeRFDR, $FDRTag, $allSampleCountInfoHsh_ref, $countType, $geneNameLenHsh_ref, $outDir, $refNameHsh_ref, $sampleInfoHsh_ref, $sampleInfoPath
#	output: $allComparisonFoldChangeHsh_ref, $allComparisonNameHsh_ref, $allResultFilePathHsh_ref, $overallCountFilePathHsh_ref
#	toCall: my ($allComparisonFoldChangeHsh_ref, $allResultFilePathHsh_ref, $allComparisonNameHsh_ref, $overallCountFilePathHsh_ref) = &runAllComparisons($sampleInfoPath, $outDir, $allSampleCountInfoHsh_ref, $DESeq1FDR, $DESeq2FDR, $EdgeRFDR, $FDRTag, $refNameHsh_ref, $geneNameLenHsh_ref, $sampleInfoHsh_ref, $countType, $BCV);
#	calledInLine: 136
#....................................................................................................................................................#
	
	my ($sampleInfoPath, $outDir, $allSampleCountInfoHsh_ref, $DESeq1FDR, $DESeq2FDR, $EdgeRFDR, $FDRTag, $refNameHsh_ref, $geneNameLenHsh_ref, $sampleInfoHsh_ref, $countType, $BCV) = @_;

	my %allComparisonFoldChangeHsh = ();
	my $allComparisonFoldChangeHsh_ref = \%allComparisonFoldChangeHsh;
	my $allComparisonNameHsh_ref = {};
	my $allResultFilePathHsh_ref = {};
	my $overallCountFilePathHsh_ref = {};
	my ($sampleInfoName,$sampleInfoDir,$sampleInfoSuffix) = fileparse($sampleInfoPath, qr/\.[^.]*/);
	my $allResultList_FH;
	open $allResultList_FH, ">", "$outDir/$sampleInfoName.allDEGList.txt";

	#---repeat exactly the same thing for all refNames
	foreach my $refName (sort {$a cmp $b} keys %{$refNameHsh_ref}) {

		system ('clear');
		print "Preparing to staring processing reference sample $refName\n";
	
		#---&generateComparisonCombination
		my ($comparisonInfoHsh_ref) = &generateComparisonCombination($sampleInfoHsh_ref, $refName, $allComparisonNameHsh_ref, $outDir);#->469

		#---generate the DESeq1 RScripts
		my $DESeq1FitType = 'parametric';
		my $DESeq1Method = 'pooled';
		my $DESeq1SharingMode = 'maximum';
		&generateDESeq1RScript($comparisonInfoHsh_ref, $DESeq1FDR, $DESeq1FitType, $DESeq1SharingMode, $DESeq1Method, $outDir);#->516

		#---generate the DESeq1 RScripts
		my $DESeq2FitType = 'parametric';
		&generateDESeq2RScript($comparisonInfoHsh_ref, $DESeq2FDR, $DESeq2FitType, $outDir);#->698
	
		#---generate the EdgeR RScripts
		&generateEdgeRRScript($comparisonInfoHsh_ref, $BCV, $outDir);#->776
	
		#---launch the cuffdiff and DESeq
		&launchAndMonitorAllCMD($comparisonInfoHsh_ref, $outDir, $allSampleCountInfoHsh_ref, $countType);#->1219

		my $combinedResultInfoHsh_ref = {};
		my $combinedResultCountHsh_ref = {};
		my $combinedCutoffOvrlpHsh_ref = {};

		foreach my $minRPKM ((0)) {
		#foreach my $minRPKM ((2)) {
			#foreach my $minLinearFC ((0, 2, 5, 10, 15)) {#----must include 0 to collect all data
			foreach my $minLinearFC ((0)) {#----must include 0 to collect all data
				#foreach my $upDnAll (('all', 'up', 'dn')) {#----must include 'all' to collect all data
				foreach my $upDnAll (('all')) {#----must include 'all' to collect all data
					&overlapGeneExpDiffResults ($comparisonInfoHsh_ref, $DESeq1FDR, $DESeq2FDR, $EdgeRFDR, $minLinearFC, $upDnAll, $combinedResultInfoHsh_ref, $combinedResultCountHsh_ref, $combinedCutoffOvrlpHsh_ref, $allComparisonFoldChangeHsh_ref, $minRPKM, $geneNameLenHsh_ref, $allSampleCountInfoHsh_ref, $allResultFilePathHsh_ref);#->1398
				}
			}
		}

		&printCombineResults($comparisonInfoHsh_ref, $combinedResultCountHsh_ref, $combinedResultInfoHsh_ref, $combinedCutoffOvrlpHsh_ref, $geneNameLenHsh_ref, "$outDir/$refName/combine", $FDRTag, $allResultList_FH, $allResultFilePathHsh_ref, $allSampleCountInfoHsh_ref, $refName, $overallCountFilePathHsh_ref);#->1951
		&printAllComparisonFoldChange($allComparisonFoldChangeHsh_ref, $geneNameLenHsh_ref, "$outDir/$refName/combine/ovrlp/", "$refName", $comparisonInfoHsh_ref, $FDRTag);#->1770

	}

	return ($allComparisonFoldChangeHsh_ref, $allResultFilePathHsh_ref, $allComparisonNameHsh_ref, $overallCountFilePathHsh_ref);
	
}
sub runAndCheckSerialTask {
#....................................................................................................................................................#
#	subroutineCategory: taskManage
#	dependOnSub: >none
#	appearInSub: overlapGeneExpDiffResults|1398
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $cmd, $errorLogPath, $grepCmd, $grepStr
#	output: none
#	toCall: &runAndCheckSerialTask($grepCmd, $grepStr, $cmd, $errorLogPath);
#	calledInLine: 1688
#....................................................................................................................................................#

	my ($grepCmd, $grepStr, $cmd, $errorLogPath) = @_;

	system (qq|$cmd 2>>$errorLogPath &|);
	my $sdout = $grepStr;
	while ($sdout =~ m/$grepStr/) {
		$sdout = `$grepCmd`;
		sleep (0.001);
	}
}
sub tripleVennDiagram {
#....................................................................................................................................................#
#	subroutineCategory: plotInR
#	dependOnSub: >none
#	appearInSub: overlapGeneExpDiffResults|1398
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $RScrLog, $colorA, $colorB, $colorC, $listAHsh_ref, $listAName, $listBHsh_ref, $listBName, $listCHsh_ref, $listCName, $pdfPathPrefix
#	output: $pdfPath, \%ABCMergeHsh, \%ABCOvrlpHsh, \%ABOvrlpHsh, \%ACOvrlpHsh, \%BCOvrlpHsh
#	toCall: my (\%ABOvrlpHsh, \%BCOvrlpHsh, \%ACOvrlpHsh, \%ABCOvrlpHsh, \%ABCMergeHsh, $pdfPath) = &tripleVennDiagram($listAHsh_ref, $listBHsh_ref, $listCHsh_ref, $listAName, $listBName, $listCName, $pdfPathPrefix, $RScrLog, $colorA, $colorB, $colorC);
#	calledInLine: 1621
#....................................................................................................................................................#

	my ($listAHsh_ref, $listBHsh_ref, $listCHsh_ref, $listAName, $listBName, $listCName, $pdfPathPrefix, $RScrLog, $colorA, $colorB, $colorC) = @_;
	
	my $listANum = keys %{$listAHsh_ref};
	my $listBNum = keys %{$listBHsh_ref};
	my $listCNum = keys %{$listCHsh_ref};
	
	my %ABOvrlpHsh;
	my %BCOvrlpHsh;
	my %ACOvrlpHsh;
	my %ABCOvrlpHsh;
	my %ABCMergeHsh;
	my %AOnlyHsh;
	my %BOnlyHsh;
	my %COnlyHsh;

	my $pdfPath = "";

	if (($listANum > 0) and ($listBNum > 0) and ($listCNum > 0)) {
		
		my %allIdHsh;
		my %pairOvrlpIDHsh;
		my %onlyIDHsh;
		foreach my $id (keys %{$listAHsh_ref}) {${$allIdHsh{'A'}}{$id}++;}
		foreach my $id (keys %{$listBHsh_ref}) {${$allIdHsh{'B'}}{$id}++;}
		foreach my $id (keys %{$listCHsh_ref}) {${$allIdHsh{'C'}}{$id}++;}

		my $listNum = keys %allIdHsh;
		
		foreach my $listRef (keys %allIdHsh) {
			
			%{$onlyIDHsh{$listRef}} = ();
			foreach my $listQry (keys %allIdHsh) {
				%{${$pairOvrlpIDHsh{$listRef}}{$listQry}} = ();
			}
			
			foreach my $id (keys %{$allIdHsh{$listRef}}) {
				$ABCMergeHsh{$id}++;
				my $inOther = 0;
				foreach my $listQry (keys %allIdHsh) {
					next if $listRef eq $listQry;
					if (exists ${$allIdHsh{$listQry}}{$id}) {
						$inOther++;
						${${$pairOvrlpIDHsh{$listRef}}{$listQry}}{$id}++;
					}
				}
				
				$ABCOvrlpHsh{$id}++ if ($inOther == $listNum-1);
				${$onlyIDHsh{$listRef}}{$id}++ if $inOther == 0;
			}
		}
		
		%AOnlyHsh = %{$onlyIDHsh{'A'}};
		%BOnlyHsh = %{$onlyIDHsh{'B'}};
		%COnlyHsh = %{$onlyIDHsh{'C'}};
		%ABOvrlpHsh = %{${$pairOvrlpIDHsh{'A'}}{'B'}};
		%BCOvrlpHsh = %{${$pairOvrlpIDHsh{'B'}}{'C'}};
		%ACOvrlpHsh = %{${$pairOvrlpIDHsh{'A'}}{'C'}};

		my $AOnlyNum = keys %AOnlyHsh;
		my $BOnlyNum = keys %BOnlyHsh;
		my $COnlyNum = keys %COnlyHsh;

		my $ABOvrlpNum = keys %ABOvrlpHsh;
		my $BCOvrlpNum = keys %BCOvrlpHsh;
		my $ACOvrlpNum = keys %ACOvrlpHsh;
		my $ABCOvrlpNum = keys %ABCOvrlpHsh;
		my $ABCMergeNum = keys %ABCMergeHsh;
		
		my %nameColorHsh;
		
		${$nameColorHsh{'A'}}{"only"} = $AOnlyNum;
		${$nameColorHsh{'B'}}{"only"} = $BOnlyNum;
		${$nameColorHsh{'C'}}{"only"} = $COnlyNum;
		${$nameColorHsh{'A'}}{"num"} = $listANum;
		${$nameColorHsh{'B'}}{"num"} = $listBNum;
		${$nameColorHsh{'C'}}{"num"} = $listCNum;
		${$nameColorHsh{'A'}}{"name"} = $listAName."[$listANum]";
		${$nameColorHsh{'B'}}{"name"} = $listBName."[$listBNum]";
		${$nameColorHsh{'C'}}{"name"} = $listCName."[$listCNum]";
		${$nameColorHsh{'A'}}{"color"} = $colorA;
		${$nameColorHsh{'B'}}{"color"} = $colorB;
		${$nameColorHsh{'C'}}{"color"} = $colorC;

		my @colorAry;
		my @nameAry;
		
		foreach my $list (sort {${$nameColorHsh{$b}}{"num"} <=> ${$nameColorHsh{$a}}{"num"}} keys %nameColorHsh) {
			push @colorAry, ${$nameColorHsh{$list}}{"color"};
			push @nameAry, ${$nameColorHsh{$list}}{"name"};
		}

		$pdfPath = "$pdfPathPrefix.ovrlp.$ABCOvrlpNum.merge.$ABCOvrlpNum.pdf";
		my $vennplotCmd = "venn.plot <- draw.triple.venn( area1=$listANum, area2=$listBNum, area3=$listCNum, n12=$ABOvrlpNum, n23=$BCOvrlpNum, n13=$ACOvrlpNum, n123=$ABCOvrlpNum, category = c(\"$nameAry[0]\", \"$nameAry[2]\", \"$nameAry[1]\"), fill = c(\"$colorAry[0]\", \"$colorAry[2]\", \"$colorAry[1]\"), lty = \"blank\", cex = 2, cat.cex = 2, cat.pos = c(0, 90, 180), cat.dist = 0.05, cat.just = list(c(0.5, 2), c(1, -9.5), c(0.5, 0)), cat.col = c(\"$colorAry[0]\", \"$colorAry[2]\", \"$colorAry[1]\"), ext.pos = 30, ext.dist = -0.05, ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = \"dashed\", scaled = FALSE, euler.d = TRUE )";
		#my $vennplotCmd = "venn.plot <- draw.triple.venn( area1=$listANum, area2=$listBNum, area3=$listCNum, n12=$ABOvrlpNum, n23=$BCOvrlpNum, n13=$ACOvrlpNum, n123=$ABCOvrlpNum, category = c(\"$nameAry[0]\", \"$nameAry[1]\", \"$nameAry[2]\"), fill = c(\"$colorA\", \"$colorB\", \"$colorC\"), lty = \"blank\", cex = 2, cat.cex = 2, cat.pos = c(0, 90, 180), cat.dist = 0.05, cat.just = list(c(0.5, 2), c(1, -9.5), c(0.5, 0)), cat.col = c(\"$colorA\", \"$colorB\", \"$colorC\"), ext.pos = 30, ext.dist = -0.05, ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = \"dashed\")";
		
		open (RSCRLOG, ">$RScrLog");
		print RSCRLOG "########## R CMD START #############\n";
		print RSCRLOG "library(VennDiagram)\n";
		print RSCRLOG "$vennplotCmd\n";
		print RSCRLOG "pdf(\"$pdfPath\")\n";
		print RSCRLOG "grid.draw(venn.plot)\n";
		print RSCRLOG "dev.off()\n";
		print RSCRLOG "########## R CMD END #############\n";
		close RSCRLOG;

		open R, "|/usr/bin/R --vanilla --slave 1>>$RScrLog 2>>$RScrLog";
		print R "library(VennDiagram)\n";
		print R "$vennplotCmd\n";
		print R "pdf(\"$pdfPath\")\n";
		print R "grid.draw(venn.plot)\n";
		print R "dev.off()\n";
		close R;
	}

	return (\%ABOvrlpHsh, \%BCOvrlpHsh, \%ACOvrlpHsh, \%ABCOvrlpHsh, \%ABCMergeHsh, $pdfPath);
}
sub userDecisionToProceed {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: checkSampleFiles|294
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 2_processInputData|108
#	input: none
#	output: none
#	toCall: &userDecisionToProceed();
#	calledInLine: 350
#....................................................................................................................................................#
	my $userDecision = "";
	while (($userDecision ne "y") and ($userDecision ne "n") ) {
		print "\nReturn <y> to proceed or <n> to quit (y/n):";
		chomp ($userDecision = <STDIN>);
		last if ($userDecision eq "y");
		die "User choose to terminate the program. Quitting\n" if ($userDecision eq "n");
	}
}

exit;
