#! perl -w
# Put your pair-end fastq files (named as XX_R1.fastq(.gz) and XX_R2.fastq(.gz)), barcode fastq file (named as XX_I1.fastq(.gz), digits like "_001" can be added after R1, R2 and I1) and map file (named as XX_map.txt) in single folders (named as dataset_1, dataset_2, dataset_3 etc.). 
# For a better view of PCA, you can group the samples by inserting a collumn of Group before the collumn Description in the map file.
# Set the parameters below and run this script.
#

######### Parameter ##############$classification="16srrna"; # 18srrna, 16srrna, fungallsu, fungalits_warcup, fungalits_unite (note their release date at: https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData).
$classification="16srrna"; # 18srrna, 16srrna, fungallsu, fungalits_warcup, fungalits_unite (note their release date at: https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData).
$rdp_date{"fungalits_warcup"}="July, 2016";
$rdp_date{"fungalits_unite"}="July, 2014";
$rdp_date{"fungallsu"}="Mar, 2014";
$rdp_date{"16srrna"}="July, 2016";


$bar_extract=1; # Need to extract barcode or not before demultiplex? Yes: 1; No: 0. If you only have two fastq files, yes you need to do it.
$bar_in_forw=1; # Is barcode linked in forward primer? Yes: 1; Reverse primer: 0, if both, make it 1. If you already have barcode fastq file (e.g. XX_I1.fastq), ignore this parameter. 
$map_rc=0; # Is barcode in fastq and barcode in mapping file reverse completentary? Yes: 1; No: 0. If you only have two fastq files, make it 0.
$quality_cutoff="20"; # Quality cutoff for quality filtering and trimming in sequence merging and demultiplexing.
$baseN=0; # Minimum N allowed in demultiplexing
$bar_mis=0; # Maximum mismatch in barcode in demultiplexing
$min_length=150; # Minimum length in pair-end merging including primers
$max_length=1500; # Maximum length in pair-end merging including primers
$overlap="50"; # The minimum overlap of pair-end reads when merging reads. Should be >30.
$pvalue=0.0001; # The p value cutoff of overlap when merging reads. Should be <0.001
$depth="10"; # Resampling depth for OTU table. Make it 'auto' if you want as many OTUs and sample as possible.
$cpu="50"; # Number of threads you want to use.
$qual_check=0; # Do you want to check sequencing quality check using FastQC? Yes: 1; No: 0.


######## Configuration ##########
$home="/home/tianrm/software/scripts/pipelines/ASAP/ASAP_1.4";
$rdp_classifier="/home/tianrm/software/rdp_classifier_2.12/dist/classifier.jar"; # Install RDP Classifier with the latest trainset: https://sourceforge.net/projects/rdp-classifier.
$fastqc="/home/tianrm/software/FastQC/fastqc";
$mafft="/home/tianrm/anaconda2/bin/mafft";
$pr_path="/raid_data2/tianrm/database/pr2/pr2_version_4.7_UTAX_Eukaryota.udb"; # If analyzing 18S data
$map="map_merged.txt"; # leave it alone

system "set -e; print_qiime_config.py; touch ready";
if(-f "ready"){
	system "rm -rf ready";
}else{
	die "Hey dude. Cannot run print_qiime_config.py. I bet you forgot to activate the environment.\n";
}

open LOG,">log.txt";

# Search dataset
opendir D,".";
@dir=readdir D;
foreach(@dir){
	if(/dataset_\d+/){
		push @dataset,$_;
	}
}
$nd=@dataset+0;
if($nd<1){
	die "Hey dude, dataset files not found!\n";
}

system "rm -rf split; mkdir split; rm -rf map_merged.txt; touch map_merged.txt";

# Extract barcode sequence
if($bar_extract==1){
	# Find fastq and map
	foreach(@dataset){
		$dataset=$_;
		chdir $dataset;
		opendir D,".";
		@dir=readdir D;
	        foreach(@dir){
	                if(/R1_?\d*\.fa?s?t?q$/){
	                        $raw_fq_r1=$_;
	                }elsif(/R2_?\d*\.fa?s?t?q$/){
	                        $raw_fq_r2=$_;
	                }elsif(/R1_?\d*\.fa?s?t?q\.gz$/){
	                        $raw_fq_r1=$_;
	                        system "gzip -d $raw_fq_r1";
	                        $raw_fq_r1=~s/\.gz$//;
	                }elsif(/R2_?\d*\.fa?s?t?q\.gz$/){
	                        $raw_fq_r2=$_;
	                        system "gzip -d $raw_fq_r2";
	                        $raw_fq_r2=~s/\.gz$//;
	                }elsif(/.*?map.*?\.txt$/){
	                        $map_file=$_;
	                }
	        }

	        unless($raw_fq_r1 & $raw_fq_r2 & $map_file){
	                die "\n\nHey dude, $dataset not detected properly. Please make sure all the required files are correctly named.\n";
	        }else{
	                print LOG "\n\nFiles of $dataset are all present.\n";
	        }
	
	
	        # Check the map file
	        print LOG "Checking map file of $dataset ...\n";
	        open IN,"$map_file";
	        $line=<IN>;
	        $line=<IN>;
	        @line=split/\t/,$line;
	        $barcode=$line[1];
	        $barcode_type=length$barcode;
	        system "rm -rf map_validation; validate_mapping_file.py -m $map_file -o map_validation > map_validation.txt";
	        open F,"map_validation.txt";
	        $map_vali=<F>;
	        if($map_vali=~/No errors or warnings were found/){
	                print LOG "\nMap file is good.\n\n";
	        }else{
	                die "\nMap file is problematic. Please see the directory map_validation";
	        }
	        system "mv map_validation.txt map_validation";

		# Extract barcode
		if($bar_in_forw==1){
			system "set -ex; extract_barcodes.py -f $raw_fq_r1 -c barcode_single_end --bc1_len $barcode_type -o extract_barcode; touch extract_barcode_done";
		        unless(-e "extract_barcode_done"){
	                die "Extraction of barcode failed.\n";
		        }else{
        		        system "rm -rf extract_barcode_done";
		        }
			system "mv extract_barcode/barcodes.fastq barcodes_I1.fastq; rm -rf extract_barcode";
		}else{
                        system "set -ex; extract_barcodes.py -f $raw_fq_r2 -c barcode_single_end --bc1_len $barcode_type -o extract_barcode; touch extract_barcode_done";
                        unless(-e "extract_barcode_done"){
                        die "Extraction of barcode failed.\n";
                        }else{
                                system "rm -rf extract_barcode_done";
                        }
			system "mv extract_barcode/barcodes.fastq barcodes_I1.fastq; rm -rf extract_barcode";
		}


		chdir "..";
	}
}


# Preprocess for each dataset
foreach(@dataset){
	$dataset=$_;
	chdir $dataset;
	opendir D,".";
	@dir=readdir D;
	foreach(@dir){
		if(/R1_?\d*\.fa?s?t?q$/){
			$raw_fq_r1=$_;
		}elsif(/R2_?\d*\.fa?s?t?q$/){
			$raw_fq_r2=$_;
		}elsif(/I1_?\d*\.fa?s?t?q$/){
			$raw_fq_barcode=$_;
                }elsif(/R1_?\d*\.fa?s?t?q\.gz$/){
                        $raw_fq_r1=$_;
			system "gzip -d $raw_fq_r1";
			$raw_fq_r1=~s/\.gz$//;
                }elsif(/R2_?\d*\.fa?s?t?q\.gz$/){
                        $raw_fq_r2=$_;
			system "gzip -d $raw_fq_r2";
			$raw_fq_r2=~s/\.gz$//;
                }elsif(/I1_?\d*\.fa?s?t?q\.gz$/){
                        $raw_fq_barcode=$_;
			system "gzip -d $raw_fq_barcode";
			$raw_fq_barcode=~s/\.gz$//;
		}elsif(/.*?map.*?\.txt$/){
			$map_file=$_;
		}
	}

	unless($raw_fq_r1 & $raw_fq_r2 & $raw_fq_barcode & $map_file){
		die "\n\nHey dude, $dataset not detected properly. Please make sure all the required files are correctly named.\n";
	}else{
		print LOG "\n\nFiles of $dataset are all present.\n";
	}


	# Check the map file
	print LOG "Checking map file of $dataset ...\n";
	open IN,"$map_file";
	$line=<IN>;
	$line=<IN>;
	@line=split/\t/,$line;
	$barcode=$line[1];
	$barcode_type=length$barcode;
	system "rm -rf map_validation; validate_mapping_file.py -m $map_file -o map_validation > map_validation.txt";
	open F,"map_validation.txt";
	$map_vali=<F>;
	if($map_vali=~/No errors or warnings were found/){
		print LOG "\nMap file is good.\n\n";
	}else{
		die "\nMap file is problematic. Please see the directory map_validation";
	}
	system "mv map_validation.txt map_validation";


	# Read quality check
	if($qual_check==1){
		print LOG "Checking read quality of $dataset ...\n";
		system "rm -rf quality_check; mkdir quality_check";
		system "set -ex; nohup $fastqc -f fastq -o quality_check -t $cpu $raw_fq_r1 > quality_check/log; touch fastqc1_done";
		system "set -ex; nohup $fastqc -f fastq -o quality_check -t $cpu $raw_fq_r2 >> quality_check/log; touch fastqc2_done";
		unless(-e "fastqc1_done" && -e "fastqc2_done"){
			die "Quality check failed.\n";
		}else{
			system "rm -rf fastqc1_done fastqc2_done";
		}
		print LOG "Quality check done. See quality_check for raw read number and quality.\n\n";
	}

	# Merge pair-end sequence files
	print LOG "Merging pair-end sequence files of $dataset ...\n";
	system "rm -rf read_merge; mkdir read_merge";
	system "set -ex; $home/bin/pear-0.9.10-bin-64 -f $raw_fq_r1 -r $raw_fq_r2 -n $min_length -m $max_length -p $pvalue -o read_merge/merged -q $quality_cutoff -j $cpu -v $overlap > read_merge/pear.log; touch merge_done";
	unless(-e "merge_done"){
		die "Merge failed.\n";
	}else{
		system "rm -rf merge_done";
	}
	system "gzip $raw_fq_r1";
	system "gzip $raw_fq_r2";
	print LOG "Merge done.\n\n";



	# Because some reads were not merged, their corresponding barcode need to be removed.
	print LOG "Extracting barcode fq of $dataset ...\n";
	system "set -ex; perl $home/bin/extract_corresponding_barcode_fq.pl $raw_fq_barcode read_merge/merged.assembled.fastq; touch extract_barcode_done";
	unless(-e "extract_barcode_done"){
	        die "Extraction of barcode fq failed.\n";
	}else{
	        system "rm -rf extract_barcode_done";
	}
	if($bar_extract==1){
		system "rm -rf $raw_fq_barcode";
	}else{
		system "gzip $raw_fq_barcode";
	}
	print LOG "Extraction of barcode fq done.\n\n";


	# Demultiplexing
	print LOG "Demultiplexing $dataset ...\n";
	$barcode_fq=$raw_fq_barcode;
	$barcode_fq=~s/\.fq$|\.fastq$//;
	$barcode_fq.="_extracted_barcode.fastq";
	system "rm -rf split";
	if($map_rc==1){
		system "set -ex; split_libraries_fastq.py -i read_merge/merged.assembled.fastq -n $baseN --max_barcode_errors $bar_mis -b $barcode_fq -m $map_file --rev_comp_mapping_barcodes --barcode_type $barcode_type -q $quality_cutoff -o split; touch split_done";
	}elsif($map_rc==0){
		system "set -ex; split_libraries_fastq.py -i read_merge/merged.assembled.fastq -n $baseN --max_barcode_errors $bar_mis -b $barcode_fq -m $map_file --barcode_type $barcode_type -q $quality_cutoff -o split; touch split_done";
	}
	
	unless(-e "split_done"){
	        die "Demultiplexing failed.\n";
	}else{
	        system "rm -rf split_done";
	}
	system "rm -rf $barcode_fq";
	print LOG "Demultiplexing done.\n\n";




	# Trim forward and reverse primer and any extra bases
	print LOG "Trimming forward and reverse primer and any extra bases from the seqs.fna of $dataset ...\n";
	system "set -ex; perl $home/bin/trim_primer.pl $map_file split/seqs.fna split/seqs_trim_primer.fna; touch trim_primer_done";
	unless(-e "trim_primer_done"){
	        die "Trim primers failed.\n";
	}else{
	        system "rm -rf trim_primer_done";
	}
	system "grep -c -e \\> split/seqs.fna > num";
	system "grep -c -e \\> split/seqs_trim_primer.fna >> num";
	open F,"num";
	@num=<F>;
	chomp $num[0];
	chomp $num[1];
	print LOG "Before primer trimming: $num[0] sequences\nAfter primer trimming: $num[1] sequences\n";

        if($num[0]<1000 || $num[1]<1000){
                die "No reads or too few reads extracted for $dataset.\n";
        }

	$ratio=$num[0]/$num[1];
	if($ratio>3){
		die "Sequence number after primer trimming is too small. Please check it.\n";
	}
	system "rm -rf num";
	print LOG "Primer trimming done.\n\n";



	# Merge sequence of different dataset
	print LOG "Merging $dataset to total dataset.\n";
	system "cat split/seqs_trim_primer.fna >> ../split/seqs_trim_primer.fna";
	system "rm -rf split/seqs*";

	# Merge map file
	print LOG "Merging map file of $dataset.\n";
	open F, $map_file;
	open FI,">>../map_merged.txt";
	while(<F>){
		@line=split/\t/,$_;
		unless($map{$line[0]}){
			print FI $_;
		}
		$map{$line[0]}++;
	}

        $raw_fq_r1=();
        $raw_fq_r2=();
        $raw_fq_barcode=();
        $map_file=();

	chdir "..";

}


# Extract fastq for data submission use
print LOG "Extracting fastq for data submission use...\n";
system "perl $home/bin/simplify_head.pl split/seqs_trim_primer.fna > split/seqs_trim_primer_simplified.fna";
system "mv split/seqs_trim_primer_simplified.fna split/seqs_trim_primer.fna";

system "set -ex; perl $home/bin/extract_fastq.pl; touch extract_done";
unless(-e "extract_done"){
        die "Extraction failed.\n";
}else{
        system "rm -rf extract_done";
}
system "gzip split/seqs_trim_primer.fastq";
print LOG "Extraction done.\n\n";
system "rm -rf dataset*/read_merge";


# Dereplicate reads and calculate abundance
print LOG "\n\nDereplicating reads and calculating abundance...\n";
system "grep -c -e \\> split/seqs_trim_primer.fna > num";
open F,"num";
$num=<F>;
chomp $num;
print LOG "In total $num sequences assigned to samples.\n";
system "rm -rf num";


# usearch 9
# system "set -ex; $home/bin/usearch_9.2.64 -fastx_uniques split/seqs_trim_primer.fna -fastaout split/dereplicated.fasta -sizeout > split/dereplicate.log; touch dereplication_done";

# usearch 7 64 bit
#system "set -ex; $home/bin/usearch_7.0_64 -derep_fulllength split/seqs_trim_primer.fna -output split/dereplicated.fasta -sizeout -minseqlength 50 -threads $cpu > split/dereplicate.log; touch dereplication_done";

# vsearch
system "set -ex; $home/bin/vsearch --derep_fulllength split/seqs_trim_primer.fna --output split/dereplicated.fasta --sizeout --threads $cpu > split/dereplicate.log; touch dereplication_done";


unless(-e "dereplication_done"){
        die "Dereplication failed.\n";
}else{
        system "rm -rf dereplication_done";
}
print LOG "Dereplication done.\n\n";


# OTU clustering, chimera and singleton removal
print LOG "Removing singleton..\n";
system "set -ex; perl $home/bin/remove_singleton.pl split/dereplicated.fasta split/dereplicated_1.fasta; touch singleton_done";
unless(-e "singleton_done"){
        die "Revomving singleton failed.\n";
}else{
        system "mv split/dereplicated_1.fasta split/dereplicated.fasta; rm -rf singleton_done";
}
print LOG "Removing singleton done.\n\n";


print LOG "UPARSE: OTU clustering, chimera and singleton removal...\n";
close LOG;
system "rm -rf OTU; mkdir OTU; mkdir OTU/log";
system "set -ex; $home/bin/usearch_10.0.240 -cluster_otus split/dereplicated.fasta -otus OTU/otus.fasta -uparseout OTU/log/out.up -relabel OTU -minsize 2 &>> log.txt; touch otu_done";

unless(-e "otu_done"){
        die "OTU clustering failed.\n";
}else{
        system "rm -rf otu_done";
}

open LOG,">>log.txt";
print LOG "OTU clustering done.\n\n";


# Modify identifier of seqs_trim_primer.fna
system "bash $home/bin/modify_identifier.sh split/seqs_trim_primer.fna";


# Map reads to OTU sequences and make OTU table
print LOG "Map reads to OTU sequences and make OTU table...\n";
#system "set -ex; $home/bin/usearch_9.2.64 -usearch_global split/seqs_trim_primer.fna -db OTU/otus.fasta -strand plus -id 0.97 -otutabout OTU/otu_table_raw.txt > OTU/log/otu_table.log; touch otu_table_done";

system "set -ex; $home/bin/vsearch --usearch_global split/seqs_trim_primer.fna --db OTU/otus.fasta --strand plus -id 0.97 --otutabout OTU/otu_table_raw.txt --threads $cpu > OTU/log/otu_table.log; touch otu_table_done";

unless(-e "otu_table_done"){
        die "OTU table failed.\n";
}else{
        system "rm -rf otu_table_done";
}
system "wc -l OTU/otu_table_raw.txt > num";
open F,"num";
$num=<F>;
$num=~s/.*?(\d+).*/$1/;
$num-=1;
print LOG "In total $num OTUs acquired.\n";
system "rm -rf num";


system "rm -rf split/dereplicated.fasta";
print LOG "OTU table done.\n\n";
system "rm -rf split/seqs_trim_primer.fna";


chdir "OTU";


# Taxonomic classification of OTU
if($classification eq "18srrna"){
	print LOG "Classification using PR database...\n";
	system "rm -rf sintax_assigned_taxonomy; mkdir sintax_assigned_taxonomy";
	system "set -ex; $home/bin/usearch_10.0.240 -sintax otus.fasta -db $pr_path -tabbedout sintax_assigned_taxonomy/tax_assignments_raw.txt -strand both -sintax_cutoff 0.8; touch sintax_done";
	unless(-f "sintax_done"){
		die "SINTAX Classification failed.\n";
	}else{
		system "rm -rf sintax_done";
	}
	open F,"sintax_assigned_taxonomy/tax_assignments_raw.txt";
	@file=<F>;
	$file=join'',@file;
	$file=~s/^(.*?)\t.*?\t[+-]\t(.*?)\n/$1\t$2\n/gm;
	$file=~s/k:/D_0__/g;
	$file=~s/d:/D_1__/g;
	$file=~s/p:/D_2__/g;
	$file=~s/c:/D_3__/g;
	$file=~s/o:/D_4__/g;
	$file=~s/f:/D_5__/g;
	$file=~s/g:/D_6__/g;
	$file=~s/s:/D_7__/g;
	open F,">sintax_assigned_taxonomy/tax_assignments_modified.txt";
	print F $file;
	print LOG "SINTAX Classification done.\n\n";

	# Add taxonomy to OTU table
        print LOG "Adding taxonomy to OTU table...\n\n";
        system "perl $home/bin/add_taxonomy_to_otu_table.pl sintax_assigned_taxonomy/tax_assignments_modified.txt otu_table_raw.txt otu_table_with_tax.txt";

}elsif($classification eq "16srrna" || $classification eq "fungallsu" || $classification eq "fungalits_warcup" || $classification eq "fungalits_unite"){
	print LOG "Classifying OTU using RDP Classifier with trainset of $rdp_date{$classification}...\n";
	system "rm -rf rdp_assigned_taxonomy; mkdir rdp_assigned_taxonomy;";
	system "set -ex; java -Xmx1g -jar $rdp_classifier classify -g $classification -f filterbyconf -o rdp_assigned_taxonomy/tax_assignments_raw.txt otus.fasta > rdp_assigned_taxonomy/rdp_cla.log; touch rdp_cla_done";
	unless(-e "rdp_cla_done"){
	        die "RDP Classifier failed.\n";
	}else{
	        system "rm -rf rdp_cla_done";
	}
	system "perl $home/bin/rdp_clssifier_format_transform.pl rdp_assigned_taxonomy/tax_assignments_raw.txt > rdp_assigned_taxonomy/tax_assignments_modified.txt";
	print LOG "RDP Classifier done.\n\n";

	# Add taxonomy to OTU table
	print LOG "Adding taxonomy to OTU table...\n\n";
	system "perl $home/bin/add_taxonomy_to_otu_table.pl rdp_assigned_taxonomy/tax_assignments_modified.txt otu_table_raw.txt otu_table_with_tax.txt";
}


# Remove Chloroplast and Mitochondria from OTU table
print LOG "Removing Chloroplast and Mitochondria from OTU table...\n";
system "perl $home/bin/remove_chloroplast_from_otu_table.pl otu_table_with_tax.txt otu_table_with_tax_prok.txt";


# Convert and resample OTU table
# convert
print LOG "Convert and resample OTU table...\n";
system "set -ex; biom convert -i otu_table_with_tax_prok.txt -o otu_table_with_tax_prok.biom --to-hdf5 --table-type=\"OTU table\" --process-obs-metadata taxonomy; touch convert_done";
unless(-e "convert_done"){
        die "Convert failed.\n";
}else{
        system "rm -rf convert_done";
}
print LOG "Convert done.\n";
# read summary
system "set -ex; biom summarize-table -i otu_table_with_tax_prok.biom > read_summary.txt; touch summarize_otu_done";
unless(-e "summarize_otu_done"){
        die "Summarize otu failed.\n";
}else{
        system "rm -rf summarize_otu_done";
}
print LOG "Summarize OTU table done. See the read_summary.txt for read statistics.\n";

# resample

if($depth eq "auto"){
	open F,"read_summary.txt";
	@file=<F>;
	$file=join'',@file;
	$file=~s/,//g;
	
	if($file=~/Total count:.*?Mean: (\d+)/sgim){
		print LOG $&;
		print LOG "\n";
	}
	
	if($file=~/Median: (\d+)/){
		$median=$1;
	}
	if($file=~/Min: (\d+)/){
	        $min=$1;
	}
	
	$depth=$median/3;
	$depth=sprintf "%.0f",$depth;
	if($depth<$min){
		$depth=$min;
	}
	$min1=$min*1.5;
	if($depth>$min && $depth<$min1){
		$depth=$min;
	}
	if($depth<8000){
		$depth=8000;
	}
	$depth=~s/\d$/0/;
}

print LOG "Resampling depth will be $depth.\n";
system "set -ex; single_rarefaction.py -i otu_table_with_tax_prok.biom -o otu_table_resampled.biom -d $depth; touch resample_otu_done";
unless(-e "resample_otu_done"){
        die "Resampling OTU table failed.\n";
}else{
        system "rm -rf resample_otu_done";
}
system "set -ex; biom convert -i otu_table_resampled.biom -o otu_table_resampled.txt --to-tsv --header-key taxonomy; touch convert_done";
unless(-e "convert_done"){
        die "Converting resampled OTU table to text format failed.\n";
}else{
        system "rm -rf convert_done";
}
system "wc -l otu_table_resampled.txt > num";
open F,"num";
$num=<F>;
$num=~s/.*?(\d+).*/$1/;
$num-=1;
print LOG "In total $num OTUs after resampling.\n";
system "rm -rf num";

print LOG "Resampling OTU table done.\n\n";


# update the otus.fasta
open FI,">otus_resampled.fasta";
open F,"otu_table_resampled.txt";
@otuid=();
$head=<F>;
$head=();
while(<F>){
	@line=split/\t/,$_;
	push @otuid,$line[0];
}
open F,"otus.fasta";
%seq=();
while(<F>){
	if(/>(.*?)[\n ]/){
		$otu=$1;
		$seq{$otu}.=$_;
	}else{
		$seq{$otu}.=$_;
	}
}
foreach(@otuid){
	print FI $seq{$_};
}
print LOG "otus_resampled.fasta also updated.\n\n";


# Alignment and filtering
print LOG "Aligning OTUs...\n";
system "rm -rf mafft_aligned_seqs; mkdir mafft_aligned_seqs";
system "set -ex; $mafft --thread $cpu --quiet otus_resampled.fasta > mafft_aligned_seqs/otus_resampled_aligned.fasta; touch align_done";
unless(-e "align_done"){
        die "Alignment failed.\n";
}else{
        system "rm -rf align_done";
}
system "$home/bin/Gblocks mafft_aligned_seqs/otus_resampled_aligned.fasta -t=d -b4=3 -b5=h > mafft_aligned_seqs/gblock.log; touch filter_alignment_done";
system "rm -rf mafft_aligned_seqs/otus_resampled_aligned.fasta-gb.htm";
system "perl -pi -e 's/ //g' mafft_aligned_seqs/otus_resampled_aligned.fasta-gb";
unless(-e "filter_alignment_done"){
        die "Filter alignment failed.\n";
}else{
        system "rm -rf filter_alignment_done";
}
print LOG "OTU alignment and filtration done.\n\n";


# Phylogenetic tree construction
print LOG "Constructing phylogentic tree...\n";
system "set -ex; make_phylogeny.py -i mafft_aligned_seqs/otus_resampled_aligned.fasta-gb -o tree.nwk; touch tree_done";
unless(-e "tree_done"){
        die "Tree failed.\n";
}else{
        system "rm -rf tree_done";
}
print LOG "Phylogenetic tree done.\n\n";

chdir "..";



# Make OTU network
print LOG "Making OTU network...\n";
system "rm -rf OTU_Network";
system "set -ex; make_otu_network.py -m $map -i OTU/otu_table_resampled.biom -o OTU_Network; touch network_done";
unless(-e "network_done"){
        die "Network failed.\n";
}else{
        system "rm -rf network_done";
}
print LOG "Making network done.\n\n";




# Make taxonomic summary chart
print LOG "Making taxonomic summary chart...\n";
system "rm -rf community_structure";
system "set -ex; summarize_taxa_through_plots.py -i OTU/otu_table_resampled.biom -o community_structure -s -m $map; touch barchart_done";
unless(-e "barchart_done"){
        die "Taxonomic summary chart failed.\n";
}else{
        system "rm -rf barchart_done";
}
print LOG "Making taxonomic summary chart done.\n\n";



# Make heatmap at phylum down to family level
print LOG "Make heatmap...\n";
chdir "community_structure";
system "rm -rf heatmap; mkdir heatmap";

for(2..5){
        $level=$_;
        system "perl -pi -e 's/# Constructed.*?\n//' otu_table_resampled_sorted_L$level.txt";
        system "Rscript $home/bin/heatmap.R otu_table_resampled_sorted_L$level.txt";
        system "mv heatmap.pdf heatmap/heatmap_level$level.pdf";
}

chdir "..";



# Alpha diversity
print LOG "Alpha diversity...\n";
system "rm -rf alpha_diversity; echo \"alpha_diversity:metrics shannon,PD_whole_tree,chao1,observed_species\" > alpha_params.txt";
system "set -ex; alpha_rarefaction.py -i OTU/otu_table_with_tax_prok.biom -m $map -o alpha_diversity -p alpha_params.txt -t OTU/tree.nwk -n 50 -a -O 24; touch alpha_done";
system "rm -rf alpha_params.txt";
unless(-e "alpha_done"){
        die "Alpha diversity failed.\n";
}else{
        system "rm -rf alpha_done";
}
system "set -ex; alpha_diversity.py -i OTU/otu_table_resampled.biom -o alpha_diversity/Alpha_diversity_indexes.txt -m PD_whole_tree,chao1,observed_otus,shannon,simpson,dominance,goods_coverage -t OTU/tree.nwk; touch index_done";
unless(-e "index_done"){
        die "Alpha index failed.\n";
}else{
        system "rm -rf index_done";
}

print LOG "Alpha diversity done.\n\n";




# Beta diversity and plots
print LOG "Beta diversity...\n";
system "rm -rf beta_diversity; mkdir beta_diversity";
system "echo 'beta_diversity:metrics  bray_curtis,euclidean,unweighted_unifrac,weighted_unifrac' > beta_diversity/parameters.txt";
system "set -ex; nohup beta_diversity_through_plots.py -p beta_diversity/parameters.txt -i OTU/otu_table_resampled.biom -f -m $map -o beta_diversity -t OTU/tree.nwk > beta_diversity/log; touch beta_done";

unless(-e "beta_done"){
        die "Beta diversity failed.\n";
}else{
        system "rm -rf beta_done";
}


# 2D plot
opendir DIR,"beta_diversity";
@dir=readdir DIR;
foreach(@dir){
        if(/(.*?)_pc.txt/){
                $matrix=$1;
                system "make_2d_plots.py -i beta_diversity/$_ -o beta_diversity/2D_plot_$matrix -m $map";
        }
}


# Jackknifed
system "set -ex; nohup jackknifed_beta_diversity.py -i OTU/otu_table_resampled.biom -t OTU/tree.nwk -m $map -e $depth -o beta_diversity/jackknife >> beta_diversity/log; touch jack_done";
unless(-e "jack_done"){
        die "Jackknifed beta diversity failed.\n";
}else{
        system "rm -rf jack_done";
}

# Dendrogram
system "make_bootstrapped_tree.py -m beta_diversity/jackknife/unweighted_unifrac/upgma_cmp/master_tree.tre -s beta_diversity/jackknife/unweighted_unifrac/upgma_cmp/jackknife_support.txt -o beta_diversity/jackknife/unweighted_unifrac/upgma_cmp/jackknife_named_nodes.pdf";
system "make_bootstrapped_tree.py -m beta_diversity/jackknife/weighted_unifrac/upgma_cmp/master_tree.tre -s beta_diversity/jackknife/weighted_unifrac/upgma_cmp/jackknife_support.txt -o beta_diversity/jackknife/weighted_unifrac/upgma_cmp/jackknife_named_nodes.pdf";
print LOG "Beta diversity done.\n\n\nPlease see the document for Materials and Methods in manuscript.\n Please use the fastq file in split/ to submit to NCBI etc.\nEnjoy it.\n\n";
system "perl $home/bin/summarize.pl";
system "cp $home/materials_and_methods.txt .";


