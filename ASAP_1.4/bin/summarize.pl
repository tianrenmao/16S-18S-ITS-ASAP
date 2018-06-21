#! perl -w
open FI,">summary.txt";

open F,"log.txt";
while(<F>){
	if(/OTUs, (\d+) chimeras\r?\n/){
		$chi=$1;
	}elsif(/In total (\d+) sequences assigned to samples/){
		$num=$1;
	}elsif(/In total (\d+) OTUs acquired/){
		$ini=$1;
	}elsif(/In total (\d+) OTUs after resampling./){
		$resotu=$1;
	}elsif(/Resampling depth will be (\d+)\./){
		$dep=$1;
	}elsif(/Total count: (\d+)/){
		$tot=$1;
	}elsif(/Min: (\d+)/){
		$min=$1;
	}elsif(/Max: (\d+)/){
		$max=$1;
	}elsif(/Median: (\d+)/){
		$median=$1;
	}elsif(/Mean: (\d+)/){
		$mean=$1;
	}
}

system "wc OTU/otu_table_with_tax.txt > num; wc OTU/otu_table_with_tax_prok.txt >> num";
open F, "num";
@file=<F>;
$file=join'',@file;
if($file=~/^ *(\d+).*?\r?\n?^ *(\d+)/sgim){
	$chl=$1-$2;
}
system "rm -rf num";

print FI "Total sequences assigned to samples: $num\nChimeric sequences removed: $chi\nsequences after removing chimeric and singlton: $tot\nread summary (see OTU/read_summary.txt for detail)\nMin: $min\nMax: $max\nMedian: $median\nMean: $mean\nInitial OTU number: $ini\nOTU of Chloroplast/Mitochondria removed: $chl\nResampling depth: $dep\nResampled OTU: $resotu\n\nPlease see the fastq file in split/ which is ready for submission to NCBI.\nPlease see the materials_and_methods.txt for manuscript writting.\nPlease see screen output or nohup.out for commands used.\n";

