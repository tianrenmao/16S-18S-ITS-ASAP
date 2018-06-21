#! perl -w
# This script extract barcode fq based on the order of a sequence fq file (usually a merged fq).

open F,"$ARGV[0]"; # barcode fq file
while(<F>){
	$id1=$_;
	$id2=<F>;
	$id3=<F>;
	$id4=<F>;
	$bar{$id1}=$id1.$id2.$id3.$id4;
}

open F,"$ARGV[1]"; # sequence fq file
$bar_file=$ARGV[0];
$bar_file=~s/\.fastq$|\.fq$//;
open FI,">$bar_file\_extracted_barcode.fastq";
while(<F>){
        $id1=$_;
        $id2=<F>;
        $id3=<F>;
        $id4=<F>;
	print FI "$bar{$id1}";
}
