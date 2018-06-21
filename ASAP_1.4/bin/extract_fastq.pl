#! perl -w
opendir D,".";
@dir=readdir D;
foreach(@dir){
	if(/dataset_\d+/){
		open F,"$_/read_merge/merged.assembled.fastq";
		while(<F>){
	        $id=$_;
	        $id=(split/[\r\n ]/,$id)[0];
	        $seq=<F>;
	        $plus=<F>;
	        $qual=<F>;
	        $seq{$id}=$seq;
	        $qual{$id}=$qual;
		}

	}
}
open FI,">split/seqs_trim_primer.fastq";
open F,"split/seqs_trim_primer.fna";
while(<F>){
	$id=(split/[\r\n ]/,$_)[1];
	$head=$_;
	$head=~s/>/\@/;
	$seq=<F>;
	$seq=~s/[\r\n]//g;
	$id="\@$id";
	
	if($seq{$id}=~/(.*?)$seq/){
		$len1=length$1;
		$len2=length$seq;
		$qual=$qual{$id};
		$qual=substr($qual,$len1,$len2);
		print FI "$head$seq\n+\n$qual\n";
	}
}
