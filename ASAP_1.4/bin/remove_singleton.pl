#! perl -w
# Input: unique sequences from dereplication with size information

open F,"$ARGV[0]";
open FI,">$ARGV[1]";

while(<F>){
	if(/>/){
		$head=$_;
		push @id,$head;
	}else{
		$seq{$head}.=$_;
	}
}


foreach(@id){
	unless(/size=1;/){
		print FI "$_$seq{$_}";
	}
}
