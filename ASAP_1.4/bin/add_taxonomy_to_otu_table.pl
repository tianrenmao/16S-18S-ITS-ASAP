#! perl -w
open F,"$ARGV[0]"; # Read the classification.
while(<F>){
	chomp;
	@line=split/\t/,$_;
	$tax{$line[0]}=$line[1];
}

open FI,">$ARGV[2]";
open F,"$ARGV[1]"; # Read the raw OTU table.
$head=<F>;
chomp $head;
$head.="\ttaxonomy\n";
print FI $head;
while(<F>){
	chomp;
	@line=split/\t/,$_;
	if($tax{$line[0]}){
		print FI "$_\t$tax{$line[0]}\n";
	}else{
		print FI "$_\tUnclassified\n";
	}
}
