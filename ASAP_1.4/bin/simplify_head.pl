#! perl -w
open F,"$ARGV[0]";
while(<F>){
	if(/>(.*?)_\d+(.*)/){
		$sample{$1}++;
		print ">$1\_$sample{$1}$2\n";
	}else{
		print $_;
	}
}
