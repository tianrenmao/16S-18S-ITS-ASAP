#! perl -w
# This script remove Chloroplast, Mitochondria and even Archaea (need revision) from input of OTU table.
open FI,">$ARGV[1]";
open F,"$ARGV[0]"; # Read the OTU table.
while(<F>){
	if(/D_2__Chloroplast/ || /Mitochondria/ ){ # If you want to remove Archaea, add "/D_0__Archaea/".
	}else{
		print FI $_;
	}
}
