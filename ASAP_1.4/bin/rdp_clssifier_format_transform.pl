#! perl -w
open F,"$ARGV[0]";
while(<F>){
	if(/(.*?)\t(.*?\n)/){
		$otu=$1;
		$tax=$2;
		$tax=~s/"//g;
		$tax=~s/^/D_0__/;
		$tax=~s/\t/; D_1__/;
		$tax=~s/\t/; D_2__/;
		$tax=~s/\t/; D_3__/;
		$tax=~s/\t/; D_4__/;
		$tax=~s/\t/; D_5__/;
		$tax=~s/\t/; D_6__/;
		print "$otu\t$tax";
	}

}
