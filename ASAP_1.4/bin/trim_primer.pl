#! perl -w
# This script trims the primer and any added bases outside the primers.
open F,"$ARGV[0]"; # read the map file and get the forward and reverse primers
while(<F>){
	if(/^#/){

	}else{
		@line=split/\t/,$_;
		$fp=$line[2];
		$rp=$line[3];
		last;
	}
}

while($fp=~/(.)/g){
        $base=$1;
	if($base=~/A|T|C|G/i){
		$fp1.=$base;
	}else{
		$fp1.=".";
	}
}


$rp=reverse $rp;
$rp=~tr/ATCGatcg/TAGCtagc/;
while($rp=~/(.)/g){
        $base=$1;
        if($base=~/A|T|C|G/i){
                $rp1.=$base;
        }else{
                $rp1.=".";
        }
}

$fp1_rc=reverse $fp1;
$fp1_rc=~tr/ATCGatcg/TAGCtagc/;
$rp1_rc=reverse $rp1;
$rp1_rc=~tr/ATCGatcg/TAGCtagc/;

open F,"$ARGV[1]"; # read the sequence file
while(<F>){
	$_=~s/[\r\n]//g;
	if(/>/){
		$id=$_;
	}else{
		$seq{$id}.=$_;
	}
}

open F,">$ARGV[2]";
foreach(keys%seq){
	$id=$_;
	$seq=$seq{$id};
	if($seq=~/^.{0,20}$fp1(.*?)$rp1.{0,20}$/si){
		print F "$id\n$1\n";
	}elsif($seq=~/^.{0,20}$rp1_rc(.*?)$fp1_rc.{0,20}$/si){
		print F "$id\n$1\n";
	}
}





