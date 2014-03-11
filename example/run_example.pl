#!/usr/bin/perl

$narg = $#ARGV+1;
if ($narg < 2)
{
	$str = "Run ldhot example.\n\n\tUsage perl run_example.pl /path/to/ldhat/ /path/to/ldhot/\n\n";
	die($str);
}
$ldhat_path = $ARGV[0];
$ldhot_path = $ARGV[1];
$out_prefix = "example-";
if ($narg >= 3)
{
	$out_prefix = $ARGV[2];
}

$seed = 1;
if ($narg >= 4)
{
	$seed = $ARGV[3];
}

$lk_path = $ldhat_path . "lk_files/";

$fin = $ldhat_path . "fin ";
$interval = $ldhat_path . "interval ";
$stat = $ldhat_path . "stat ";
$lkgen = $ldhat_path . "lkgen ";
$ldhot = $ldhot_path . "ldhot ";
$ldhot_summary = $ldhot_path . "ldhot_summary ";

$L = 100;
$rmap = "fin_rmap.txt";
$n_samp = 20;

$bpen = 5;
$its = 1500000;
$samp = 1000;
$burn = 500;

$len = 0;
open(RMAP, '<', $rmap) or die("Could not open $rmap\n");
$header = <RMAP>;
while($line = <RMAP>)
{
    $len++;
}
close(RMAP);

$theta_per_kb = 1.0;
$theta_per_fin_site = $theta_per_kb * $L / $len;

# First simulate some data using fin
$str = $fin . " -nsamp $n_samp -len $len -theta $theta_per_fin_site -r $rmap -i -p -L $L -s $seed -prefix $out_prefix";
print "$str\n";
system($str);

$loc_file = $out_prefix . "sim.loc";
$seq_file = $out_prefix . "sim.seq";
$lk_file = "lk_n" . $n_samp . "_t0.001 ";

# Now run LDhat
$str = $interval . " -seq $seq_file -loc $loc_file -lk $lk_file -bpen $bpen -its $its -samp $samp -concise -seed $seed -prefix $out_prefix ";
print "$str\n";
system($str);

# Now run stat
$rates_file = $out_prefix . "rates.txt";
$str = $stat . " -input $rates_file -burn $burn -loc $loc_file -prefix $out_prefix ";
print "$str\n";
system($str);

# Rate file no longer needed
unlink($rates_file);

$res_file = $out_prefix . "res.txt";

# Now run LDhot
$out_prefix = substr($out_prefix, 0, -1);
$str = $ldhot . " --seq $seq_file --loc $loc_file --lk $lk_file --res $res_file --nsim 1000 --seed $seed --out $out_prefix ";
print "$str\n";
system($str);

$hot = $out_prefix . ".hotspots.txt";
$str = $ldhot_summary . " --res $res_file --hot $hot --out $out_prefix ";
print "$str\n";
system($str);
