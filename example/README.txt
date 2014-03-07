# An example of running LDhat and LDhot.
# The following document lists the steps used to generate the example output shown in this folder.
# Note that these steps can be run using the run_example.pl perl script.
#
# The example generates a test dataset containing a large, central hotspot.
#
# Step 1. Simulate some data using the LDhat fin program.
#
/ldhat_path/fin  -nsamp 20 -len 5000 -theta 0.02 -r fin_rmap.txt -i -p -L 100 -s 1 -prefix example-
#
# Step 2. Estimate recombination rates using the LDhat interval program.
#
/ldhat_path/interval -seq example-sim.seq -loc example-sim.loc -lk lk_n20_t0.001 -bpen 5 -its 1500000 -samp 1000 -concise -seed 1 -prefix example-
#
# Step 3. Summarize the interval output to generate a 'res' file.
#
/ldhat_path/stat -input example-rates.txt -burn 500 -loc example-sim.loc -prefix example-
# 
# Step 4. The rate file can be large, but is no longer needed, and so can be deleted
#
rm example-rates.txt
#
# Step 5. Run my_ldhot
#
/ldhot_path/my_ldhot --seq example-sim.seq --loc example-sim.loc --lk lk_n20_t0.001 --res example-res.txt --nsim 1000 --seed 1 --out example
#
# Step 6. Summarize the my_lhot_output
#
/ldhot_path/my_ldhot_summary --res example-res.txt --hot example.hotspots.txt --out example
#
# Done!
#
