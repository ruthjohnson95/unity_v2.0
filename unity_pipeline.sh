#!/usr/bin/env sh 

input_yaml=$1 

# parse arguments 
TRAIT_1=$(cat $input_yaml | grep "TRAIT 1" | cut -d':' -f2 | tr -d " \t\n\r" )
TRAIT_2=$(cat $input_yaml | grep "TRAIT 2" | cut -d':' -f2 | tr -d " \t\n\r" )
GWAS_FILE_1=$(cat $input_yaml | grep "GWAS 1" | cut -d':' -f2 | tr -d " \t\n\r" )
GWAS_FILE_2=$(cat $input_yaml | grep "GWAS 2" | cut -d':' -f2 | tr -d " \t\n\r" )
N1=$(cat $input_yaml | grep "N1" | cut -d':' -f2 | tr -d " \t\n\r" )
N2=$(cat $input_yaml | grep "N2" | cut -d':' -f2 | tr -d " \t\n\r" )
SNP_HEAD=$(cat $input_yaml | grep "SNP HEAD" | cut -d':' -f2 | tr -d " \t\n\r" )
FREQ_HEAD=$(cat $input_yaml | grep "FREQ HEAD" | cut -d':' -f2 | tr -d " \t\n\r" )
Z_HEAD=$(cat $input_yaml | grep "Z HEAD" | cut -d':' -f2 | tr -d " \t\n\r" )
N_HEAD=$(cat $input_yaml | grep "N HEAD" | cut -d':' -f2 | tr -d " \t\n\r" )
MAF_THRESH=$(cat $input_yaml | grep "MAF THRESH" | cut -d':' -f2 | tr -d " \t\n\r" )
BP_HEAD=$(cat $input_yaml | grep "BP HEAD" | cut -d':' -f2 | tr -d " \t\n\r" )
WINDOW=$(cat $input_yaml | grep "LD WINDOW" | cut -d':' -f2 | tr -d " \t\n\r" )
LDSC_DIR=$(cat $input_yaml | grep "LDSC DIR" | cut -d':' -f2 | tr -d " \t\n\r" )
HM3_PATH=$(cat $input_yaml | grep "HM3 PATH" | cut -d':' -f2 | tr -d " \t\n\r" )
REF_LD_PATH=$(cat $input_yaml | grep "REF LD PATH" | cut -d':' -f2 | tr -d " \t\n\r" )
SEED=$(cat $input_yaml | grep "SEED" | cut -d':' -f2 | tr -d " \t\n\r" )

# make results directory for processed sumstats 
SUMSTATS_DIR=results/$TRAIT_1'_'$TRAIT_2
mkdir -p $SUMSTATS_DIR

# filter by MAF, covert zscores to betas, overlap sumstats 
Rscript scripts/preprocess_sumstats.R --gwas_file_1 $GWAS_FILE_1 --gwas_file_2 $GWAS_FILE_2 --snp_head $SNP_HEAD --freq_head $FREQ_HEAD --z_head $Z_HEAD --N_head $N_HEAD --trait_1 $TRAIT_1 --trait_2 $TRAIT_2 --maf_thresh $MAF_THRESH --outdir $SUMSTATS_DIR

# perform LD pruning 
OVERLAP_GWAS_FILE_1=$SUMSTATS_DIR'/'$TRAIT_1'_'$TRAIT_2'.'processed
OVERLAP_GWAS_FILE_2=$SUMSTATS_DIR'/'$TRAIT_2'_'$TRAIT_1'.'processed
python scripts/prune_ld.py --gwas_file $OVERLAP_GWAS_FILE_1 --window $WINDOW --bp_head $BP_HEAD
python scripts/prune_ld.py --gwas_file $OVERLAP_GWAS_FILE_2 --window $WINDOW --bp_head $BP_HEAD

# run CT-LDSC 
bash scripts/run_ldsc.sh $LDSC_DIR $HM3_PATH $REF_LD_PATH $TRAIT_1 $TRAIT_2 $OVERLAP_GWAS_FILE_1 $OVERLAP_GWAS_FILE_2 $N1 $N2 $SUMSTATS_DIR 

# parse ldsc files 
H1=$(cat $SUMSTATS_DIR/'ldsc'.$TRAIT_1'_'$TRAIT_2.$TRAIT_2'_'$TRAIT_1.log | grep "Total Observed scale h2:" | head -n 1 | cut -d':' -f2 | cut -d'(' -f1 | awk '$1=$1')
H2=$(cat $SUMSTATS_DIR/'ldsc'.$TRAIT_1'_'$TRAIT_2.$TRAIT_2'_'$TRAIT_1.log | grep "Total Observed scale h2:" | tail -n 1 | cut -d':' -f2 | cut -d'(' -f1 | awk '$1=$1')
RHO=$(cat $SUMSTATS_DIR/'ldsc'.$TRAIT_1'_'$TRAIT_2.$TRAIT_2'_'$TRAIT_1.log | grep "Genetic Correlation:" | cut -d':' -f2 | cut -d'(' -f1 | awk '$1=$1' )

# run UNITY 
python src/main.py  --file1 $OVERLAP_GWAS_FILE_1  --file2 $OVERLAP_GWAS_FILE_2 --N1 $N1 --N2 $N2 --id $TRAIT_1'_'$TRAIT_2 --s $SEED --H1 $H1 --H2 $H2 --rho $RHO --outdir $SUMSTATS_DIR





