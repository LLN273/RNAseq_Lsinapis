##################################################################################################################################################################### 
# 
# Pipeline used in:
#
# "Gene expression profiling across ontogenetic stages in wood white (Leptidea sinapis) reveals pathways linked to butterfly diapause regulation."
# Luis Leal, Venkat Talla, Thomas Källman, Magne Friberg, Christer Wiklund, Vlad Dincă, Roger Vila, Niclas Backström
# Mol. Eco. (2018)
#
#
##################################################################################################################################################################### 
#  
# Luis Leal
# Uppsala University, Uppsala, Sweden, 2017
#
##################################################################################################################################################################### 




# ********************************************************************************************* 
#
# Quality Control I: Filter each library using TrimGalore, v. 0.4.1 (Babraham Bioinformatics)
# 
# ********************************************************************************************* 

trim_galore \
            --quality 30 \
            --paired \
            --illumina \
            --phred33 \
            --stringency 1 \
            -e 0.1 \
            --clip_R1 12 \
            --clip_R2 12 \
            --three_prime_clip_R1 3 \
            --three_prime_clip_R2 3 \
            --length 30 \
            --gzip \
            --output_dir $OUTDIR \
            --fastqc \
            $READ1_raw \
            $READ2_raw 






# ********************************************************************************************* 
#
# Quality Control II: Mask each library using Fastq Masker (FASTX, v. 0.0.14 toolkit)
# 
# ********************************************************************************************* 

fastq_masker \
        -v \
        -q 20 \
        -r N \
        -z \
        -i $READ1_tg \
        -o $READ1_masked


fastq_masker \
        -v \
        -q 20 \
        -r N \
        -z \
        -i $READ2_tg \
        -o $READ2_masked






# ********************************************************************************************* 
#
# Quality Control III: Filter each library using prinseq, v. 0.20.4 (Schmieder & Edwards 2011) 
# 
# *********************************************************************************************

prinseq \
           -verbose \
           -fastq $READ1_masked \
           -fastq2 $READ2_masked \
           -out_good $READ_polyA_OUT \
           -out_bad $READ_REJECTED \
           -log prinseq.log \
           -graph_data prinseq.gd \
           -graph_stats ld,gc,qd, ns, pt, ts, aq, de, sc, dn \
           -qual_noscale \
           -exact_only \
           -min_len 30 \
           -ns_max_p 10 \
           -trim_tail_left 2 \
           -trim_tail_right 2 \
           -lc_method dust \
           -lc_threshold 7 \
           -trim_ns_left 1 \
           -trim_ns_right 1 \
           -trim_qual_left 30 \
           -trim_qual_right 30






# ********************************************************************************************* 
#
# Quality Control IV: Filter each library using condetri, v. 2.3 (Smeds & Künstner 2011)   
# 
# *********************************************************************************************

perl condetri_v2.3.pl \
            -fastq1=$READ_polyA_trim1 \
            -fastq2=$READ_polyA_trim2 \
            -prefix=condetri- \
            -cutfirst=0 \
            -cutlast=0 \
            -rmN \
            -hq=30 \
            -lq=0 \
            -frac=0.8 \
            -minlen=30 \
            -mh=1 \
            -ml=1 \
            -sc=33






# ********************************************************************************************* 
#
# Quality Control V: Screen each library using fastQ Screen, v. 0.9.2 (Babraham Bioinformatics) 
# Dependencies: fastq_screen_i.conf (fastq screen configuration file)   
# 
# *********************************************************************************************

fastq_screen \
                --threads 3 \
                --aligner bowtie2 \
                --conf fastq_screen_i.conf \
                --subset 100000000 \
                --filter '-00000' \
                --outdir ./ \
                --tag \
                $condetri-READ1

fastq_screen \
                --threads 3 \
                --aligner bowtie2 \
                --conf fastq_screen_i.conf \
                --subset 100000000 \
                --filter '-00000' \
                --outdir ./ \
                --tag \
                $condetri-READ2






# *********************************************************************************************** 
#
# Reconcile paired-end files after 'fastQ Screen' filtering   
# 
# ***********************************************************************************************

python3 adjust_paired_files.py $READ1_fqs_filtered $READ2_fqs_filtered










