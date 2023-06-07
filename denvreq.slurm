#!/bin/bash
# file: denvreq.slurm
#SBATCH --partition=normal
#SBATCH --job-name=denvreq
#SBATCH --output=./logs_%j.out
#SBATCH --error=./logs_%j.err
#SBATCH --exclude=c001
#SBATCH -n 72
#SBATCH -N 1-1
#SBATCH --ntasks-per-node=72
#SBATCH --mem=400G 
#SBATCH --output=slurm-%j_denvreq.out
#SBATCH --error=slurm-%j_denvreq.err

# Enable hyperthreading. Set it after modules load.
export SLURM_HINT=multithread
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
export PGENV="/share/home/dw_user0001/program/denvreq_env"
export DBS=$PGENV/DBs
export BINS=$PGENV/tools
export EGGNOG_DATA_DIR=$DBS/eggNOGdb
export TRIN="/share/apps/anaconda3/2021.05/envs/denvreq/pkgs/trinity-2.8.5-h8b12597_3/opt/trinity-2.8.5"
export PATH=$PATH:/share/apps/anaconda3/2021.05/envs/denvreq/bin
echo SLURM job id : $SLURM_JOB_ID
eval "$(conda shell.bash hook)"

# run the application

echo "Time: $(date)"
start=$(date +%s)
conda activate denvreq

echo "---------- Begin Rawdata Download Job -------------"

#mkdir -p fastq_data
#cp list.txt ./fastq_data && cd fastq_data
#conda activate rna_seq
#cat list.txt | parallel -j5 "kingfisher get -r {} -m aws-cp aws-http prefetch --download-threads 5"
#rm -rf list.txt
#cd ..

echo "---------- Finish Rawdata Download Job $(date)  -------------"

echo "---------- Begin Sequence data QC $(date)-------------"

mkdir -p qc_result
cd qc_result
ls ../fastq_data/*.fq.gz | xargs -i echo nohup fastqc -o ./ --nogroup {} \& >fastqc.sh
bash fastqc.sh
cd ..

echo "---------- Finish QC Job $(date)  -------------"

echo "---------- Begin Sequence data cleanning $(date)-------------"

mkdir -p clean2qc_result
mkdir -p clean_data
ln -s $(pwd)/fastq_data/*.fq.gz $(pwd)/clean_data
cd clean_data
for n in *R1.fq.gz;do
k=${n%%R1.fq.gz}
trim_galore --fastqc --gzip \
--quality 25 \
--trim-n \
--length 25 \
--stringency 3 \
--cores $OMP_NUM_THREADS \
--paired ${k}R1.fq.gz ${k}R2.fq.gz \
--fastqc_args "--nogroup --outdir ../clean2qc_result" \
-o ./ >> cleaning.log
done
cd ..

echo "---------- Finish Sequence data cleanning $(date)-------------"

echo "---------- Begin Sequence data Assembling $(date)-------------"

mkdir -p assembly
ln -s $(pwd)/clean_data/*_val*.fq.gz $(pwd)/assembly
cd assembly
for i in *R1_val_1.fq.gz; do echo -e "${i%.R1_val_1.fq.gz}\t$i\t${i%R1_val_1.fq.gz}R2_val_2.fq.gz"; done >tmp.txt
awk -F'\t' 'NR==FNR{a[$1]=$4;}NR!=FNR && a[$1] {print a[$1]"\t" $0}' ../metadata.txt tmp.txt > samples.txt
cp ../samples2.txt ./samples2.txt

Trinity --seqType fq \
--CPU $OMP_NUM_THREADS \
--max_memory $((SLURM_MEM_PER_NODE/1024))G \
--min_contig_length 150 \
--samples_file samples2.txt \
--SS_lib_type RF

cd-hit-est -i trinity_out_dir.Trinity.fasta -c 0.95 -n 10 -M 0 -T $OMP_NUM_THREADS -d 0 -o Trinity.fasta
get_Trinity_gene_to_trans_map.pl Trinity.fasta >Trinity.fasta.gene_trans_map
cd ..

echo "---------- Finish Sequence data Assembling $(date)-------------"

echo "---------- Begin Sequence data Analyzing $(date)-------------"

mkdir -p $(pwd)/annotation
ln -s $(pwd)/assembly/Trinity.fasta $(pwd)/annotation
ln -s $(pwd)/assembly/Trinity.fasta.gene_trans_map $(pwd)/annotation
cd annotation

mkdir -p annotation_final

TrinityStats.pl Trinity.fasta >Assembly_Stat.txt
perl $TRIN/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta >longest_isoform.fasta

busco -i longest_isoform.fasta \
-l $DBS/busco_db/mammalia_odb10 \
-o buscoresult \
-m trans \
--cpu $OMP_NUM_THREADS \
-f \
--offline

TransDecoder.LongOrfs -S -t Trinity.fasta

diamond blastp --query Trinity.fasta.transdecoder_dir/longest_orfs.pep \
--db $DBS/sqlite_db/uniprot_sprot.pep \
--max-target-seqs 1 \
--outfmt 6 \
--evalue 1e-5 \
--threads $OMP_NUM_THREADS >Trinity.fasta.transdecoder_dir/longest_orfs.pep.blastp.outfmt6

hmmscan --cpu $OMP_NUM_THREADS \
--domtblout Trinity.fasta.transdecoder_dir/longest_orfs.pep.pfam.domtblout \
$DBS/hmmerdb/Pfam-A.hmm \
Trinity.fasta.transdecoder_dir/longest_orfs.pep >pfam.log

TransDecoder.Predict -t Trinity.fasta \
--retain_pfam_hits Trinity.fasta.transdecoder_dir/longest_orfs.pep.pfam.domtblout \
--retain_blastp_hits Trinity.fasta.transdecoder_dir/longest_orfs.pep.blastp.outfmt6

diamond blastx --query Trinity.fasta \
--db $DBS/sqlite_db/uniprot_sprot.pep \
--threads $OMP_NUM_THREADS \
--max-target-seqs 1 \
--outfmt 6 >blastx.outfmt6

diamond blastp --query Trinity.fasta.transdecoder.pep \
--db $DBS/sqlite_db/uniprot_sprot.pep \
--threads $OMP_NUM_THREADS \
--max-target-seqs 1 \
--outfmt 6 >blastp.outfmt6

hmmscan --cpu $OMP_NUM_THREADS \
--domtblout TrinotatePFAM.out \
$DBS/sqlite_db/Pfam-A.hmm \
Trinity.fasta.transdecoder.pep >pfam.log

$BINS/signalp-4.1/signalp -f short -n signalp.out Trinity.fasta.transdecoder.pep
$BINS/tmhmm-2.0c/bin/tmhmm --short < Trinity.fasta.transdecoder.pep >tmhmm.out
sed -i 's@/util/@/@g' `find /share/apps/anaconda3/2021.05/envs/denvreq/bin -name "RnammerTranscriptome.pl"`
RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer $BINS/rnammer_env/rnammer
mv Trinity.fasta.rnammer.gff rnammer.gff

cp $DBS/sqlite_db/Trinotate20230101.sqlite ./Trinotate.sqlite

Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map \
--transcript_fasta Trinity.fasta \
--transdecoder_pep Trinity.fasta.transdecoder.pep

Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.out
Trinotate Trinotate.sqlite LOAD_rnammer rnammer.gff
Trinotate Trinotate.sqlite report -E 10e-10 > trinotate_annotation_report.xls

cp trinotate_annotation_report.xls ./annotation_final/trinotate_annotation_report.xls

mkdir -p eggnog_out
emapper.py \
-i Trinity.fasta.transdecoder.pep \
--output_dir eggnog_out \
--itype proteins \
--cpu $OMP_NUM_THREADS \
--sensmode ultra-sensitive \
--override \
-m diamond \
--target_orthologs one2one \
--report_orthologs \
--output trinity_eggnog

cp ./eggnog_out/trinity_eggnog.emapper.annotations.xlsx ./annotation_final

conda activate greedyFAS
sed -i.bak 's|\*||g' Trinity.fasta.transdecoder.pep

interproscan.sh -appl Pfam \
-iprlookup \
-goterms \
-dp \
-pa \
-dra \
-t p \
-f tsv \
--cpu $OMP_NUM_THREADS \
-i Trinity.fasta.transdecoder.pep \
-T interproscan_output

conda activate denvreq
cp Trinity.fasta.transdecoder.pep.tsv ./annotation_final/interproscan_annotation_report.tsv

cd annotation_final
 
cut -f 1,8,9,11-14 trinotate_annotation_report.xls | sed "s/#gene/gene/g" >Annotation.txt

cut -f 1,4 Annotation.txt|awk '{if(length($2)>1){print $0}}' >eggNOG.anno
if [[ ! -s eggNOG.anno ]];then
perl $PGENV/SRC/matchEggNOG.pl eggNOG.anno >eggNOG.anno.plotdata
Rscript $PGENV/SRC/simpleBarplot.R eggNOG.anno.plotdata 6 12
fi

perl $PGENV/SRC/trino2pfam.pl Annotation.txt >PFAM.anno

perl $PGENV/SRC/trino2go.pl Annotation.txt >GO.anno
perl $PGENV/SRC/go_classification_v3.pl GO.anno GO.anno.output

cut -f 1,5 Annotation.txt|awk '{if(length($2)>1){print $0}}' >Kegg.anno
#perl $PGENV/SRC/matchKegg.pl Kegg.anno >Kegg.anno.plotdata
#perl $PGENV/SRC/kegg_classification_v3.pl Kegg.anno.plotdata Kegg.anno.output
#perl $PGENV/SRC/formatAnno.pl Annotation.txt| perl $PGENV/SRC/addAnno.pl eggNOG.anno.plotdata Kegg.anno.plotdata - >Annotation.xls

cd ../../

echo "---------- Finish Sequence data Analyzing $(date)-------------"

echo "---------- Begin Sequence data quantification $(date)-------------"

mkdir -p qnt_result
cd assembly

##Quantification
align_and_estimate_abundance.pl \
--transcripts Trinity.fasta \
--samples_file samples.txt \
--seqType fq \
--SS_lib_type RF \
--est_method RSEM \
--aln_method bowtie2 \
--thread_count $OMP_NUM_THREADS \
--prep_reference \
--gene_trans_map Trinity.fasta.gene_trans_map

ls */RSEM.isoforms.results >isoforms.quant_files.txt

abundance_estimates_to_matrix.pl \
--est_method RSEM \
--gene_trans_map Trinity.fasta.gene_trans_map \
--cross_sample_norm TMM \
--name_sample_by_basedir \
--quant_files isoforms.quant_files.txt

ln -s $(pwd)/RSEM* $(pwd)/../qnt_result
cp $(pwd)/samples.txt $(pwd)/../qnt_result
cp $(pwd)/Trinity.fasta $(pwd)/../qnt_result
cp $(pwd)/Trinity.fasta.gene_trans_map $(pwd)/../qnt_result
cd ../qnt_result

ln -s $PGENV/SRC/obo $(pwd)/
cp ../annotation/GO.anno ./
$TRIN/util/misc/fasta_seq_length.pl Trinity.fasta > Trinity.fasta.seq_lens
$TRIN/util/misc/TPM_weighted_gene_length.py --gene_trans_map Trinity.fasta.gene_trans_map \
--trans_lengths Trinity.fasta.seq_lens \
--TPM_matrix RSEM.isoform.TMM.EXPR.matrix >Trinity.gene_lengths.txt

extract_GO_assignments_from_Trinotate_xls.pl \
--Trinotate_xls ../annotation/trinotate_annotation_report.xls \
-G --include_ancestral_terms >go_annotations.txt

Trinotate_get_feature_name_encoding_attributes.pl ../annotation/trinotate_annotation_report.xls >annot_feature_map.txt
$TRIN/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl RSEM.gene.counts.matrix annot_feature_map.txt >RSEM.gene.counts.wAnnot.matrix

## DE_analysis
run_DE_analysis.pl \
--matrix RSEM.gene.counts.matrix \
--method DESeq2 \
--samples_file samples.txt \
--contrasts ../contrasts.txt

## QC && Enrichment
cd DESeq2.*.dir

## QC Samples
# PCA analysis amount samples
PtR --matrix ../RSEM.gene.TMM.EXPR.matrix \
--samples ../samples.txt \
--log2 --CPM --prin_comp 2 --min_rowSums 10

# Compare gene counts across samples
PtR --matrix ../RSEM.gene.TMM.EXPR.matrix \
--min_rowSums 10 \
-s ../samples.txt --log2 --CPM --sample_cor_matrix

# Compare gene counts between biological samples
PtR --matrix ../RSEM.gene.TMM.EXPR.matrix \
--samples ../samples.txt --log2 --CPM \
--min_rowSums 10 \
--compare_replicates

mkdir -p qntqcresult && mv *TMM.EXPR.matrix.* *rep_compare.pdf qntqcresult

# Enrichment

#awk '($4>1 || $4<-1) && $7<0.05 && NR>1{print $1}' *.DE_results >DESEQ.DE.lst
# Set fold difference to 2-fold (ie. -C 1 = 2^1)
# P value <= 0.05

analyze_diff_expr.pl \
-P 0.05 -C 1 \
--samples ../samples.txt \
--matrix ../RSEM.gene.TMM.EXPR.matrix \
--examine_GO_enrichment \
--GO_annots ../go_annotations.txt \
--gene_lengths ../Trinity.gene_lengths.txt \
--max_genes_clust 50000

for k in diffExpr.*.matrix.RData; do
define_clusters_by_cutting_tree.pl --Ktree 8 -R $k
done

echo "---------- Finish Sequence data quantification $(date)-------------"

echo "Time: $(date)"
end=$(date +%s)
runtime=$((end-start))
echo "Runtime was $runtime"

