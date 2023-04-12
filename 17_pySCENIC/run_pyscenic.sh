pwd
# 不同物种的数据库不一样，这里是人类是human 
dir=/data/nas1/luchunlin/project/BJTC-406-12/17_pySCENIC #改成自己的目录
tfs=$dir/hs_hgnc_tfs.txt
feather=$dir/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl 
# 一定要保证上面的数据库文件完整无误哦
input_loom=pbmc_all.loom
ls $tfs  $feather  $tbl  

# pyscenic 的3个步骤之 grn
pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
$input_loom  $tfs 

#pyscenic 的3个步骤之 cistarget
pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 20  \
--mask_dropouts

#pyscenic 的3个步骤之 AUCell
pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 20 