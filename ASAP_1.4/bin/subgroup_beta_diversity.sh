#! /bin/bash
# put subgroup map files in subgroup_map

set -ex
rm -rf subgroup_beta_diversity
mkdir subgroup_beta_diversity
echo 'beta_diversity:metrics  bray_curtis,euclidean,unweighted_unifrac,weighted_unifrac' > subgroup_beta_diversity/parameters.txt
for i in `ls subgroup_map`; do beta_diversity_through_plots.py --suppress_emperor_plots -p subgroup_beta_diversity/parameters.txt -i OTU/otu_table_resampled.biom -f -m subgroup_map/$i -o subgroup_beta_diversity/`echo $i | perl -pi -e 's/.txt//'` -t OTU/tree.nwk; make_2d_plots.py -i subgroup_beta_diversity/`echo $i | perl -pi -e 's/.txt//'`/bray_curtis_pc.txt -o subgroup_beta_diversity/`echo $i | perl -pi -e 's/.txt//'`/2D_plot_bray_curtis -m subgroup_map/$i; make_2d_plots.py -i subgroup_beta_diversity/`echo $i | perl -pi -e 's/.txt//'`/euclidean_pc.txt -o subgroup_beta_diversity/`echo $i | perl -pi -e 's/.txt//'`/2D_plot_euclidean -m subgroup_map/$i; make_2d_plots.py -i subgroup_beta_diversity/`echo $i | perl -pi -e 's/.txt//'`/unweighted_unifrac_pc.txt -o subgroup_beta_diversity/`echo $i | perl -pi -e 's/.txt//'`/2D_plot_unweighted_unifrac -m subgroup_map/$i; make_2d_plots.py -i subgroup_beta_diversity/`echo $i | perl -pi -e 's/.txt//'`/weighted_unifrac_pc.txt -o subgroup_beta_diversity/`echo $i | perl -pi -e 's/.txt//'`/2D_plot_weighted_unifrac -m subgroup_map/$i; done




