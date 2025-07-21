# compare inferred backbone tree distances

require(phangorn);

# full has every gene (n = 53) modelled with GTR+G
# part has optimal partitioning strategy with idiosyncratic substitution models
# - both optimal partitioning strategies have 24 partitions, but different constituents and models
# gappy contains a lot of empty sites (total: 100764)
# clean has these removed (total: 78795)
clean_full <- read.nexus("backbone_53genes_cleaned_ML_rooted.tre");
clean_part <- read.nexus("backbone_53genes_cleaned_partfind_ml_rooted.tre"); # this should probably be preferred
gappy_full <- read.nexus("backbone_53genes_gappy_ML_rooted.tre");
gappy_part <- read.nexus("backbone_53genes_gappy_partfind_ml_rooted.tre");

# clean vs. gappy
treedist(clean_full, gappy_full);
# 60
treedist(clean_part, gappy_part);
# 26. much closer

# clean: influence of optimal partitioning
treedist(clean_full, clean_part);
# 34

# gappy: influence of optimal partitioning
treedist(gappy_full, gappy_part);
# 56