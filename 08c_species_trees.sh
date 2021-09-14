# first, concatenate all the trees from raxml into a single file (certhia_50kbp.trees or certhia_100kbp.trees)

# maximum clade credibility tree (simple summarization) from all gene trees using dendropy:
# gives info about which nodes have support from what proportion of gene trees
sumtrees.py --output=certhia_50kbp.tre --min-clade-freq=0.05 certhia_50kbp.trees 
sumtrees.py --output=certhia_100kbp.tre --min-clade-freq=0.05 certhia_100kbp.trees 

# coalescent tree of all gene trees using ASTRAL III
# automatically calculates local branch support using quartets, described here: https://doi.org/10.1093/molbev/msw079
java -jar ~/Astral/astral.5.6.3.jar -i certhia_50kbp.trees -o certhia_50kbp_astral.tre 2> certhia_50kbp_astral.log
java -jar ~/Astral/astral.5.6.3.jar -i certhia_100kbp.trees -o certhia_100kbp_astral.tre 2> certhia_100kbp_astral.log
