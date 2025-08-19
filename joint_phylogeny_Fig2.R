# Script for generating the joint phylogenetic tree on Figure 2 of Chacón-Duque et al. 2025
# Requires: BEAST tree output and regions.txt file with sample metadata
# Output: Annotated tree figures (PNG)

# Load required libraries for data manipulation and tree visualization
library(tidyverse)
library(ape)
library(ggtree)
library(treeio) 
library(tidytree)
library(ggnewscale)

# Read BEAST tree output (phylogenetic tree with posterior probabilities)
beasttree <- read.beast("Chacon-Duque-2025.nodating.100M.combined.trees.out")

# Read sample metadata (regions, tip labels, etc.)
regions <- read.delim("regions.txt")

# Check that all tip labels in metadata are present in the tree
regions$tip.label %in% beasttree@phylo$tip.label
# and this to make sure that there's no missing data (list missing tips)
beasttree@phylo$tip.label[!regions$tip.label %in% beasttree@phylo$tip.label]

# Visualize node labels to identify clades and nodes of interest
ggtree(beasttree) + geom_text(aes(label=node), hjust=-.3)

# Exclude elephant outgroup by subsetting tree at node 248
treesub <- tree_subset(beasttree,248,levels_back = 0,root_edge=FALSE)

# Re-check node labels after subsetting
ggtree(treesub) + geom_text(aes(label=node), hjust=-.3)

# List of nodes to annotate with posterior probabilities (update as needed)
nodes <- c("226","227","228","229","230","381","432",
           "433","434","435","436","437","445","446","447",
           "448","449","382","383","384","385","386",
           "387","388","390","391","392","399","401",
           "402","403","404","231","232","321","322",
           "323","380","326","327","325","351","352",
           "353","354","361","362","364","233","234",
           "236","309","237","293","238")

# List of tips to highlight in the plot (update as needed)
tips_to_highlight <- c("P037_M.tro_NAS_ND", "P035_M.tro_NAS_ND", 
                       "P045_M.sp_NAS_ND","P033_M.tro_NAS_ND",
                       "P043_M.pri_NAS_ND","L171_M.pri_NAS_ND",
                       "L169_M.pri_NAS_ND","L168_M.pri_NAS_ND",
                       "P048_M.pri_NAS_ND","MD228_M.pri_NNA_ND",
                       "L172_M.pri_NAS_ND","EID18_M.sp_NNA_ND",
                       "GM20_M.pri_NAS_ND","MD227_M.pri_NAS_ND",
                       "SP1145_M.pri_NAS_ND","TB161_M.pri_NAS_ND",
                       "L263_M.pri_NAS_ND")

tips_to_plot <- c("P037_M.tro_NAS_ND","P035_M.tro_NAS_ND","P043_M.pri_NAS_ND",
                  "P045_M.sp_NAS_ND","P033_M.tro_NAS_ND")

# Define color palette for regions
cols <- c("Siberia" = "#d73027", "North America" = "#2166ac", "Europe" = "#5aae61")

# Build the main tree plot with time axis, clade labels, and highlights
test_plot <- ggtree(treesub) %>% flip(382, 432) %<+% regions + 
  # Add vertical line for population bottleneck
  geom_vline(xintercept = -281000, colour= "#636363", linetype = "dashed", alpha=.5) + 
  # Annotate and colour time periods (Pleistocene, etc.)
  annotate("text", x=-316000, y=180, label="Population bottleneck", angle=90, colour= "#636363") +
  annotate("rect", xmin = -122000, xmax = 4000, ymin = 0, ymax = 225, 
           fill='white', alpha = .3) +
  annotate("rect", xmin = -776000, xmax = -122000, ymin = 0, ymax = 225, 
           fill='#e0e0e0', alpha = .3) +
  annotate("text", x = -546000, y = 222, label = "Middle Pleistocene") + 
  annotate("rect", xmin = -776000, xmax = -1746000, ymin = 0, ymax = 225, 
           fill='white', alpha = .3) + 
  annotate("text", x = -1610000, y = 222, label = "Early Pleistocene")

test_plot <- revts(test_plot)

height95 = t(matrix(unlist(test_plot$data$height_0.95_HPD),nrow=2))
height95 = height95 - 4000 #min offset of the text

bar_height95 = as.data.frame(height95) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(test_plot$data, y))

test_plot <- test_plot + 
# Add 95% HPD bars for node ages
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_height95, 
               color="#969696",
               alpha = 0.5,
               size = 2)

test_plot <- test_plot + 
  # Highlight selected tips and add region colors
  geom_tippoint(aes(subset=(label %in% tips_to_highlight),
                    color=Region), shape=17, size = 2) + 
  scale_colour_manual(values = cols) + 
  geom_text2(aes(subset=(label %in% tips_to_plot),label = shortlabel),
             size=3,nudge_y=2,angle=90,hjust = 0) +
  theme_tree2(legend.position=c(0.2, 0.6),
              legend.background = element_blank(),
              legend.box.background = element_rect(fill = "white", color = "black")) +
  labs(color = "Deep-time specimens (by region)", 
       x = "Millions of years before present (Ma)") + 
  scale_x_continuous(breaks = c(-1746000,-1496000,-1246000,-996000,-746000,-496000,-246000,4000),
                     labels = c (1.75,1.50,1.25,1.00,0.75,0.50,0.25,0)) + 
  geom_point2(aes(subset=(node %in% nodes), 
                fill=cut(posterior, c(0, 0.7, 0.9, 1))), 
            shape=21, size=1.5) +
  scale_fill_manual(values=c("black", "#bdbdbd", "white"), guide='legend', 
                    name='Node support (posterior probability)', 
                    breaks=c("(0.9,1]", "(0.7,0.9]", "(0,0.7]"), 
                    labels=c("Higher than 0.9","Between 0.7 and 0.9","Lower than 0.7"),
                    limits =c("(0.9,1]", "(0.7,0.9]", "(0,0.7]")) +
  # Add clade labels
  geom_cladelabel(node=382, label="Clade 3", colour="#636363",
                  align=TRUE,offset=4000,offset.text=30000,hjust=0.5,angle=90) + 
  geom_cladelabel(node=432, label="Clade 2", colour="#636363",
                  align=TRUE,offset=4000,offset.text=30000,hjust=0.5,angle=90) +
  geom_cladelabel(node=231, label="Clade 1", colour="#636363",
                  align=TRUE,offset=4000,offset.text=30000,hjust=0.5,angle=90)

test_plot

###saving the plot
ggsave(plot=test_plot, file="figure2_v2.png", width=6.5, height=8, dpi=600)
