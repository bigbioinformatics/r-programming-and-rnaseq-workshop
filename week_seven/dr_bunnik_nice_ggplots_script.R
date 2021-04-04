# set text element for axis labels
text <- element_text(size = 14)

# plot for paired observations
p <- ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) + geom_point(size = 4) + geom_line(size = 1,linetype = 2)  + theme_classic()
p <- p + scale_color_manual(name = "Cell line", values=c("orangered", "limegreen", "deepskyblue", "deeppink"))
p <- p + theme(axis.title = text) + theme(axis.text.x = element_text(color="black", size = 12), axis.text.y = element_text(color="black", size = 12))
p <- p + scale_x_discrete(name="Treatment", limits=c("untrt", "trt"), labels=c("None", "Dex")) + scale_y_log10(name="Normalized count") + ggtitle("ENSG00000127954")
p
ggsave(file="paired_dots.png", height=5, width=5)

# boxplot for unpaired observations
p <- ggplot(geneCounts, aes(x = dex, y = count, fill = dex)) + geom_boxplot(outlier.shape = NA) + theme_classic() + geom_point(size = 3, position=position_jitterdodge(),alpha=0.5)
p <- p + scale_fill_manual(name = "Treatment", limits=c("untrt", "trt"), labels = c("None", "Dex"), values=c("orangered", "deepskyblue"))
p <- p + theme(axis.title = text) + theme(axis.text.x = element_text(color="black", size = 12), axis.text.y = element_text(color="black", size = 12))
p <- p + scale_x_discrete(name="Treatment", limits=c("untrt", "trt"), labels=c("None", "Dex")) + scale_y_log10(name="Normalized count")  + ggtitle("ENSG00000127954")
p
ggsave(file="boxplot_with_dots.png", height=5, width=5)
