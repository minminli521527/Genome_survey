rm(list = ls())
setwd("C:/Users/99039/Desktop")

#导入 ggplot2 包
library(ggplot2)

#读取文件
depth_base <- read.delim('depth_base.stat.txt')

#(1)
#ggplot2 散点图
depth_GC <- ggplot(depth_base, aes(GC, depth)) +
  geom_point(color = 'gray', alpha = 0.6, pch = 19, size = 0.5) +
  theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(x = 'GC %', y = 'Depth')
#以输入文件中的GC列为横坐标，depth列为纵坐标，绘制散点图（geom_point()）。见要设置点的颜色为灰色并半透明，设置背景框样式，添加坐标轴标签。输出作图结果“depth_GC_plot.png”。
ggsave('depth_GC_plot.png', depth_GC, width = 5.5, height = 5.5)

#由于PCR的偏好性扩增等影响，导致基因组中个别序列出现了非特异性扩增，测序深度远大于平均值。因此图中出现了纵坐标较大的偏离值点。由于这些序列所占比例极少，在作图时可考虑将其过滤掉。
#统计平均测序深度（百分比）
depth_mean <- round(mean(depth_base$depth), 2)
#，并屏蔽掉11倍于平均测序深度的滑窗区间。使用剩余的滑窗序列的统计信息进行作图展示。
depth_base <- subset(depth_base, depth <= 11 * depth_mean)


#统计平均 GC 含量（百分比）
GC_mean <- round(mean(depth_base$GC), 2)
#ggplot2 散点图
depth_GC <- ggplot(depth_base, aes(GC, depth)) +
  geom_point(color = 'gray', alpha = 0.6, pch = 19, size = 0.5) +
  theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_vline(xintercept = GC_mean, color = 'red', lty = 2, lwd = 0.5) + 
  geom_hline(yintercept = depth_mean, color = 'red', lty = 2, lwd = 0.5) +
  labs(x= paste('GC % (Average :', GC_mean, '%)'), y = paste('Depth (Average :', depth_mean, 'X)')) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12))
ggsave('depth_GC_plot.png', depth_GC, width = 5.5, height = 5.5)

#对于正常的组装基因组序列来讲，无论GC含量还是测序深度，都会在某一区域大量集中。对于上图来讲，虽然能够看到散点主要分布的区域，但效果不很明显。因此我们可考虑在图中点的集中区域，使用易于区分的颜色展示点的密度。
#添加颜色密度
#在前述作图结果“depth_GC”样式的基础上，使用stat_density2d()命令，根据散点在图中的密度分布继续添加颜色密度标记，并使用scale_fill_gradientn()指定颜色梯度，根据散点分布密度由低到高在图中标记为不同颜色。输出作图结果“depth_GC_plot.png”。。
depth_GC <- depth_GC +
  stat_density2d(aes(fill = ..density.., alpha = ..density..), geom = 'tile', contour = FALSE, n = 500, show.legend = FALSE) +
  scale_fill_gradientn(colors = c('transparent', 'gray', 'yellow', 'red'))

ggsave('depth_GC_plot.png', depth_GC, width = 5.5, height = 5.5)
#到这里，基因组Depth-GC密度分布图基本上就算完成了。通过该图可以获得我们想要得知的信息，判断组装基因组的质量，以及测序质量等。




#(2)depth 深度直方密度图
#此外，还可考虑统计滑窗GC含量，或者滑窗测序reads覆盖度，绘制二者的分布直方图，添加在散点图的两侧，以更直观展示数据分布情况。
#geom_histogram()用于绘制频数直方图，并使用binwidth参数划分为100个频数统计窗格（即无论测序深度怎样分布，均根据深度频数划分为100个“柱子”展示）。额外使用geom_rug()在频数直方图的下方添加地毯图，作为展示每个基因组滑窗序列测序深度的位置。输出作图结果“depth_hist_plot.png”。
depth_hist <- ggplot(depth_base, aes(depth)) +
  geom_histogram(binwidth = (max(depth_base$depth) - min(depth_base$depth))/100, fill = 'gray', color = 'gray40', size = 0.1) +
  geom_rug(color = 'gray', alpha = 0.6) +
  theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(axis.line = element_line(color = 'black', size = 0.3), axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
  labs(x = '', y = 'Numbers') +
  geom_vline(xintercept = depth_mean, color = 'red', lty = 2, lwd = 0.5)

ggsave('depth_hist_plot.png', depth_hist, width = 5.5, height = 2.5)




#(3)GC 含量直方密度图
GC_hist <- ggplot(depth_base, aes(GC)) +
  geom_histogram(binwidth = (max(depth_base$GC) - min(depth_base$GC))/100, fill = 'gray', color = 'gray40', size = 0.1) +
  geom_rug(color = 'gray', alpha = 0.6) +
  theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(axis.line = element_line(color = 'black', size = 0.3), axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
  labs(x = '', y = 'Numbers') +
  geom_vline(xintercept = GC_mean, color = 'red', lty = 2, lwd = 0.5)

ggsave('GC_hist_plot.png', GC_hist, width = 5.5, height = 2.5)




#(4)组合图
#最后可使用grid包，将上述已完成的散点图与两个频数直方图组合在一起，得到最终的结果。
#组合后，输出的作图结果“depth_GC_all_plot.png”。
#导入 grid 包
library(grid)
#此处首先将上述 depth 深度直方密度图旋转 90 度
depth_hist <- depth_hist + coord_flip()
#grid 组合图形
png('depth_GC_all_plot.png', width = 4000, height = 4000, res = 600, units = 'px')
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 3)))
print(depth_GC, vp = viewport(layout.pos.row = 2:3, layout.pos.col = 1:2))
print(GC_hist, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(depth_hist, vp = viewport(layout.pos.row = 2:3, layout.pos.col = 3))
dev.off()

