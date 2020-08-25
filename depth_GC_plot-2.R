rm(list = ls())
setwd("C:/Users/99039/Desktop")

#���� ggplot2 ��
library(ggplot2)

#��ȡ�ļ�
depth_base <- read.delim('depth_base.stat.txt')

#(1)
#ggplot2 ɢ��ͼ
depth_GC <- ggplot(depth_base, aes(GC, depth)) +
  geom_point(color = 'gray', alpha = 0.6, pch = 19, size = 0.5) +
  theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(x = 'GC %', y = 'Depth')
#�������ļ��е�GC��Ϊ�����꣬depth��Ϊ�����꣬����ɢ��ͼ��geom_point()������Ҫ���õ����ɫΪ��ɫ����͸�������ñ�������ʽ�������������ǩ�������ͼ�����depth_GC_plot.png����
ggsave('depth_GC_plot.png', depth_GC, width = 5.5, height = 5.5)

#����PCR��ƫ����������Ӱ�죬���»������и������г����˷��������������������Զ����ƽ��ֵ�����ͼ�г�����������ϴ��ƫ��ֵ�㡣������Щ������ռ�������٣�����ͼʱ�ɿ��ǽ�����˵���
#ͳ��ƽ��������ȣ��ٷֱȣ�
depth_mean <- round(mean(depth_base$depth), 2)
#�������ε�11����ƽ��������ȵĻ������䡣ʹ��ʣ��Ļ������е�ͳ����Ϣ������ͼչʾ��
depth_base <- subset(depth_base, depth <= 11 * depth_mean)


#ͳ��ƽ�� GC �������ٷֱȣ�
GC_mean <- round(mean(depth_base$GC), 2)
#ggplot2 ɢ��ͼ
depth_GC <- ggplot(depth_base, aes(GC, depth)) +
  geom_point(color = 'gray', alpha = 0.6, pch = 19, size = 0.5) +
  theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_vline(xintercept = GC_mean, color = 'red', lty = 2, lwd = 0.5) + 
  geom_hline(yintercept = depth_mean, color = 'red', lty = 2, lwd = 0.5) +
  labs(x= paste('GC % (Average :', GC_mean, '%)'), y = paste('Depth (Average :', depth_mean, 'X)')) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12))
ggsave('depth_GC_plot.png', depth_GC, width = 5.5, height = 5.5)

#������������װ��������������������GC�������ǲ�����ȣ�������ĳһ����������С�������ͼ��������Ȼ�ܹ�����ɢ����Ҫ�ֲ������򣬵�Ч���������ԡ�������ǿɿ�����ͼ�е�ļ�������ʹ���������ֵ���ɫչʾ����ܶȡ�
#������ɫ�ܶ�
#��ǰ����ͼ�����depth_GC����ʽ�Ļ����ϣ�ʹ��stat_density2d()�������ɢ����ͼ�е��ܶȷֲ�����������ɫ�ܶȱ�ǣ���ʹ��scale_fill_gradientn()ָ����ɫ�ݶȣ�����ɢ��ֲ��ܶ��ɵ͵�����ͼ�б��Ϊ��ͬ��ɫ�������ͼ�����depth_GC_plot.png������
depth_GC <- depth_GC +
  stat_density2d(aes(fill = ..density.., alpha = ..density..), geom = 'tile', contour = FALSE, n = 500, show.legend = FALSE) +
  scale_fill_gradientn(colors = c('transparent', 'gray', 'yellow', 'red'))

ggsave('depth_GC_plot.png', depth_GC, width = 5.5, height = 5.5)
#�����������Depth-GC�ܶȷֲ�ͼ�����Ͼ�������ˡ�ͨ����ͼ���Ի��������Ҫ��֪����Ϣ���ж���װ��������������Լ����������ȡ�




#(2)depth ���ֱ���ܶ�ͼ
#���⣬���ɿ���ͳ�ƻ���GC���������߻�������reads���Ƕȣ����ƶ��ߵķֲ�ֱ��ͼ��������ɢ��ͼ�����࣬�Ը�ֱ��չʾ���ݷֲ������
#geom_histogram()���ڻ���Ƶ��ֱ��ͼ����ʹ��binwidth��������Ϊ100��Ƶ��ͳ�ƴ��񣨼����۲�����������ֲ������������Ƶ������Ϊ100�������ӡ�չʾ��������ʹ��geom_rug()��Ƶ��ֱ��ͼ���·����ӵ�̺ͼ����Ϊչʾÿ�������黬�����в�����ȵ�λ�á������ͼ�����depth_hist_plot.png����
depth_hist <- ggplot(depth_base, aes(depth)) +
  geom_histogram(binwidth = (max(depth_base$depth) - min(depth_base$depth))/100, fill = 'gray', color = 'gray40', size = 0.1) +
  geom_rug(color = 'gray', alpha = 0.6) +
  theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(axis.line = element_line(color = 'black', size = 0.3), axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
  labs(x = '', y = 'Numbers') +
  geom_vline(xintercept = depth_mean, color = 'red', lty = 2, lwd = 0.5)

ggsave('depth_hist_plot.png', depth_hist, width = 5.5, height = 2.5)




#(3)GC ����ֱ���ܶ�ͼ
GC_hist <- ggplot(depth_base, aes(GC)) +
  geom_histogram(binwidth = (max(depth_base$GC) - min(depth_base$GC))/100, fill = 'gray', color = 'gray40', size = 0.1) +
  geom_rug(color = 'gray', alpha = 0.6) +
  theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(axis.line = element_line(color = 'black', size = 0.3), axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
  labs(x = '', y = 'Numbers') +
  geom_vline(xintercept = GC_mean, color = 'red', lty = 2, lwd = 0.5)

ggsave('GC_hist_plot.png', GC_hist, width = 5.5, height = 2.5)




#(4)���ͼ
#����ʹ��grid��������������ɵ�ɢ��ͼ������Ƶ��ֱ��ͼ�����һ�𣬵õ����յĽ����
#��Ϻ��������ͼ�����depth_GC_all_plot.png����
#���� grid ��
library(grid)
#�˴����Ƚ����� depth ���ֱ���ܶ�ͼ��ת 90 ��
depth_hist <- depth_hist + coord_flip()
#grid ���ͼ��
png('depth_GC_all_plot.png', width = 4000, height = 4000, res = 600, units = 'px')
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 3)))
print(depth_GC, vp = viewport(layout.pos.row = 2:3, layout.pos.col = 1:2))
print(GC_hist, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(depth_hist, vp = viewport(layout.pos.row = 2:3, layout.pos.col = 3))
dev.off()
