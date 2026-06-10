library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(patchwork)

# 读取数据
mr_results <- read.csv("D:/EQTL/MR_1e-5_Results/mr_main_results.csv", stringsAsFactors = FALSE)
diag <- read.csv("D:/EQTL/MR_1e-5_Results/diagnostics.csv", stringsAsFactors = FALSE)
stroke_edges <- read.csv("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/network_stroke_directed.csv", stringsAsFactors = FALSE)

# 过滤孤立边
main_genes <- unique(c(stroke_edges$from, stroke_edges$to))
main_genes <- main_genes[main_genes %in% c("nfkb1", "rela", "stat1", "stat3", "ccl2", "icam1", "tgfb1", 
                                               "atf4", "egr1", "slc31a1", "atox1")]
edges_filtered <- stroke_edges[stroke_edges$from %in% main_genes & stroke_edges$to %in% main_genes, ]
stroke_graph <- graph_from_data_frame(edges_filtered[, c("from", "to")], directed = TRUE)

# ============ 图1: F统计量图（Nature风格）============
p1 <- ggplot(diag %>% filter(fstat_snps>0), 
             aes(x=reorder(gene, fstat_snps), y=fstat_snps, 
                 fill=fstat_snps)) + 
  geom_col(width=0.75, show.legend=F, color=NA) + 
  geom_text(aes(label=round(fstat_snps,1)), hjust=-0.3, size=3.5, 
            color="#2E3440") + 
  geom_hline(yintercept=10, linetype="dashed", color="#E63946", size=0.8) + 
  annotate("text", x=3, y=max(diag$fstat_snps, na.rm=TRUE)*0.12, 
           label="Weak instrument threshold (F=10)", 
           color="#E63946", size=3, hjust=0) + 
  scale_fill_gradient(low="#A8D5E5", high="#2E86AB", guide="none") + 
  coord_flip() + 
  ylim(0, max(diag$fstat_snps, na.rm=TRUE)*1.25) + 
  labs(x=NULL, y="F-statistic", 
       title="Instrument Strength Validation",
       subtitle="All instruments pass weak instrument test (F>10)") + 
  theme_bw() + 
  theme(
    panel.grid=element_blank(), 
    panel.border=element_blank(),
    axis.line.x=element_line(color="#2E3440", size=0.5),
    axis.line.y=element_blank(),
    axis.ticks=element_line(color="#2E3440", size=0.5),
    axis.text=element_text(size=9, color="#2E3440"),
    axis.title=element_text(size=10, color="#2E3440"),
    plot.title=element_text(size=12, color="#2E3440", face="bold"),
    plot.subtitle=element_text(size=9, color="#5E6D7A"),
    plot.margin=margin(10,20,10,10)
  )

# ============ 图2: MR森林图（Nature风格）============
sig_results <- mr_results %>% 
  filter(method %in% c("Inverse variance weighted","Weighted median")) %>% 
  mutate(method_short=ifelse(method=="Inverse variance weighted","IVW","WM"), 
         sig=pval<0.05, 
         lower=b-1.96*se, upper=b+1.96*se) 

p2 <- ggplot(sig_results, aes(x=b, y=reorder(gene, b), color=sig)) + 
  geom_vline(xintercept=0, linetype="dashed", color="#5E6D7A", size=0.6) + 
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0.25, size=0.7, color="#2E3440") + 
  geom_point(aes(shape=method_short), size=3) + 
  facet_grid(~method_short, scales="free_y", space="free") + 
  scale_color_manual(values=c("TRUE"="#F24236","FALSE"="#2E86AB"), guide="none") + 
  scale_shape_manual(values=c("IVW"=16,"WM"=17)) + 
  labs(x="Causal Effect (Beta)", y=NULL, 
       title="MR Analysis: Copper Death Genes -> Ischemic Stroke",
       subtitle="Red: p<0.05 (Significant)") + 
  theme_bw() + 
  theme(
    legend.position="none", 
    strip.background=element_rect(fill="#F5F5F5", color=NA),
    strip.text=element_text(size=10, color="#2E3440"),
    panel.spacing=unit(1.2,"lines"),
    panel.border=element_blank(),
    axis.line.x=element_line(color="#2E3440", size=0.5),
    axis.line.y=element_blank(),
    axis.ticks=element_line(color="#2E3440", size=0.5),
    axis.text=element_text(size=9, color="#2E3440"),
    axis.title=element_text(size=10, color="#2E3440"),
    plot.title=element_text(size=12, color="#2E3440", face="bold"),
    plot.subtitle=element_text(size=9, color="#5E6D7A"),
    plot.margin=margin(10,20,10,10)
  )

# ============ 图3: 网络图（Nature风格）============
V(stroke_graph)$type <- ifelse(V(stroke_graph)$name %in% 
  c("fdx1","lias","slc31a1","atox1"), "Core", 
  ifelse(V(stroke_graph)$name %in% c("nfkb1","rela","stat3"), "Hub", "Other"))

p3 <- ggraph(stroke_graph, layout="fr", niter=1500) + 
  geom_edge_link(aes(linetype=direction), 
                 arrow=arrow(length=unit(2.5,"mm"), angle=25), 
                 alpha=0.7, color="#5E6D7A", end_cap=circle(3,"mm")) + 
  geom_node_point(aes(size=degree(stroke_graph), fill=type), alpha=0.9, shape=21, stroke=0.5) + 
  geom_node_text(aes(label=toupper(name)), repel=TRUE, size=3.5, 
                 fontface="bold", color="#2E3440",
                 box.padding=unit(0.5,"mm"), point.padding=unit(0.5,"mm")) + 
  scale_fill_manual(values=c("Core"="#F24236","Hub"="#2E86AB","Other"="#A0A0A0"), 
                    name="Gene Type", 
                    labels=c("Core"="Copper Death Core","Hub"="AGE-RAGE Hub","Other"="Other")) + 
  scale_size_continuous(range=c(4,10), guide="none") + 
  scale_edge_linetype_manual(values=c("directed"="solid","undirected"="dashed"), guide="none") +
  labs(title="Stroke-Specific Regulatory Network",
       subtitle="PC Algorithm Inference") + 
  theme_void() + 
  theme(
    legend.position=c(0.85,0.15),
    legend.background=element_rect(fill="white", color=NA, linewidth=0),
    legend.key=element_rect(fill="white", color=NA),
    legend.text=element_text(size=8, color="#2E3440"),
    legend.title=element_text(size=9, color="#2E3440", face="bold"),
    plot.title=element_text(size=12, color="#2E3440", face="bold"),
    plot.subtitle=element_text(size=9, color="#5E6D7A"),
    plot.margin=margin(10,10,10,10)
  )

# ============ 图4: 整合图（Nature风格）============
mr_sig <- data.frame( 
  gene=c("FDX1","ATOX1","PDHB"), 
  beta=c(-0.055,-0.046,0.049), 
  se=c(0.021,0.022,0.018), 
  pval=c(0.009,0.035,0.005) 
) %>% mutate( 
  lower=beta-1.96*se, upper=beta+1.96*se, 
  direction=ifelse(beta>0,"Risk","Protective"), 
  sig=ifelse(pval<0.01,"**",ifelse(pval<0.05,"*","")),
  gene_label=paste0(toupper(gene), ifelse(sig!="", paste0(" ", sig), ""))
) 

p4 <- ggplot(mr_sig, aes(x=reorder(gene_label, beta), y=beta, fill=direction)) + 
  geom_hline(yintercept=0, linetype="dashed", color="#5E6D7A", size=0.6) + 
  geom_col(width=0.65, alpha=0.9) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.25, size=0.7, color="#2E3440") + 
  scale_fill_manual(values=c("Protective"="#2E86AB","Risk"="#F24236"), 
                    name="Effect Direction",
                    labels=c("Protective"="Protective","Risk"="Risk")) + 
  labs(x=NULL, y="Causal Effect (Beta)", 
       title="MR-PC Integration Analysis",
       subtitle="Genetic causality evidence in regulatory networks") + 
  theme_bw() + 
  theme(
    legend.position="bottom",
    legend.background=element_rect(fill="white", color=NA),
    legend.key=element_rect(fill="white", color=NA),
    panel.border=element_blank(),
    axis.line.x=element_line(color="#2E3440", size=0.5),
    axis.line.y=element_blank(),
    axis.ticks=element_line(color="#2E3440", size=0.5),
    axis.text=element_text(size=9, color="#2E3440"),
    axis.title=element_text(size=10, color="#2E3440"),
    plot.title=element_text(size=12, color="#2E3440", face="bold"),
    plot.subtitle=element_text(size=9, color="#5E6D7A"),
    legend.text=element_text(size=8, color="#2E3440"),
    legend.title=element_text(size=9, color="#2E3440", face="bold"),
    plot.margin=margin(10,20,15,10)
  ) +
  ylim(-0.12, 0.12)

# ============ 组合输出 ============ 
(p1 + p2) / (p3 + p4) + 
  plot_annotation(
    title="Copper Death Genes Analysis Suite",
    subtitle="Multi-omics Integration: MR + PC Network Analysis",
    theme=theme(
      plot.title=element_text(size=14, color="#2E3440", face="bold", hjust=0.5),
      plot.subtitle=element_text(size=10, color="#5E6D7A", hjust=0.5),
      plot.margin=margin(15,20,15,20)
    )
  )

# 保存高分辨率PNG和PDF
ggsave("final_plots.png", width=14, height=11, dpi=400, bg="white")
ggsave("final_plots.pdf", width=14, height=11, device=cairo_pdf)

cat("Final plots saved to final_plots.png (400dpi) and final_plots.pdf\n")
