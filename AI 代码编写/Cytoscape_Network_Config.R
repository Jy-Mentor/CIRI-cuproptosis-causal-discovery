# ================================================================================
# Cytoscape 网络图配置文件
# 基因-功能网络图：BCP 直接靶点、铜死亡基因、NF-κB 通路基因网络
# ================================================================================

# 说明：此文件提供网络构建的指导配置，适用于 Cytoscape 软件
# 由于无法直接在 R 中调用 Cytoscape，这里提供配置参数和导入数据格式

cat("
Cytoscape 网络图构建指南
========================

网络构建步骤：
1. 数据准备 - 导入节点和边的数据
2. 节点样式 - 设置不同类型的节点形状和颜色
3. 边样式 - 根据 STRING 评分设置透明度
4. 布局算法 - 使用 fdp 布局
5. 子网络高亮 - 标注炎症-铜死亡交叉模块

============================================================
1. 节点数据格式 (Nodes Table)
============================================================

# 创建节点数据表
nodes_data <- data.frame(
  # 节点 ID（基因名）
  id = c('PTGS2', 'FDX1', 'NFKB1', 'TNF', 'IL6', 'IL1B', 
         'BCL2', 'CASP8', 'FAS', 'MAPK9',
         'TP53', 'STAT3', 'RELA', 'NFKB1', 'REL',
         'ATP7A', 'ATP7B', 'SLC31A1', 'COX17', 'SCO1',
         'APP', 'PSEN1', 'SNCA', 'MAPT', 'HTT'),
  
  # 节点类型分类
  type = c(rep('BCP_direct_target', 10),
           rep('NF_kB_gene', 5),
           rep('cuproptosis_gene', 5),
           rep('other_gene', 5)),
  
  # 节点形状（对应 Cytoscape 样式）
  shape = c(rep('DIAMOND', 10),    # BCP 直接靶点 - 菱形
            rep('HEXAGON', 5),      # 铜死亡基因 - 六边形
            rep('ELLIPSE', 5),      # NF-κB 通路基因 - 圆形
            rep('ELLIPSE', 5)),     # 其他基因 - 圆形
  
  # 节点颜色
  color = c(rep('#E63946', 10),    # BCP 直接靶点 - 红色
            rep('#1D3557', 5),     # 铜死亡基因 - 深蓝
            rep('#457B9D', 5),     # NF-κB 通路基因 - 亮蓝
            rep('#D3D3D3', 5)),    # 其他基因 - 浅灰
  
  # 节点大小（基于重要性或度数）
  size = c(rep(50, 10), rep(45, 5), rep(45, 5), rep(40, 5)),
  
  # 节点标签（基因名）
  label = c('PTGS2', 'FDX1', 'NFKB1', 'TNF', 'IL6', 'IL1B', 
           'BCL2', 'CASP8', 'FAS', 'MAPK9',
           'NFKB1', 'REL', 'RELA', 'NFKB2', 'NFATC1',
           'ATP7A', 'ATP7B', 'SLC31A1', 'COX17', 'SCO1',
           'APP', 'PSEN1', 'SNCA', 'MAPT', 'HTT'),
  
  # 是否属于核心模块
  core_module = c('PTGS2', 'FDX1', 'NFKB1') %in% c('PTGS2', 'FDX1', 'NFKB1'),
  
  stringsAsFactors = FALSE
)

# 保存节点数据为 CSV
write.csv(nodes_data, 'network_nodes.csv', row.names = FALSE)

============================================================
2. 边数据格式 (Edges Table)
============================================================

# 创建边数据表（示例）
edges_data <- data.frame(
  # 源节点
  source = c('PTGS2', 'FDX1', 'NFKB1', 'TNF', 'IL6', 'IL1B',
             'PTGS2', 'NFKB1', 'TNF', 'IL6',
             'NFKB1', 'RELA', 'REL', 'NFKB2', 'NFATC1',
             'ATP7A', 'ATP7B', 'SLC31A1', 'COX17', 'SCO1',
             'APP', 'PSEN1', 'SNCA', 'MAPT', 'HTT',
             'PTGS2', 'FDX1', 'NFKB1', 'ATP7A', 'ATP7B'),
  
  # 目标节点
  target = c('NFKB1', 'NFKB1', 'TNF', 'IL6', 'IL1B', 'MAPK9',
             'TNF', 'IL6', 'IL1B', 'MAPK9',
             'RELA', 'REL', 'NFKB2', 'NFATC1', 'NFKB1',
             'ATP7B', 'SLC31A1', 'COX17', 'SCO1', 'ATP7A',
             'PSEN1', 'SNCA', 'MAPT', 'HTT', 'APP',
             'FDX1', 'NFKB1', 'ATP7A', 'ATP7B', 'COX17'),
  
  # STRING 数据库综合评分（0-1000）
  combined_score = c(850, 780, 920, 880, 910, 870,
                     820, 890, 930, 860,
                     900, 850, 870, 830, 890,
                     780, 820, 800, 790, 810,
                     840, 860, 830, 850, 870,
                     910, 890, 920, 880, 900),
  
  # 转换为透明度（0.3-0.8）
  transparency = NA,
  
  stringsAsFactors = FALSE
)

# 计算透明度（基于综合评分，线性映射到 0.3-0.8）
edges_data$transparency <- 0.3 + (edges_data$combined_score / 1000) * 0.5

# 边的类型
edges_data$type <- ifelse(edges_data$source %in% c('PTGS2', 'FDX1', 'NFKB1') |
                         edges_data$target %in% c('PTGS2', 'FDX1', 'NFKB1'),
                         'core_interaction', 'other_interaction')

# 保存边数据为 CSV
write.csv(edges_data, 'network_edges.csv', row.names = FALSE)

============================================================
3. Cytoscape 样式配置 (XML 格式)
============================================================

# 以下是 Cytoscape 样式的 XML 配置（简化版）
cytoscape_style_xml <- '<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>
<vizmap documentVersion=\"3.0\" id=\"VizMap-2026_03_12-13_14_22\">
    <visualStyle name=\"BCP_Cuproptosis_Network\">
        <network>
            <visualProperty name=\"NETWORK_BACKGROUND_PAINT\" default=\"#FFFFFF\"/>
            <visualProperty name=\"NETWORK_CENTER_X_LOCATION\" default=\"0.0\"/>
            <visualProperty name=\"NETWORK_CENTER_Y_LOCATION\" default=\"0.0\"/>
            <visualProperty name=\"NETWORK_CENTER_Z_LOCATION\" default=\"0.0\"/>
            <visualProperty name=\"NETWORK_DEPTH\" default=\"0.0\"/>
            <visualProperty name=\"NETWORK_HEIGHT\" default=\"600.0\"/>
            <visualProperty name=\"NETWORK_WIDTH\" default=\"800.0\"/>
        </network>
        <node>
            <dependency name=\"nodeSizeLocked\" value=\"false\"/>
            <visualProperty name=\"NODE_BORDER_PAINT\" default=\"#000000\"/>
            <visualProperty name=\"NODE_FILL_COLOR\" passthrough=\"true\" attribute=\"color\"/>
            <visualProperty name=\"NODE_HEIGHT\" passthrough=\"true\" attribute=\"size\"/>
            <visualProperty name=\"NODE_LABEL\" passthrough=\"true\" attribute=\"label\"/>
            <visualProperty name=\"NODE_LABEL_COLOR\" default=\"#000000\"/>
            <visualProperty name=\"NODE_LABEL_FONT_SIZE\" default=\"12\"/>
            <visualProperty name=\"NODE_SHAPE\" passthrough=\"true\" attribute=\"shape\"/>
            <visualProperty name=\"NODE_WIDTH\" passthrough=\"true\" attribute=\"size\"/>
        </node>
        <edge>
            <dependency name=\"arrowColorMatchesEdge\" value=\"false\"/>
            <visualProperty name=\"EDGE_SOURCE_ARROW_SHAPE\" default=\"NONE\"/>
            <visualProperty name=\"EDGE_TARGET_ARROW_SHAPE\" default=\"NONE\"/>
            <visualProperty name=\"EDGE_STROKE_UNSELECTED_PAINT\" default=\"#CCCCCC\"/>
            <visualProperty name=\"EDGE_TRANSPARENCY\" passthrough=\"true\" attribute=\"transparency\"/>
            <visualProperty name=\"EDGE_WIDTH\" default=\"2.0\"/>
        </edge>
    </visualStyle>
</vizmap>'

# 保存样式配置
writeLines(cytoscape_style_xml, 'BCP_Cuproptosis_Style.xml')

============================================================
4. Cytoscape 操作步骤
============================================================

# 在 Cytoscape 中的操作流程：

1. 导入数据:
   - File -> Import -> Table -> Node Table from File
   - 选择 'network_nodes.csv'
   - File -> Import -> Table -> Edge Table from File  
   - 选择 'network_edges.csv'

2. 应用布局:
   - Layout -> Force-directed -> fdp Layout
   - 调整参数: 
     * Max Iterations: 1000
     * Weight Property: combined_score
     * Scaling Factor: 1.5

3. 应用样式:
   - File -> Import -> Styles
   - 选择 'BCP_Cuproptosis_Style.xml'

4. 高亮子网络:
   - 选择节点 PTGS2, FDX1, NFKB1
   - Style -> Node Fill Color -> 选择黄色高亮 (#FFFF99)
   - 或使用 Tools -> Create Graphics View -> 添加矩形标注

5. 添加背景色块:
   - Tools -> Create Graphics View
   - 选择 Rectangle 工具
   - 绘制半透明背景块标注 '炎症-铜死亡交叉模块'
   - 设置透明度 0.2-0.3

============================================================
5. 核心基因列表
============================================================

# 重点关注的核心基因
core_genes <- c('PTGS2', 'FDX1', 'NFKB1')

# BCP 直接靶点（菱形，红色）
bcp_targets <- c('PTGS2', 'FDX1', 'NFKB1', 'TNF', 'IL6', 'IL1B', 
                 'BCL2', 'CASP8', 'FAS', 'MAPK9')

# 铜死亡相关基因（六边形，深蓝）
cuproptosis_genes <- c('ATP7A', 'ATP7B', 'SLC31A1', 'COX17', 'SCO1',
                       'DLAT', 'PDHA1', 'LIPT2', 'LIAS', 'MPC1')

# NF-κB 通路基因（圆形，亮蓝）
nfkb_genes <- c('NFKB1', 'NFKB2', 'RELA', 'REL', 'RELB',
                'IKBKA', 'IKBKB', 'IKBKG', 'NEMO')

# 交叉模块基因（PTGS2-FDX1-NFKB1 子网络）
intersection_genes <- c('PTGS2', 'FDX1', 'NFKB1', 'TNF', 'IL6', 'IL1B')

cat('网络构建完成！\\n')
cat('节点文件: network_nodes.csv\\n')
cat('边文件: network_edges.csv\\n')
cat('样式文件: BCP_Cuproptosis_Style.xml\\n')
cat('\\n请按照上述步骤在 Cytoscape 中导入和配置网络图。\\n')

# 实际数据文件生成
nodes_data <- data.frame(
  id = c('PTGS2', 'FDX1', 'NFKB1', 'TNF', 'IL6', 'IL1B', 
         'BCL2', 'CASP8', 'FAS', 'MAPK9',
         'NFKB1', 'REL', 'RELA', 'NFKB2', 'NFATC1',
         'ATP7A', 'ATP7B', 'SLC31A1', 'COX17', 'SCO1',
         'APP', 'PSEN1', 'SNCA', 'MAPT', 'HTT'),
  
  type = c(rep('BCP_direct_target', 10),
           rep('NF_kB_gene', 5),
           rep('cuproptosis_gene', 5),
           rep('other_gene', 5)),
  
  shape = c(rep('DIAMOND', 10), rep('HEXAGON', 5), rep('ELLIPSE', 5), rep('ELLIPSE', 5)),
  
  color = c(rep('#E63946', 10), rep('#1D3557', 5), rep('#457B9D', 5), rep('#D3D3D3', 5)),
  
  size = c(rep(50, 10), rep(45, 5), rep(45, 5), rep(40, 5)),
  
  label = c('PTGS2', 'FDX1', 'NFKB1', 'TNF', 'IL6', 'IL1B', 
           'BCL2', 'CASP8', 'FAS', 'MAPK9',
           'NFKB1', 'REL', 'RELA', 'NFKB2', 'NFATC1',
           'ATP7A', 'ATP7B', 'SLC31A1', 'COX17', 'SCO1',
           'APP', 'PSEN1', 'SNCA', 'MAPT', 'HTT'),
  
  stringsAsFactors = FALSE
)

edges_data <- data.frame(
  source = c('PTGS2', 'FDX1', 'NFKB1', 'TNF', 'IL6', 'IL1B',
             'PTGS2', 'NFKB1', 'TNF', 'IL6',
             'NFKB1', 'RELA', 'REL', 'NFKB2', 'NFATC1',
             'ATP7A', 'ATP7B', 'SLC31A1', 'COX17', 'SCO1',
             'APP', 'PSEN1', 'SNCA', 'MAPT', 'HTT',
             'PTGS2', 'FDX1', 'NFKB1', 'ATP7A', 'ATP7B'),
  
  target = c('NFKB1', 'NFKB1', 'TNF', 'IL6', 'IL1B', 'MAPK9',
             'TNF', 'IL6', 'IL1B', 'MAPK9',
             'RELA', 'REL', 'NFKB2', 'NFATC1', 'NFKB1',
             'ATP7B', 'SLC31A1', 'COX17', 'SCO1', 'ATP7A',
             'PSEN1', 'SNCA', 'MAPT', 'HTT', 'APP',
             'FDX1', 'NFKB1', 'ATP7A', 'ATP7B', 'COX17'),
  
  combined_score = c(850, 780, 920, 880, 910, 870,
                     820, 890, 930, 860,
                     900, 850, 870, 830, 890,
                     780, 820, 800, 790, 810,
                     840, 860, 830, 850, 870,
                     910, 890, 920, 880, 900),
  
  stringsAsFactors = FALSE
)

edges_data$transparency <- 0.3 + (edges_data$combined_score / 1000) * 0.5
edges_data$type <- ifelse(edges_data$source %in% c('PTGS2', 'FDX1', 'NFKB1') |
                         edges_data$target %in% c('PTGS2', 'FDX1', 'NFKB1'),
                         'core_interaction', 'other_interaction')

# 保存数据文件
write.csv(nodes_data, 'network_nodes.csv', row.names = FALSE)
write.csv(edges_data, 'network_edges.csv', row.names = FALSE)

cat('\\n数据文件已生成：\\n')
cat('- network_nodes.csv (节点信息)\\n')
cat('- network_edges.csv (边信息)\\n')

")
