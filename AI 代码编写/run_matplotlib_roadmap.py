import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(16, 10))
ax.set_xlim(0, 16)
ax.set_ylim(0, 10)
ax.axis('off')

# 配色方案 (JMIR Pastel)
colors = {
    'p1_bg': '#E3F2FD',      # 淡蓝背景
    'p1_border': '#1976D2',  # 深蓝边框
    'p2_bg': '#FFF8E1',      # 淡黄背景  
    'p2_border': '#F9A825',  # 深黄边框
    'p3_bg': '#E8F5E9',      # 淡绿背景
    'p3_border': '#43A047',  # 深绿边框
    'p4_bg': '#FFE0B2',      # 淡橙背景
    'p4_border': '#FB8C00',  # 深橙边框
    'text': '#263238',
    'arrow': '#78909C'
}

def draw_phase_box(ax, y_bottom, y_top, color_bg, color_border, label_text, label_color):
    """绘制阶段背景框和左侧标签"""
    # 背景大框
    bg = FancyBboxPatch((1.0, y_bottom), 14.0, y_top-y_bottom,
                        boxstyle="round,pad=0.02,rounding_size=0.3",
                        facecolor=color_bg,
                        edgecolor=color_border,
                        linewidth=3)
    ax.add_patch(bg)
    
    # 左侧标签（竖排文字）
    label_bg = FancyBboxPatch((0.3, y_bottom+0.3), 0.6, y_top-y_bottom-0.6,
                              boxstyle="round,pad=0.02,rounding_size=0.1",
                              facecolor=color_border,
                              edgecolor='none')
    ax.add_patch(label_bg)
    
    # 竖排文字
    ax.text(0.6, (y_bottom+y_top)/2, label_text,
            ha='center', va='center', fontsize=11, weight='bold', color='white',
            rotation=90)

def draw_process_box(ax, x, y, width, height, text, border_color):
    """绘制流程步骤框"""
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.3",
                         facecolor='white',
                         edgecolor=border_color,
                         linewidth=2)
    ax.add_patch(box)
    
    ax.text(x + width/2, y + height/2, text,
            ha='center', va='center', fontsize=10, color=colors['text'],
            linespacing=1.3, weight='medium')

def draw_arrow_with_text(ax, x1, y1, x2, y2, arrow_text=None):
    """绘制带标注的箭头"""
    # 主箭头
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=2))
    
    # 文字标注（在箭头上方）
    if arrow_text:
        mid_x, mid_y = (x1+x2)/2, (y1+y2)/2
        ax.text(mid_x, mid_y+0.15, arrow_text, ha='center', va='bottom',
                fontsize=9, color='#455A64', style='italic')

# ==================== 标题 ====================
ax.text(8, 9.5, '桂艾β-石竹烯调控AGE-RAGE-铜死亡通路干预CIRI技术路线图',
        ha='center', va='center', fontsize=18, weight='bold', color='#1A237E')

# ==================== 阶段一：生物信息学预测 (y: 7.2-8.8) ====================
draw_phase_box(ax, 7.2, 8.8, colors['p1_bg'], colors['p1_border'], 
               '阶段一\n生物信息学预测', colors['p1_border'])

# 横向流程框
p1_boxes = [
    (2.0, 8.0, 'BCP靶点挖掘\nSwissTargetPrediction'),
    (4.5, 8.0, 'CIRI差异基因 ∩\n铜死亡27基因'),
    (7.0, 8.0, 'PPI网络构建\nGO/KEGG富集'),
    (9.5, 8.0, 'PC算法共定位\nMR因果推断'),
    (12.0, 8.0, '分子对接验证\nBCP-RAGE/FDX1')
]

for i, (x, y, text) in enumerate(p1_boxes):
    draw_process_box(ax, x, y, 2.2, 0.6, text, colors['p1_border'])
    if i < len(p1_boxes)-1:
        ax.annotate('', xy=(x+2.4, y+0.3), xytext=(x+2.2, y+0.3),
                   arrowprops=dict(arrowstyle='->', color=colors['p1_border'], lw=2))

# ==================== 阶段二：体内表型验证 (y: 5.2-6.8) ====================
draw_phase_box(ax, 5.2, 6.8, colors['p2_bg'], colors['p2_border'], 
               '阶段二\n体内表型验证', colors['p2_border'])

p2_boxes = [
    (2.5, 6.0, 'MCAO/R大鼠\n线栓法建模'),
    (5.2, 6.0, 'mNSS评分\nTTC梗死体积'),
    (7.9, 6.0, '铜稳态检测\n血清/脑ICP-MS'),
    (10.6, 6.0, 'TTM挽救实验\n铜死亡特异性验证')
]

for i, (x, y, text) in enumerate(p2_boxes):
    draw_process_box(ax, x, y, 2.4, 0.6, text, colors['p2_border'])
    if i < len(p2_boxes)-1:
        ax.annotate('', xy=(x+2.6, y+0.3), xytext=(x+2.4, y+0.3),
                   arrowprops=dict(arrowstyle='->', color=colors['p2_border'], lw=2))

# 阶段一到阶段二连接
draw_arrow_with_text(ax, 8.0, 7.2, 8.0, 6.8, '筛选核心靶点(RAGE/FDX1)')

# ==================== 阶段三：体外机制阐明 (y: 3.2-4.8) ====================
draw_phase_box(ax, 3.2, 4.8, colors['p3_bg'], colors['p3_border'], 
               '阶段三\n体外机制阐明', colors['p3_border'])

p3_boxes = [
    (2.5, 4.0, 'SH-SY5Y细胞\nOGD/R造模'),
    (5.2, 4.0, 'FPS-ZM1阻断\nsiRNA敲低'),
    (7.9, 4.0, '线粒体功能\nROS/ATP/JC-1'),
    (10.6, 4.0, '铜死亡标志物\nFDX1/DLAT WB')
]

for i, (x, y, text) in enumerate(p3_boxes):
    draw_process_box(ax, x, y, 2.4, 0.6, text, colors['p3_border'])
    if i < len(p3_boxes)-1:
        ax.annotate('', xy=(x+2.6, y+0.3), xytext=(x+2.4, y+0.3),
                   arrowprops=dict(arrowstyle='->', color=colors['p3_border'], lw=2))

# 阶段二到阶段三连接
draw_arrow_with_text(ax, 8.0, 5.2, 8.0, 4.8, '验证假设一(体内→体外)')

# ==================== 阶段四：数据整合产出 (y: 1.2-2.8) ====================
draw_phase_box(ax, 1.2, 2.8, colors['p4_bg'], colors['p4_border'], 
               '阶段四\n数据整合产出', colors['p4_border'])

p4_boxes = [
    (3.5, 2.0, '多组学整合分析\n(转录组+蛋白组)'),
    (6.8, 2.0, '机制图绘制\nGraphical Abstract'),
    (10.1, 2.0, '论文撰写与发表\n挑战杯竞赛申报')
]

for i, (x, y, text) in enumerate(p4_boxes):
    draw_process_box(ax, x, y, 2.8, 0.6, text, colors['p4_border'])
    if i < len(p4_boxes)-1:
        ax.annotate('', xy=(x+3.0, y+0.3), xytext=(x+2.8, y+0.3),
                   arrowprops=dict(arrowstyle='->', color=colors['p4_border'], lw=2))

# 阶段三到阶段四连接
draw_arrow_with_text(ax, 8.0, 3.2, 8.0, 2.8, '验证假设二(机制阐明)')

OUTPUT_DIR = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\AI 代码编写"

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}\\tech_roadmap_v2.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f'{OUTPUT_DIR}\\tech_roadmap_v2.svg', format='svg', bbox_inches='tight', facecolor='white')
plt.close()
print(f"已生成 {OUTPUT_DIR}\\tech_roadmap_v2.png")
print(f"已生成 {OUTPUT_DIR}\\tech_roadmap_v2.svg")