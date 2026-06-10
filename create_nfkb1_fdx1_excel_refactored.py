#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NFKB1到FDX1最短路径分析Excel报告生成器（重构版）
使用ExcelReportBuilder工具类，消除重复代码

原文件: create_nfkb1_fdx1_excel.py
重构日期: 2025-01-24
"""

import os
import sys

# 添加utils目录到路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'utils'))

from excel_report_builder import ExcelReportBuilder


def create_nfkb1_fdx1_excel():
    """
    创建NFKB1-FDX1通路分析Excel报告（重构版）
    使用ExcelReportBuilder工具类简化代码
    """
    # 配置路径
    work_dir = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
    result_dir = os.path.join(work_dir, "String_Network_Systematic_Analysis")
    output_file = os.path.join(work_dir, "NFKB1_FDX1_Pathway_Analysis.xlsx")
    
    # 创建报告构建器
    builder = ExcelReportBuilder(output_file, theme='default')
    
    # ========== Sheet 1: 分析概述 ==========
    ws_summary = builder.create_sheet("分析概述")
    
    # 使用add_info_section快速添加信息区块
    summary_data = [
        ["分析背景", ""],
        ["分析网络", "g_k3 (K>=3核心子网络)"],
        ["网络节点数", 124],
        ["网络边数", 1924],
        ["", ""],
        ["最短路径统计", ""],
        ["源节点", "NFKB1 (铜死亡通路关键调控因子)"],
        ["目标节点", "FDX1 (铜死亡核心基因)"],
        ["最短路径数", 8],
        ["最短路径长度", "2跳 (3个节点)"],
        ["", ""],
        ["桥接节点分析", ""],
        ["桥接节点总数", 2],
        ["HSPA5通路", "4条路径 (50%)"],
        ["HMOX1通路", "4条路径 (50%)"],
        ["", ""],
        ["关键发现", ""],
        ["通路1", "NFKB1 → HSPA5 → FDX1"],
        ["通路2", "NFKB1 → HMOX1 → FDX1"],
        ["生物学意义", "NFKB1可能通过HSPA5或HMOX1间接调控FDX1"]
    ]
    
    # 使用add_table简化数据写入
    builder.add_title(ws_summary, "NFKB1到FDX1最短路径分析报告", row=1, col=1)
    builder.add_table(ws_summary, summary_data, start_row=3, start_col=1)
    
    # 为子标题行应用样式
    subtitle_rows = [3, 9, 15, 21]  # 子标题所在的行
    for row in subtitle_rows:
        ws_summary.cell(row=row, column=1).font = builder.subtitle_font
        ws_summary.cell(row=row, column=1).fill = builder.subtitle_fill
    
    # 设置列宽
    builder.set_column_widths(ws_summary, {
        'A': 20,
        'B': 45
    })
    
    # ========== Sheet 2: 8条最短路径 ==========
    ws_paths = builder.create_sheet("8条最短路径")
    builder.add_title(ws_paths, "NFKB1到FDX1的8条最短路径", row=1, col=1)
    
    # 简化显示：去重后的路径
    unique_paths = [
        ["路径1-4", "NFKB1 → HSPA5 → FDX1", "HSPA5通路"],
        ["路径5-8", "NFKB1 → HMOX1 → FDX1", "HMOX1通路"]
    ]
    
    builder.add_subtitle(ws_paths, "去重后的通路（共2条）", row=3)
    headers = ['路径ID', '通路', "通路类型"]
    path_data = [headers] + unique_paths
    builder.add_table(ws_paths, path_data, start_row=4, has_header=True)
    
    # 详细路径说明（静态数据）
    builder.add_subtitle(ws_paths, "详细路径信息", row=8)
    detailed_info = [
        ["路径1", "NFKB1 → HSPA5 → FDX1", "HSPA5作为桥接节点"],
        ["路径2", "NFKB1 → HMOX1 → FDX1", "HMOX1作为桥接节点"]
    ]
    detail_headers = ['路径ID', '通路', '说明']
    detail_data = [detail_headers] + detailed_info
    builder.add_table(ws_paths, detail_data, start_row=9, has_header=True)
    
    builder.set_column_widths(ws_paths, {
        'A': 15,
        'B': 30,
        'C': 25
    })
    
    # ========== Sheet 3: 桥接节点分析 ==========
    ws_bridge = builder.create_sheet("桥接节点分析")
    builder.add_title(ws_bridge, "桥接节点频率分析", row=1, col=1)
    
    # 桥接节点数据
    bridge_data = [
        ["桥接节点", "出现次数", "占比"],
        ["HSPA5", 4, "50%"],
        ["HMOX1", 4, "50%"]
    ]
    builder.add_table(ws_bridge, bridge_data, start_row=2, has_header=True)
    
    # 桥接节点说明
    builder.add_subtitle(ws_bridge, "桥接节点功能说明", row=6)
    bridge_info = [
        ["基因", "别名", "功能"],
        ["HSPA5", "GRP78/BiP", "内质网应激分子伴侣，参与蛋白折叠和应激反应"],
        ["HMOX1", "血红素加氧酶-1", "抗氧化应激关键酶，参与氧化还原调控"]
    ]
    builder.add_table(ws_bridge, bridge_info, start_row=7, has_header=True)
    
    builder.set_column_widths(ws_bridge, {
        'A': 15,
        'B': 15,
        'C': 50
    })
    
    # ========== Sheet 4: 生物学意义解读 ==========
    ws_bio = builder.create_sheet("生物学意义解读")
    builder.add_title(ws_bio, "NFKB1到FDX1通路的生物学意义", row=1, col=1)
    
    bio_content = [
        ["1. 通路概述"],
        ["NFKB1（核因子κB1）是炎症反应和细胞存活的关键转录因子。"],
        ["FDX1（铁氧还蛋白1）是铜死亡（Cuproptosis）的核心调控基因。"],
        ["本分析发现NFKB1到FDX1的最短距离仅为2跳，表明两者存在紧密的功能联系。"],
        [""],
        ["2. 两条调控通路"],
        ["通路1: NFKB1 → HSPA5 → FDX1"],
        ["  • HSPA5（GRP78/BiP）是内质网应激的关键分子伴侣"],
        ["  • 铜死亡诱导内质网应激，HSPA5可能作为信号中介"],
        ["  • NFKB1可能通过调控HSPA5影响FDX1介导的铜死亡"],
        [""],
        ["通路2: NFKB1 → HMOX1 → FDX1"],
        ["  • HMOX1（血红素加氧酶-1）是重要的抗氧化酶"],
        ["  • 铜死亡涉及氧化应激，HMOX1可能参与氧化还原调控"],
        ["  • NFKB1可能通过HMOX1调控FDX1相关的氧化应激反应"],
        [""],
        ["3. 研究意义"],
        ["• 这两条通路为理解BCP（β-石竹烯）调控铜死亡的分子机制提供了线索"],
        ["• BCP可能通过抑制NFKB1，进而影响HSPA5/HMOX1，最终调控FDX1介导的铜死亡"],
        ["• HSPA5和HMOX1可作为潜在的药物靶点进行验证"],
        [""],
        ["4. 后续实验建议"],
        ["• 验证NFKB1对HSPA5和HMOX1的转录调控"],
        ["• 检测HSPA5和HMOX1对FDX1表达的影响"],
        ["• 在铜死亡模型中验证这两条通路的功能重要性"],
        ["• 评估BCP对这些节点表达的影响"]
    ]
    
    builder.add_table(ws_bio, bio_content, start_row=3, start_col=1)
    
    # 为子标题行应用样式
    subtitle_texts = ["1. 通路概述", "2. 两条调控通路", "3. 研究意义", "4. 后续实验建议"]
    for row_idx, row_data in enumerate(bio_content, 3):
        if row_data[0] in subtitle_texts:
            ws_bio.cell(row=row_idx, column=1).font = builder.subtitle_font
            ws_bio.cell(row=row_idx, column=1).fill = builder.subtitle_fill
    
    builder.set_column_width(ws_bio, 'A', 80)
    
    # 保存文件
    builder.save()
    
    print(f"✅ NFKB1-FDX1通路分析汇总Excel已生成: {output_file}")
    print("包含以下sheet:")
    for sheet in builder.get_sheetnames():
        print(f"  - {sheet}")


if __name__ == "__main__":
    create_nfkb1_fdx1_excel()
