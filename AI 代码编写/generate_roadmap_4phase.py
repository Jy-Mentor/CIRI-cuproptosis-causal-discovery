# -*- coding: utf-8 -*-
from PIL import Image, ImageDraw, ImageFont
import os

OUTPUT_PATH = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\AI 代码编写\技术路线图_四阶段.png"

WIDTH, HEIGHT = 2400, 1400
MARGIN = 60

def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

colors = {
    'bg': '#FFFFFF',
    'title': '#1a237e',
    'title_bg': '#ffffff',
    'title_border': '#37474f',
    'phase1_fill': '#e3f2fd',
    'phase1_border': '#1976d2',
    'phase1_label': '#1565c0',
    'phase2_fill': '#fff8e1',
    'phase2_border': '#f9a825',
    'phase2_label': '#f57f17',
    'phase3_fill': '#e8f5e9',
    'phase3_border': '#43a047',
    'phase3_label': '#2e7d32',
    'phase4_fill': '#ffe0b2',
    'phase4_border': '#fb8c00',
    'phase4_label': '#ef6c00',
    'box_fill': '#ffffff',
    'box_border': '#546e7a',
    'box_text': '#263238',
    'label_text': '#ffffff',
    'arrow_main': '#455a64',
    'arrow_dash': '#78909c'
}

def draw_rounded_rect(draw, box, radius, fill, outline, width=2):
    draw.rounded_rectangle([box[0], box[1], box[2], box[3]], radius=radius, fill=fill, outline=outline, width=width)

def create_roadmap():
    img = Image.new('RGB', (WIDTH, HEIGHT), colors['bg'])
    draw = ImageDraw.Draw(img)
    
    try:
        title_font = ImageFont.truetype("msyh.ttc", 32)
        label_font = ImageFont.truetype("msyh.ttc", 20)
        box_title_font = ImageFont.truetype("msyh.ttc", 16)
        box_detail_font = ImageFont.truetype("msyh.ttc", 13)
    except:
        title_font = ImageFont.load_default()
        label_font = ImageFont.load_default()
        box_title_font = ImageFont.load_default()
        box_detail_font = ImageFont.load_default()
    
    title_y = 50
    draw_rounded_rect(draw, [40, title_y - 15, WIDTH - 40, title_y + 45], 10, 
                     hex_to_rgb(colors['title_bg']), hex_to_rgb(colors['title_border']), 3)
    draw.text((WIDTH//2, title_y + 15), "桂艾β-石竹烯调控AGE-RAGE-铜死亡通路干预CIRI技术路线图", 
              font=title_font, fill=hex_to_rgb(colors['title']), anchor="mm")
    
    phase_y = 130
    phase_height = 280
    box_width = 130
    box_height = 85
    box_gap = 20
    phase_gap = 80
    
    phases = [
        {
            'name': '阶段一',
            'subtitle': '生物信息学预测',
            'fill': colors['phase1_fill'],
            'border': colors['phase1_border'],
            'label_color': colors['phase1_label'],
            'items': [
                ('BCP靶点挖掘', 'SwissTargetPrediction'),
                ('CIRI差异基因 ∩\n铜死亡27基因', ''),
                ('PPI网络构建\nGO/KEGG富集', ''),
                ('PC算法共定位\nMR因果推断', ''),
                ('分子对接\nBCP-RAGE/FDX1', '')
            ]
        },
        {
            'name': '阶段二',
            'subtitle': '体内表型验证',
            'fill': colors['phase2_fill'],
            'border': colors['phase2_border'],
            'label_color': colors['phase2_label'],
            'items': [
                ('MCAO/R大鼠\n线栓法建模', ''),
                ('mNSS评分\nTTC梗死体积', ''),
                ('铜稳态检测\n血清/脑ICP-MS', ''),
                ('TTM挽救实验\n铜死亡特异性', '')
            ]
        },
        {
            'name': '阶段三',
            'subtitle': '体外机制阐明',
            'fill': colors['phase3_fill'],
            'border': colors['phase3_border'],
            'label_color': colors['phase3_label'],
            'items': [
                ('SH-SY5Y细胞\nOGD/R造模', ''),
                ('FPS-ZM1阻断\nsiRNA敲低', ''),
                ('线粒体功能\nROS/ATP/JC-1', ''),
                ('铜死亡标志物\nFDX1/DLAT WB', '')
            ]
        },
        {
            'name': '阶段四',
            'subtitle': '数据整合产出',
            'fill': colors['phase4_fill'],
            'border': colors['phase4_border'],
            'label_color': colors['phase4_label'],
            'items': [
                ('多组学整合\n转录组+蛋白组', ''),
                ('机制图绘制\nGraphical Abstract', ''),
                ('论文撰写\n挑战杯申报', '')
            ]
        }
    ]
    
    phase_start_x = 60
    
    for phase_idx, phase in enumerate(phases):
        phase_x = phase_start_x + phase_idx * (5 * box_width + 4 * box_gap + phase_gap)
        
        draw_rounded_rect(draw, [phase_x - 20, phase_y - 40, phase_x + 5 * box_width + 4 * box_gap + 20, phase_y + phase_height], 
                        15, hex_to_rgb(phase['fill']), hex_to_rgb(phase['border']), 3)
        
        label_x = phase_x
        draw_rounded_rect(draw, [label_x, phase_y - 35, label_x + 120, phase_y + 5], 
                        8, hex_to_rgb(phase['label_color']), hex_to_rgb(phase['label_color']), 0)
        draw.text((label_x + 60, phase_y - 15), f"{phase['name']}\n{phase['subtitle']}", 
                  font=label_font, fill=hex_to_rgb(colors['label_text']), anchor="mm", align="center")
        
        for i, item in enumerate(phase['items']):
            box_x = phase_x + i * (box_width + box_gap)
            box_y = phase_y + 20
            
            draw_rounded_rect(draw, [box_x, box_y, box_x + box_width, box_y + box_height], 
                            10, hex_to_rgb(colors['box_fill']), hex_to_rgb(colors['box_border']), 1)
            
            lines = item[0].split('\n')
            text_y = box_y + 15
            for line in lines:
                draw.text((box_x + box_width//2, text_y), line, 
                          font=box_title_font, fill=hex_to_rgb(colors['box_text']), anchor="mm")
                text_y += 18
            
            if i < len(phase['items']) - 1:
                arrow_x = box_x + box_width + box_gap // 2
                draw.line([(arrow_x, box_y + box_height//2), 
                          (arrow_x + box_gap // 2 + 5, box_y + box_height//2)], 
                          fill=hex_to_rgb(colors['arrow_main']), width=2)
                draw.polygon([
                    (arrow_x + box_gap // 2 + 5, box_y + box_height//2),
                    (arrow_x + box_gap // 2 - 5, box_y + box_height//2 - 6),
                    (arrow_x + box_gap // 2 - 5, box_y + box_height//2 + 6)
                ], fill=hex_to_rgb(colors['arrow_main']))
    
    connect_y = phase_y + phase_height + 40
    
    text_y = connect_y
    draw.text((60, text_y), "筛选核心靶点", font=box_title_font, fill=hex_to_rgb(colors['arrow_dash']))
    
    dash_x1 = 60 + 100
    dash_x2 = phase_start_x + (5 * box_width + 4 * box_gap) + phase_gap + 30
    draw.line([(dash_x1, text_y + 25), (dash_x2, text_y + 25)], 
              fill=hex_to_rgb(colors['arrow_dash']), width=2)
    
    p2_items_end_x = phase_start_x + 1 * (5 * box_width + 4 * box_gap + phase_gap) + 4 * (box_width + box_gap) - box_gap
    p2_items_start_x = phase_start_x + 1 * (5 * box_width + 4 * box_gap + phase_gap)
    draw.line([(p2_items_end_x, text_y + 25), (p2_items_end_x, text_y + 60)], 
              fill=hex_to_rgb(colors['arrow_dash']), width=2)
    
    text_y2 = connect_y + 70
    draw.text((60, text_y2), "验证假设一", font=box_title_font, fill=hex_to_rgb(colors['arrow_dash']))
    
    p3_items_end_x = phase_start_x + 2 * (5 * box_width + 4 * box_gap + phase_gap) + 4 * (box_width + box_gap) - box_gap
    draw.line([(dash_x1, text_y2 + 25), (dash_x2, text_y2 + 25)], 
              fill=hex_to_rgb(colors['arrow_dash']), width=2)
    draw.line([(p3_items_end_x, text_y2 + 25), (p3_items_end_x, text_y2 + 60)], 
              fill=hex_to_rgb(colors['arrow_dash']), width=2)
    
    text_y3 = connect_y + 140
    draw.text((60, text_y3), "验证假设二", font=box_title_font, fill=hex_to_rgb(colors['arrow_dash']))
    
    p4_items_end_x = phase_start_x + 3 * (5 * box_width + 4 * box_gap + phase_gap) + 3 * (box_width + box_gap) - box_gap
    draw.line([(dash_x1, text_y3 + 25), (dash_x2, text_y3 + 25)], 
              fill=hex_to_rgb(colors['arrow_dash']), width=2)
    draw.line([(p4_items_end_x, text_y3 + 25), (p4_items_end_x, text_y3 + 60)], 
              fill=hex_to_rgb(colors['arrow_dash']), width=2)
    
    legend_y = HEIGHT - 80
    
    legend_items = [
        ('淡蓝-计算', colors['phase1_fill'], colors['phase1_border']),
        ('淡黄-动物', colors['phase2_fill'], colors['phase2_border']),
        ('淡绿-细胞', colors['phase3_fill'], colors['phase3_border']),
        ('淡橙-产出', colors['phase4_fill'], colors['phase4_border'])
    ]
    
    legend_x = WIDTH // 2 - 400
    for i, (label, fill, border) in enumerate(legend_items):
        draw_rounded_rect(draw, [legend_x + i * 120, legend_y, legend_x + i * 120 + 100, legend_y + 30], 
                        5, hex_to_rgb(fill), hex_to_rgb(border), 1)
        draw.text((legend_x + i * 120 + 50, legend_y + 15), label, 
                  font=box_detail_font, fill=hex_to_rgb(colors['box_text']), anchor="mm")
    
    img.save(OUTPUT_PATH)
    print(f"四阶段技术路线图已保存至: {OUTPUT_PATH}")
    return img

if __name__ == "__main__":
    create_roadmap()
