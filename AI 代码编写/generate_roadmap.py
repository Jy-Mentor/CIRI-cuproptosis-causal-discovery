# -*- coding: utf-8 -*-
from PIL import Image, ImageDraw, ImageFont
import os

OUTPUT_PATH = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\AI 代码编写\技术路线图.png"

WIDTH = 1400
HEIGHT = 2200
MARGIN = 100
BOX_WIDTH = 1100
BOX_HEIGHT = 200

def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def hex_to_rgba(hex_color, alpha=255):
    rgb = hex_to_rgb(hex_color)
    return rgb + (alpha,)

def lighten_color(hex_color, factor=0.3):
    rgb = hex_to_rgb(hex_color)
    return tuple(min(255, int(c + (255 - c) * factor)) for c in rgb)

colors = {
    'bg': '#FFFFFF',
    'bg_gray': '#F7FAFC',
    'title': '#1a365d',
    'phase1': '#1a365d',
    'phase2': '#1e40af',
    'phase3': '#2563eb',
    'phase4': '#3b82f6',
    'accent': '#0d9488',
    'text': '#FFFFFF',
    'subtext': '#2d3748',
    'subtext_light': '#718096',
    'line': '#4a5568',
    'subline': '#a0aec0',
    'border': '#e2e8f0'
}

phases = [
    {
        'title': '第一阶段：生物信息学预测',
        'title_en': 'Phase 1: Bioinformatics Prediction',
        'content': [
            '• 分子对接预测β-石竹烯与RAGE蛋白结合位点',
            '• GEO数据库挖掘CIRI差异表达基因',
            '• 构建PPI网络筛选核心靶点',
            '• 通路富集分析AGE-RAGE/铜死亡轴'
        ],
        'color': colors['phase1']
    },
    {
        'title': '第二阶段：体内模型验证',
        'title_en': 'Phase 2: In Vivo Validation',
        'content': [
            '• 构建大鼠CIRI模型',
            '• 桂艾β-石竹烯干预处理',
            '• 神经功能评分与组织病理学检测',
            '• Western Blot检测RAGE/铜死亡通路蛋白'
        ],
        'color': colors['phase2']
    },
    {
        'title': '第三阶段：体外机制阐明',
        'title_en': 'Phase 3: In Vitro Mechanism',
        'content': [
            '• 细胞实验：OGD/R模型建立',
            '• siRNA干扰与过表达验证',
            '• ROS/线粒体功能检测',
            '• 铜死亡关键指标检测（铜离子、FDX1等）'
        ],
        'color': colors['phase3']
    },
    {
        'title': '第四阶段：数据整合输出',
        'title_en': 'Phase 4: Data Integration',
        'content': [
            '• 多组学数据整合分析',
            '• 构建调控网络图谱',
            '• 撰写研究报告与论文',
            '• 申报材料整理与成果转化'
        ],
        'color': colors['phase4']
    }
]

def create_roadmap():
    img = Image.new('RGBA', (WIDTH, HEIGHT), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)
    
    for y in range(250):
        factor = y / 250
        r = int(247 + (255 - 247) * factor)
        g = int(250 + (255 - 250) * factor)
        b = int(252 + (255 - 252) * factor)
        draw.line([(0, y), (WIDTH, y)], fill=(r, g, b, 255))
    
    for x in range(0, WIDTH, 40):
        draw.line([(x, 0), (x, 10)], fill=hex_to_rgb('#cbd5e0') + (255,))
    
    try:
        title_font = ImageFont.truetype("msyh.ttc", 48)
        phase_font = ImageFont.truetype("msyh.ttc", 32)
        content_font = ImageFont.truetype("msyh.ttc", 20)
        en_font = ImageFont.truetype("msyh.ttc", 16)
        small_font = ImageFont.truetype("msyh.ttc", 18)
    except:
        title_font = ImageFont.load_default()
        phase_font = ImageFont.load_default()
        content_font = ImageFont.load_default()
        en_font = ImageFont.load_default()
        small_font = ImageFont.load_default()
    
    title_y = 80
    draw.text((WIDTH//2, title_y), "桂艾β-石竹烯调控AGE-RAGE/铜死亡轴", 
              font=title_font, fill=hex_to_rgb(colors['title']) + (255,), anchor="mm")
    
    subtitle_y = title_y + 60
    draw.text((WIDTH//2, subtitle_y), "干预CIRI的机制研究", 
              font=phase_font, fill=hex_to_rgb(colors['phase2']) + (255,), anchor="mm")
    
    y = 280
    for i, phase in enumerate(phases):
        box_y = y
        
        shadow = Image.new('RGBA', (BOX_WIDTH + 20, BOX_HEIGHT + 20), (0, 0, 0, 0))
        shadow_draw = ImageDraw.Draw(shadow)
        shadow_draw.rounded_rectangle([0, 0, BOX_WIDTH + 18, BOX_HEIGHT + 18], radius=25, 
                                       fill=(0, 0, 0, 25))
        img.paste(shadow, (MARGIN - 10, box_y - 10), shadow)
        
        draw.rounded_rectangle(
            [MARGIN, box_y, MARGIN + BOX_WIDTH, box_y + BOX_HEIGHT],
            radius=20, fill=hex_to_rgb(phase['color']) + (255,), 
            outline=lighten_color(phase['color'], 0.2) + (255,), width=2
        )
        
        phase_num = f"0{i+1}"
        draw.text((MARGIN + 30, box_y + 25), phase_num, 
                  font=phase_font, fill=(255, 255, 255, 128))
        
        draw.text((MARGIN + 110, box_y + 30), phase['title'], 
                  font=phase_font, fill=hex_to_rgb(colors['text']) + (255,))
        
        draw.text((MARGIN + 110, box_y + 75), phase['title_en'], 
                  font=en_font, fill=(255, 255, 255, 180))
        
        content_x = MARGIN + 40
        content_y = box_y + 115
        for line in phase['content']:
            draw.text((content_x, content_y), line, 
                      font=content_font, fill=hex_to_rgb(colors['text']) + (255,))
            content_y += 26
        
        if i < len(phases) - 1:
            arrow_y = box_y + BOX_HEIGHT + 25
            
            center_x = WIDTH // 2
            draw.line([(center_x, arrow_y), (center_x, arrow_y + 30)], 
                     fill=hex_to_rgb(colors['accent']) + (255,), width=4)
            
            draw.polygon([
                (center_x - 15, arrow_y + 30),
                (center_x + 15, arrow_y + 30),
                (center_x, arrow_y + 55)
            ], fill=hex_to_rgb(colors['accent']) + (255,))
            
            y = box_y + BOX_HEIGHT + 90
        else:
            y = box_y + BOX_HEIGHT + 60
    
    total_height = y + 80
    
    final_img = Image.new('RGBA', (WIDTH, total_height), (255, 255, 255, 255))
    final_img.paste(img, (0, 0), img)
    draw = ImageDraw.Draw(final_img)
    
    center_x = WIDTH // 2
    
    footer_y = total_height - 100
    draw.line([(MARGIN, footer_y - 20), (WIDTH - MARGIN, footer_y - 20)], 
              fill=hex_to_rgb(colors['border']) + (255,), width=1)
    
    draw.text((WIDTH//2, footer_y + 15), "研究思路", 
              font=small_font, fill=hex_to_rgb(colors['subtext']) + (255,), anchor="mm")
    
    flow_y = footer_y + 50
    steps = ["生物信息学预测", "体内模型验证", "体外机制阐明", "数据整合输出"]
    step_colors = [colors['phase1'], colors['phase2'], colors['phase3'], colors['phase4']]
    
    total_width = 0
    for step in steps:
        bbox = draw.textbbox((0, 0), step, font=small_font)
        total_width += bbox[2] - bbox[0]
    total_width += 60 * 3
    
    start_x = (WIDTH - total_width) // 2
    
    x_pos = start_x
    for j, (step, color) in enumerate(zip(steps, step_colors)):
        draw.text((x_pos, flow_y), step, font=small_font, fill=hex_to_rgb(color) + (255,))
        
        bbox = draw.textbbox((0, 0), step, font=small_font)
        text_width = bbox[2] - bbox[0]
        
        x_pos += text_width + 60
        
        if j < len(steps) - 1:
            draw.line([(x_pos - 30, flow_y + 10), (x_pos, flow_y + 10)], 
                     fill=hex_to_rgb(colors['subline']) + (255,), width=2)
            draw.polygon([
                (x_pos, flow_y + 10),
                (x_pos - 8, flow_y + 6),
                (x_pos - 8, flow_y + 14)
            ], fill=hex_to_rgb(colors['subline']) + (255,))
    
    rgb_img = Image.new('RGB', final_img.size, (255, 255, 255))
    rgb_img.paste(final_img, mask=final_img.split()[3])
    
    rgb_img.save(OUTPUT_PATH)
    print(f"技术路线图已保存至: {OUTPUT_PATH}")
    return rgb_img

if __name__ == "__main__":
    create_roadmap()
