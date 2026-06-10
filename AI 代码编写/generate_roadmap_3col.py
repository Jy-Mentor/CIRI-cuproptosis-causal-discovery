# -*- coding: utf-8 -*-
from PIL import Image, ImageDraw, ImageFont
import os

OUTPUT_PATH = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\AI 代码编写\技术路线图_三栏.png"

WIDTH, HEIGHT = 2000, 1600
MARGIN = 60

COL_L_X = 100
COL_C_X = 700
COL_R_X = 1400

BOX_WIDTH = 500
BOX_HEIGHT = 120
BOX_GAP = 40
ARROW_LEN = 35

def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

colors = {
    'bg': '#FFFFFF',
    'frame_fill': '#f8f9fa',
    'frame_stroke': '#1e3a5f',
    'content_fill': '#ffffff',
    'content_stroke': '#b87333',
    'method_fill': '#f1f3f4',
    'method_stroke': '#64748b',
    'title_fill': '#1e3a5f',
    'text_dark': '#1e3a5f',
    'text_light': '#ffffff',
    'text_gray': '#334155',
    'arrow': '#1e3a5f',
    'dash_arrow': '#64748b'
}

left_items = [
    {'title': '科学问题提出', 'subtitle': '桂艾BCP如何调控\n铜死亡防治卒中？', 'type': 'frame'},
    {'title': '机制假说构建', 'subtitle': 'AGE-RAGE→SLC31A1\n→FDX1通路', 'type': 'frame'},
    {'title': '实验验证', 'subtitle': '计算→体内→体外', 'type': 'frame'},
    {'title': '机制阐明', 'subtitle': '论文/专利产出', 'type': 'frame'}
]

center_items = [
    {'title': '阶段一', 'subtitle': '生物信息学预测', 'detail': 'SwissTargetPrediction + PC算法 + MR验证', 'type': 'content'},
    {'title': '阶段二', 'subtitle': 'MCAO/R大鼠验证', 'detail': '线栓法 + mNSS + TTM挽救实验', 'type': 'content'},
    {'title': '阶段三', 'subtitle': 'OGD/R细胞机制', 'detail': 'SH-SY5Y + FPS-ZM1 + siRNA干预', 'type': 'content'},
    {'title': '阶段四', 'subtitle': '数据整合产出', 'detail': '多组学 + Western Blot + 论文撰写', 'type': 'content'}
]

right_items = [
    {'title': '网络药理学', 'subtitle': '靶点预测 & PPI网络', 'type': 'method'},
    {'title': '因果推断', 'subtitle': 'PC共定位 & FinnGen MR', 'type': 'method'},
    {'title': '分子生物学', 'subtitle': '分子对接 & ICP-MS', 'type': 'method'},
    {'title': '细胞功能', 'subtitle': 'CCK-8 & 免疫荧光', 'type': 'method'}
]

def get_box_style(box_type):
    if box_type == 'frame':
        return {'fill': colors['frame_fill'], 'stroke': colors['frame_stroke']}
    elif box_type == 'content':
        return {'fill': colors['content_fill'], 'stroke': colors['content_stroke']}
    elif box_type == 'method':
        return {'fill': colors['method_fill'], 'stroke': colors['method_stroke']}
    return {'fill': '#ffffff', 'stroke': '#000000'}

def draw_rounded_rect(draw, box, radius, fill, outline, width=2):
    draw.rounded_rectangle([box[0], box[1], box[2], box[3]], radius=radius, fill=fill, outline=outline, width=width)

def create_roadmap():
    img = Image.new('RGB', (WIDTH, HEIGHT), colors['bg'])
    draw = ImageDraw.Draw(img)
    
    try:
        title_font = ImageFont.truetype("msyh.ttc", 36)
        phase_font = ImageFont.truetype("msyh.ttc", 24)
        subtitle_font = ImageFont.truetype("msyh.ttc", 18)
        detail_font = ImageFont.truetype("msyh.ttc", 14)
    except:
        title_font = ImageFont.load_default()
        phase_font = ImageFont.load_default()
        subtitle_font = ImageFont.load_default()
        detail_font = ImageFont.load_default()
    
    title_y = 40
    draw.text((WIDTH//2, title_y), "桂艾β-石竹烯调控AGE-RAGE-铜死亡通路干预脑缺血再灌注损伤机制研究", 
              font=title_font, fill=hex_to_rgb(colors['title_fill']), anchor="mm")
    
    col_y_start = 120
    
    draw.text((COL_L_X + BOX_WIDTH//2, col_y_start), "研究框架", 
              font=phase_font, fill=hex_to_rgb(colors['frame_stroke']), anchor="mm")
    draw.text((COL_C_X + BOX_WIDTH//2, col_y_start), "研究内容", 
              font=phase_font, fill=hex_to_rgb(colors['content_stroke']), anchor="mm")
    draw.text((COL_R_X + BOX_WIDTH//2, col_y_start), "研究方法", 
              font=phase_font, fill=hex_to_rgb(colors['method_stroke']), anchor="mm")
    
    y = col_y_start + 40
    
    for i in range(4):
        box_y = y + i * (BOX_HEIGHT + BOX_GAP)
        
        style = get_box_style(left_items[i]['type'])
        draw_rounded_rect(draw, [COL_L_X, box_y, COL_L_X + BOX_WIDTH, box_y + BOX_HEIGHT], 
                        15, hex_to_rgb(style['fill']), hex_to_rgb(style['stroke']), 3)
        
        draw.text((COL_L_X + 25, box_y + 15), left_items[i]['title'], 
                  font=subtitle_font, fill=hex_to_rgb(colors['text_dark']))
        
        draw.text((COL_L_X + 25, box_y + 55), left_items[i]['subtitle'], 
                  font=detail_font, fill=hex_to_rgb(colors['text_gray']))
    
    for i in range(4):
        box_y = y + i * (BOX_HEIGHT + BOX_GAP)
        
        style = get_box_style(center_items[i]['type'])
        draw_rounded_rect(draw, [COL_C_X, box_y, COL_C_X + BOX_WIDTH, box_y + BOX_HEIGHT], 
                        15, hex_to_rgb(style['fill']), hex_to_rgb(style['stroke']), 3)
        
        draw.text((COL_C_X + 25, box_y + 12), center_items[i]['title'], 
                  font=subtitle_font, fill=hex_to_rgb(colors['content_stroke']))
        
        draw.text((COL_C_X + 25, box_y + 38), center_items[i]['subtitle'], 
                  font=subtitle_font, fill=hex_to_rgb(colors['text_dark']))
        
        draw.text((COL_C_X + 25, box_y + 70), center_items[i]['detail'], 
                  font=detail_font, fill=hex_to_rgb(colors['text_gray']))
    
    for i in range(4):
        box_y = y + i * (BOX_HEIGHT + BOX_GAP)
        
        style = get_box_style(right_items[i]['type'])
        draw_rounded_rect(draw, [COL_R_X, box_y, COL_R_X + BOX_WIDTH, box_y + BOX_HEIGHT], 
                        15, hex_to_rgb(style['fill']), hex_to_rgb(style['stroke']), 2)
        
        draw.text((COL_R_X + 25, box_y + 25), right_items[i]['title'], 
                  font=subtitle_font, fill=hex_to_rgb(colors['text_gray']))
        
        draw.text((COL_R_X + 25, box_y + 60), right_items[i]['subtitle'], 
                  font=detail_font, fill=hex_to_rgb(colors['text_gray']))
    
    center_x = WIDTH // 2
    
    for i in range(3):
        from_y = y + (i + 1) * (BOX_HEIGHT + BOX_GAP) - BOX_GAP
        to_y = y + (i + 2) * (BOX_HEIGHT + BOX_GAP) - BOX_GAP
        
        draw.line([(COL_L_X + BOX_WIDTH//2, from_y + BOX_HEIGHT), 
                  (COL_L_X + BOX_WIDTH//2, to_y)], 
                  fill=hex_to_rgb(colors['arrow']), width=3)
        
        draw.polygon([
            (COL_L_X + BOX_WIDTH//2, to_y),
            (COL_L_X + BOX_WIDTH//2 - 8, to_y - 12),
            (COL_L_X + BOX_WIDTH//2 + 8, to_y - 12)
        ], fill=hex_to_rgb(colors['arrow']))
    
    for i in range(3):
        from_y = y + (i + 1) * (BOX_HEIGHT + BOX_GAP) - BOX_GAP
        to_y = y + (i + 2) * (BOX_HEIGHT + BOX_GAP) - BOX_GAP
        
        draw.line([(COL_C_X + BOX_WIDTH//2, from_y + BOX_HEIGHT), 
                  (COL_C_X + BOX_WIDTH//2, to_y)], 
                  fill=hex_to_rgb(colors['content_stroke']), width=3)
        
        draw.polygon([
            (COL_C_X + BOX_WIDTH//2, to_y),
            (COL_C_X + BOX_WIDTH//2 - 8, to_y - 12),
            (COL_C_X + BOX_WIDTH//2 + 8, to_y - 12)
        ], fill=hex_to_rgb(colors['content_stroke']))
    
    for i in range(3):
        from_y = y + (i + 1) * (BOX_HEIGHT + BOX_GAP) - BOX_GAP
        to_y = y + (i + 2) * (BOX_HEIGHT + BOX_GAP) - BOX_GAP
        
        draw.line([(COL_R_X + BOX_WIDTH//2, from_y + BOX_HEIGHT), 
                  (COL_R_X + BOX_WIDTH//2, to_y)], 
                  fill=hex_to_rgb(colors['method_stroke']), width=2)
        
        draw.polygon([
            (COL_R_X + BOX_WIDTH//2, to_y),
            (COL_R_X + BOX_WIDTH//2 - 6, to_y - 10),
            (COL_R_X + BOX_WIDTH//2 + 6, to_y - 10)
        ], fill=hex_to_rgb(colors['method_stroke']))
    
    draw.line([(COL_L_X + BOX_WIDTH, y + BOX_HEIGHT//2), 
              (COL_C_X - 10, y + BOX_HEIGHT//2)], 
              fill=hex_to_rgb(colors['dash_arrow']), width=2)
    draw.line([(COL_C_X + BOX_WIDTH, y + BOX_HEIGHT//2), 
              (COL_R_X - 10, y + BOX_HEIGHT//2)], 
              fill=hex_to_rgb(colors['dash_arrow']), width=2)
    
    draw.line([(COL_L_X + BOX_WIDTH, y + (BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - BOX_GAP), 
              (COL_C_X - 10, y + (BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - BOX_GAP)], 
              fill=hex_to_rgb(colors['dash_arrow']), width=2)
    draw.line([(COL_C_X + BOX_WIDTH, y + (BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - BOX_GAP), 
              (COL_R_X - 10, y + (BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - BOX_GAP)], 
              fill=hex_to_rgb(colors['dash_arrow']), width=2)
    
    draw.line([(COL_L_X + BOX_WIDTH, y + 2*(BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - 2*BOX_GAP), 
              (COL_C_X - 10, y + 2*(BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - 2*BOX_GAP)], 
              fill=hex_to_rgb(colors['dash_arrow']), width=2)
    draw.line([(COL_C_X + BOX_WIDTH, y + 2*(BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - 2*BOX_GAP), 
              (COL_R_X - 10, y + 2*(BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - 2*BOX_GAP)], 
              fill=hex_to_rgb(colors['dash_arrow']), width=2)
    
    draw.line([(COL_L_X + BOX_WIDTH, y + 3*(BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - 3*BOX_GAP), 
              (COL_C_X - 10, y + 3*(BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - 3*BOX_GAP)], 
              fill=hex_to_rgb(colors['dash_arrow']), width=2)
    draw.line([(COL_C_X + BOX_WIDTH, y + 3*(BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - 3*BOX_GAP), 
              (COL_R_X - 10, y + 3*(BOX_HEIGHT + BOX_GAP) + BOX_HEIGHT//2 - 3*BOX_GAP)], 
              fill=hex_to_rgb(colors['dash_arrow']), width=2)
    
    img.save(OUTPUT_PATH)
    print(f"三栏式技术路线图已保存至: {OUTPUT_PATH}")
    return img

if __name__ == "__main__":
    create_roadmap()
