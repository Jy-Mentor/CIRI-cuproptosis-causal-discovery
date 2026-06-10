import gzip
import os

src_dir = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创"

gz_files = [f for f in os.listdir(src_dir) if f.endswith('.gz')]
print(f"找到 {len(gz_files)} 个 .gz 文件")

for i, gz_file in enumerate(gz_files):
    out_file = gz_file.replace('.gz', '')
    gz_path = os.path.join(src_dir, gz_file)
    out_path = os.path.join(src_dir, out_file)

    with gzip.open(gz_path, 'rb') as f_in:
        with open(out_path, 'wb') as f_out:
            f_out.write(f_in.read())

    if (i + 1) % 10 == 0:
        print(f"已解压 {i + 1}/{len(gz_files)}")

print(f"解压完成！共 {len(gz_files)} 个文件")
