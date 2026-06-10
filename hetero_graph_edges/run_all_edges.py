# -*- coding: utf-8 -*-
"""
异构图多维度边提取主控脚本
=====================================================
一键运行所有边提取和数据构建步骤, 生成异构图所需的全部数据文件。

输出目录: D:/反向网络药理学/GAT拓展维度

执行顺序:
  1. extract_tf_target_edges.py      → tf_target_edges.txt + tf_target_nodes.csv
  2. extract_gene_pathway_edges.py   → gene_pathway_edges.txt + pathway_nodes.csv
  3. compute_pathway_features.py     → pathway_features.npy (更新 pathway_nodes.csv)
  4. build_disease_features.py       → disease_features.csv + disease_features.npy
  5. extract_gene_methylation_edges.py → gene_methylation_edges.txt (可选)
  6. extract_gene_mirna_edges.py     → gene_mirna_edges.txt (可选)
  7. extract_archs4_coexp_edges.py   → gene_coexp_edges.txt (ARCHS4共表达)

使用方式:
  python run_all_edges.py              # 运行全部步骤
  python run_all_edges.py --skip-optional  # 跳过甲基化和miRNA (可选边)
  python run_all_edges.py --step 2     # 只运行步骤2
  python run_all_edges.py --step 7     # 只运行ARCHS4共表达查询

作者: 优化版 v2.0
日期: 2026-05-31
"""

import os
import sys
import subprocess
import argparse
import time
from pathlib import Path

# ============================================================
# 0. 配置
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")

STEPS = [
    {
        'id': 1,
        'name': 'TF-靶基因调控边',
        'script': 'extract_tf_target_edges.py',
        'outputs': ['tf_target_edges.txt', 'tf_target_nodes.csv'],
        'required': True,
        'desc': '从 TRRUST 数据库提取 TF→靶基因调控关系',
    },
    {
        'id': 2,
        'name': '基因-通路关联边',
        'script': 'extract_gene_pathway_edges.py',
        'outputs': ['gene_pathway_edges.txt', 'pathway_nodes.csv'],
        'required': True,
        'desc': '从 MSigDB/Reactome/KEGG 提取基因→通路关联',
    },
    {
        'id': 3,
        'name': '通路节点特征',
        'script': 'compute_pathway_features.py',
        'outputs': ['pathway_features.npy', 'pathway_feature_names.txt'],
        'required': True,
        'desc': '计算通路内基因嵌入均值作为通路特征',
        'depends_on': [2],
    },
    {
        'id': 4,
        'name': '疾病语义嵌入',
        'script': 'build_disease_features.py',
        'outputs': ['disease_features.csv', 'disease_features.npy'],
        'required': True,
        'desc': '构建 CIRI 疾病维度的语义嵌入 (基因均值法)',
    },
    {
        'id': 5,
        'name': '基因-甲基化关联',
        'script': 'extract_gene_methylation_edges.py',
        'outputs': ['gene_methylation_edges.txt'],
        'required': False,
        'desc': '从 EWAS Atlas 提取基因↔CpG甲基化关联 (需手动下载)',
    },
    {
        'id': 6,
        'name': 'miRNA-靶基因关联',
        'script': 'extract_gene_mirna_edges.py',
        'outputs': ['gene_mirna_edges.txt'],
        'required': False,
        'desc': '从 miRTarBase 提取 miRNA→靶基因调控关系 (需手动下载)',
    },
    {
        'id': 7,
        'name': 'ARCHS4 共表达边',
        'script': 'extract_archs4_coexp_edges.py',
        'outputs': ['gene_coexp_edges.txt'],
        'required': True,
        'desc': '从 ARCHS4 数据库查询基因共表达相关性 (|corr|>0.7)',
    },
]


def check_outputs(step):
    """检查步骤的输出文件是否已存在"""
    for output_file in step['outputs']:
        output_path = OUTPUT_DIR / output_file
        if not output_path.exists():
            return False, f"缺少: {output_file}"
    return True, "全部就绪"


def run_step(step, force=False):
    """运行单个步骤"""
    print(f"\n{'#'*70}")
    print(f"# 步骤 {step['id']}: {step['name']}")
    print(f"# {step['desc']}")
    print(f"{'#'*70}")

    ready, msg = check_outputs(step)
    if ready and not force:
        print(f"[SKIP] 输出文件已存在: {msg}")
        print(f"       使用 --force 强制重新运行")
        return True

    if not ready and force:
        print(f"[INFO] 重新运行 (--force)")

    script_path = SCRIPT_DIR / step['script']
    if not script_path.exists():
        print(f"[ERROR] 脚本不存在: {script_path}")
        return False

    print(f"[RUN] {script_path}")
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(SCRIPT_DIR),
            capture_output=False,
            timeout=1800,
        )
        if result.returncode != 0:
            print(f"[FAIL] 步骤 {step['id']} 返回码: {result.returncode}")
            if step['required']:
                return False
            else:
                print(f"[WARN] 可选步骤失败, 跳过")
                return True
    except subprocess.TimeoutExpired:
        print(f"[FAIL] 步骤 {step['id']} 超时 (30分钟)")
        if step['required']:
            return False
        else:
            return True
    except Exception as e:
        print(f"[FAIL] 步骤 {step['id']} 异常: {e}")
        if step['required']:
            return False
        else:
            return True

    ready, msg = check_outputs(step)
    if ready:
        print(f"[OK] 步骤 {step['id']} 完成: {msg}")
    else:
        print(f"[WARN] 步骤 {step['id']} 运行完毕但输出缺失: {msg}")

    return ready


def print_summary_header():
    """打印运行摘要表头"""
    print(f"\n{'='*70}")
    print(f"异构图多维度边提取 - 运行摘要")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"开始时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*70}")


def print_summary_footer(results, start_time):
    """打印运行摘要结果"""
    elapsed = time.time() - start_time

    print(f"\n{'='*70}")
    print(f"运行完成")
    print(f"耗时: {elapsed:.0f} 秒 ({elapsed/60:.1f} 分钟)")
    print(f"{'='*70}")

    print(f"\n{'步骤':<4} {'名称':<25} {'状态':<10} {'输出文件'}")
    print(f"{'-'*70}")
    for step, ok in results:
        status = '✓ 完成' if ok else '✗ 失败'
        outputs = ', '.join(step['outputs'][:3])
        if len(step['outputs']) > 3:
            outputs += f' +{len(step["outputs"])-3}'
        optional = ' (可选)' if not step['required'] else ''
        print(f"{step['id']:<4} {step['name']+optional:<25} {status:<10} {outputs}")

    n_ok = sum(1 for _, ok in results if ok)
    n_req = sum(1 for s, _ in results if s['required'])
    n_req_ok = sum(1 for s, ok in results if s['required'] and ok)
    print(f"\n必需步骤: {n_req_ok}/{n_req} 完成")
    print(f"可选步骤: {n_ok - n_req_ok}/{len(results) - n_req} 完成")

    print(f"\n生成的文件:")
    for f in sorted(OUTPUT_DIR.glob('*')):
        if f.is_file():
            size_kb = f.stat().st_size / 1024
            print(f"  {f.name:<40} {size_kb:>8.1f} KB")

    if n_req_ok == n_req:
        print(f"\n[OK] 所有必需步骤完成! 可以开始构建异构图。")
    else:
        failed = [s['name'] for s, ok in results if s['required'] and not ok]
        print(f"\n[WARN] 以下必需步骤失败: {', '.join(failed)}")


def main():
    parser = argparse.ArgumentParser(
        description='异构图多维度边提取 - 一键运行',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python run_all_edges.py                    # 运行全部步骤
  python run_all_edges.py --skip-optional     # 跳过甲基化/miRNA
  python run_all_edges.py --step 1 --step 2  # 只运行步骤1和2
  python run_all_edges.py --force            # 强制重新运行所有步骤
        """
    )
    parser.add_argument('--skip-optional', action='store_true',
                        help='跳过可选步骤 (甲基化 + miRNA)')
    parser.add_argument('--step', type=int, action='append', dest='steps',
                        help='只运行指定步骤 (可多次使用)')
    parser.add_argument('--force', action='store_true',
                        help='强制重新运行, 即使输出文件已存在')
    parser.add_argument('--dry-run', action='store_true',
                        help='仅显示将要执行的步骤, 不实际运行')
    parser.add_argument('--list', action='store_true',
                        help='列出所有步骤及其状态')
    args = parser.parse_args()

    if args.list:
        print(f"\n{'步骤':<4} {'名称':<25} {'必需':<6} {'状态'}")
        print(f"{'-'*60}")
        for step in STEPS:
            ready, msg = check_outputs(step)
            status = '✓ 已就绪' if ready else f'✗ {msg}'
            required = '是' if step['required'] else '否'
            print(f"{step['id']:<4} {step['name']:<25} {required:<6} {status}")
        print(f"\n输出目录: {OUTPUT_DIR}")
        return

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print(f"\n[DRY RUN] 将执行以下步骤:")
        for step in STEPS:
            if args.skip_optional and not step['required']:
                print(f"  [SKIP] 步骤 {step['id']}: {step['name']} (可选, 已跳过)")
                continue
            if args.steps and step['id'] not in args.steps:
                print(f"  [SKIP] 步骤 {step['id']}: {step['name']} (未选中)")
                continue
            ready, _ = check_outputs(step)
            force_mark = ' [FORCE]' if args.force and ready else ''
            skip_mark = ' [ALREADY DONE]' if ready and not args.force else ''
            print(f"  [RUN]  步骤 {step['id']}: {step['name']}{force_mark}{skip_mark}")
        return

    start_time = time.time()
    print_summary_header()

    results = []
    for step in STEPS:
        if args.skip_optional and not step['required']:
            print(f"\n[SKIP] 步骤 {step['id']}: {step['name']} (可选, 已跳过)")
            results.append((step, True))
            continue

        if args.steps and step['id'] not in args.steps:
            print(f"\n[SKIP] 步骤 {step['id']}: {step['name']} (未选中)")
            results.append((step, True))
            continue

        ok = run_step(step, force=args.force)
        results.append((step, ok))

        if not ok and step['required']:
            print(f"\n[ABORT] 必需步骤 {step['id']} 失败, 停止执行")
            break

    print_summary_footer(results, start_time)


if __name__ == "__main__":
    main()