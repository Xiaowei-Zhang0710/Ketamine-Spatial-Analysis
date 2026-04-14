import os, sys
sys.path.append('/data2/usr/yangmy_conda/Tangram')
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Wedge
from matplotlib.cm import get_cmap
from matplotlib.lines import Line2D
import torch
import tangram as tg
import gc

# ---------------- 参数设置 ----------------
deg_file = "/data2/Project/anding/test_20250923/seurat_cluster_deg.csv"
sc_h5ad = "/data2/Project/anding/test_20250923/subtype.h5ad"
sp_h5ad_dir = "/data2/Project/anding/test_20250923/ref_h5ad/"
output_dir = "/data2/Project/anding/test_20250923/output/"
#skip_keywords = ["", ""]  # 带有这些字符串的文件将被跳过
scatter_point_size = 20
pie_radius = 60

os.makedirs(output_dir, exist_ok=True)

# ---------------- 读取 DEG 基因 ----------------
df_genes_con = pd.read_csv(deg_file, index_col=0)
markers_con = df_genes_con.values.flatten().tolist()

# ---------------- 读取单细胞参考 ----------------
ad_sc1 = sc.read_h5ad(sc_h5ad)
sc.pp.normalize_total(ad_sc1)

# ---------------- 第一阶段：Tangram 映射 + 散点图 ----------------
print("=== 第一阶段：Tangram 映射 + 散点图 ===")
for fname in os.listdir(sp_h5ad_dir):
    if not fname.endswith(".h5ad"):
        continue
    if "subtype" in fname:  # 避免重复读 ad_sc1
        continue
    #if any(kw in fname for kw in skip_keywords):
    #    print(f"跳过 {fname}")
    #    continue

    basename = os.path.splitext(fname)[0]
    out_dir = os.path.join(output_dir, basename)
    os.makedirs(out_dir, exist_ok=True)

    print(f"处理 {fname} ...")
    ad_sp1 = sc.read_h5ad(os.path.join(sp_h5ad_dir, fname))
    ad_sp1.var_names = ad_sp1.var_names.astype(str).str.capitalize()

    # Tangram 映射
    tg.pp_adatas(ad_sc1, ad_sp1, genes=markers_con)
    mapping = tg.map_cells_to_space(ad_sc1, ad_sp1, mode='cells')

    if mapping.shape[0] != ad_sp1.shape[0]:
        mapping = mapping.T
    mapping.obs_names = ad_sp1.obs_names

    # 计算细胞坐标
    cell_to_spot_idx = mapping.X.argmax(axis=0)
    assigned_spots = mapping.obs_names[cell_to_spot_idx]
    spot_coords_df = ad_sp1.obs[["centroid_1", "centroid_2"]].copy()
    cell_coords = spot_coords_df.loc[assigned_spots].copy()
    cell_coords["cell"] = mapping.var_names
    cell_coords.index = mapping.var_names
    cell_coords["celltype.stim"] = ad_sc1.obs.loc[cell_coords.index, "celltype2"]

    # 保存散点图数据和坐标
    cell_coords.to_csv(os.path.join(out_dir, f"{basename}_scatter_data.csv"), index=True)
    cell_coords.to_csv(os.path.join(out_dir, f"{basename}_cells_to_spcoord.txt"), sep="\t", index=False)

    # 绘制散点图
    plt.figure(figsize=(12, 10))
    sns.scatterplot(
        data=cell_coords,
        x="centroid_1", y="centroid_2",
        hue="celltype.stim",
        s=scatter_point_size, alpha=0.6, palette="tab20"
    )
    plt.gca().invert_yaxis()
    plt.title(f"Mapped cells - {basename}")
    plt.legend(markerscale=2, bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{basename}_scatter.pdf"), dpi=300, bbox_inches="tight")
    plt.close()

    # 释放内存
    del ad_sp1, mapping, cell_coords
    gc.collect()

print("=== 第一阶段完成 ===\n")

# ---------------- 第二阶段：饼图绘制 ----------------
print("=== 第二阶段：饼图绘制 ===")
for fname in os.listdir(output_dir):
    sample_dir = os.path.join(output_dir, fname)
    if not os.path.isdir(sample_dir):
        continue

    pie_data_file = os.path.join(sample_dir, f"{fname}_scatter_data.csv")
    if not os.path.exists(pie_data_file):
        continue

    print(f"绘制饼图：{fname}")
    cell_coords = pd.read_csv(pie_data_file, index_col=0)

    # 按空间点统计比例
    grouped = cell_coords.groupby(['centroid_1','centroid_2','celltype.stim']).size().reset_index(name='count')
    total_counts = grouped.groupby(['centroid_1','centroid_2'])['count'].transform('sum')
    grouped['pct'] = grouped['count'] / total_counts

    # 转为宽表格
    count_wide = grouped.pivot_table(index=['centroid_1','centroid_2'], columns='celltype.stim', values='count', fill_value=0)
    pct_wide = grouped.pivot_table(index=['centroid_1','centroid_2'], columns='celltype.stim', values='pct', fill_value=0)
    count_wide.columns = [f'count_{c}' for c in count_wide.columns]
    pct_wide.columns = [f'pct_{c}' for c in pct_wide.columns]
    final_df = pd.concat([count_wide, pct_wide], axis=1).reset_index()
    final_df_filtered = final_df[(count_wide.sum(axis=1) != 0)].copy()

    # 保存饼图数据
    final_df_filtered.to_csv(os.path.join(sample_dir, f"{fname}_pie_data.csv"), index=False)

    # 绘制饼图
    pct_cols = [c for c in final_df_filtered.columns if c.startswith('pct_')]
    fig, ax = plt.subplots(figsize=(12, 10))
    cmap = get_cmap('tab10', len(pct_cols))
    color_map = {col: cmap(i) for i, col in enumerate(pct_cols)}

    for _, row in final_df_filtered.iterrows():
        x, y = row['centroid_1'], row['centroid_2']
        sizes = np.array([row[c] for c in pct_cols])
        sizes = sizes / sizes.sum() if sizes.sum() > 0 else sizes
        start_angle = 0
        for frac, col in zip(sizes, pct_cols):
            if frac == 0:
                continue
            theta1 = start_angle * 360
            theta2 = (start_angle + frac) * 360
            wedge = Wedge((x, y), pie_radius, theta1, theta2, color=color_map[col], edgecolor='black', linewidth=0.5)
            ax.add_patch(wedge)
            start_angle += frac

    ax.set_xlim(final_df_filtered['centroid_1'].min()-pie_radius*2, final_df_filtered['centroid_1'].max()+pie_radius*2)
    ax.set_ylim(final_df_filtered['centroid_2'].min()-pie_radius*2, final_df_filtered['centroid_2'].max()+pie_radius*2)
    ax.set_aspect('equal')
    ax.set_xlabel('centroid_1')
    ax.set_ylabel('centroid_2')
    ax.set_title(f"Spatial pie charts - {fname}")

    legend_elements = [Line2D([0],[0], marker='o', color='w', label=col,
                              markerfacecolor=color_map[col], markersize=10)
                       for col in pct_cols]
    ax.legend(handles=legend_elements, title='Celltype Percent', bbox_to_anchor=(1.05,1), loc='upper left')

    plt.savefig(os.path.join(sample_dir, f"{fname}_spatial_piecharts.pdf"), bbox_inches='tight')
    plt.close()
    del cell_coords, final_df, final_df_filtered
    gc.collect()

print("=== 所有样本饼图绘制完成 ===")
