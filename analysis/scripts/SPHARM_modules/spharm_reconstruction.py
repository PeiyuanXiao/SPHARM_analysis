# ==============================================================================
# spharm_reconstruction.py
#
# 功能：
#   1. 全谱重建          — EXP轴端标本 + IM + SDG + EXP全部逐件
#   2. 分阶重建          — EXP轴端标本
#   3. 四象限轴端对比图  — EXP轴端标本
#   4. 分类型平均重建    — EXP均值 + IM逐个对比
#   5. 轴连续变化轨迹    — Axis1/2 × 形态/疤痕，多视角渲染
# ==============================================================================

import os
os.environ['PYVISTA_OFF_SCREEN'] = 'true'

import numpy as np
import pandas as pd
import pyshtools as pysh
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap

pv.global_theme.window_size = [1200, 1200]
pv.global_theme.background = 'white'
pv.global_theme.font.color = 'black'

# ── 路径配置 ──────────────────────────────────────────────────────────────────
MORPH_CSV  = "analysis/data/derived_data/SPHARM_morphology.csv"
SCAR_CSV   = "analysis/data/derived_data/SPHARM_direction.csv"
COORD_CSV  = "analysis/data/derived_data/EXP_CIA_coords_full.csv"
OUT_DIR    = "analysis/output/figures/reconstruction"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, "EXP_individual"), exist_ok=True)

# ── 全局常量 ──────────────────────────────────────────────────────────────────
LMAX = 20

AXIS_EXTREMES = {
    'Axis1+\n(EXP30 Lev.Nubian)':     'EXP30_Levallois Nubian',
    'Axis1-\n(EXP45 Conical uni.)':    'EXP45_Conical unidirectional',
    'Axis2+\n(EXP06 Lev.convergent)':  'EXP06_Levallois convergent',
    'Axis2-\n(EXP31 Lev.laminar)':     'EXP31_Levallois laminar',
}

DEGREES_TO_SHOW = [1, 2, 3, 4, 5]

TYPOLOGY_ORDER = [
    'Unidirectional', 'Bidirectional', 'Levallois',
    'Discoid', 'Multiplatform'
]

TRAJECTORY_N_STEPS   = 7
TRAJECTORY_BANDWIDTH = 0.4
TRAJECTORY_PERCENTILE = 10

# 多视角相机位置
CAMERA_VIEWS = {
    'iso':   'iso',
    'top':   'zy',
    'front': 'xz',
}

EDGE_COLOR  = '#333333'
EDGE_WIDTH  = 0.22

CMAP = LinearSegmentedColormap.from_list(
    'custom_spharm',
    ['#5C7F71', '#F5EDDC', '#802520'],
    N=256
)

# ── 辅助函数 ──────────────────────────────────────────────────────────────────

def load_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df['ID'] = df['ID'].str.strip()
    return df


def get_coeff_array(row: pd.Series) -> np.ndarray:
    coeff_cols = sorted(
        [c for c in row.index if c.startswith('coeff_')],
        key=lambda x: int(x.split('_')[1])
    )
    vals = row[coeff_cols].values.astype(float)
    cilm = vals.reshape(2, LMAX + 1, LMAX + 1)
    return cilm


def cilm_to_grid(cilm: np.ndarray, lmax_out: int = None) -> np.ndarray:
    if lmax_out is None:
        lmax_out = LMAX
    clm = pysh.SHCoeffs.from_array(cilm, normalization='4pi',
                                    csphase=1, lmax=LMAX)
    clm_pad = clm.pad(lmax=lmax_out)
    grid = clm_pad.expand(grid='DH')
    return np.real(grid.data)


def cilm_single_degree(cilm: np.ndarray, l_only: int) -> np.ndarray:
    masked = np.zeros_like(cilm)
    masked[:, l_only, :l_only + 1] = cilm[:, l_only, :l_only + 1]
    return masked


def smooth_poles(grid_data: np.ndarray, n_rows: int = 4) -> np.ndarray:
    g = grid_data.copy()
    n_lat = g.shape[0]
    for i in range(1, n_rows + 1):
        w = i / (n_rows + 1)
        g[i - 1, :]     = (1 - w) * g[i - 1, :].mean()     + w * g[i - 1, :]
        g[n_lat - i, :] = (1 - w) * g[n_lat - i, :].mean() + w * g[n_lat - i, :]
    return g


def grid_to_mesh(grid_data: np.ndarray):
    grid_data = smooth_poles(grid_data, n_rows=4)
    n_lat, n_lon = grid_data.shape
    theta = np.linspace(0, np.pi,   n_lat, endpoint=True)
    phi   = np.linspace(0, 2*np.pi, n_lon, endpoint=True)
    theta_g, phi_g = np.meshgrid(theta, phi, indexing='ij')
    r = grid_data
    x = (r * np.sin(theta_g) * np.cos(phi_g)).T
    y = (r * np.sin(theta_g) * np.sin(phi_g)).T
    z = (r * np.cos(theta_g)).T
    struct = pv.StructuredGrid(x, y, z)
    tri    = struct.triangulate()
    r_flat = r.flatten()
    return tri, r_flat


def add_mesh_with_wireframe(pl, mesh, scalars, smin=None, smax=None):
    """
    实体面 + 独立线框叠加，解决 off_screen 下 line_width 不生效的问题。
    """
    if smin is None:
        smin = scalars.min()
    if smax is None:
        smax = scalars.max()
    # 实体面，不显示边线
    pl.add_mesh(
        mesh,
        scalars=scalars,
        clim=[smin, smax],
        cmap=CMAP,
        show_edges=False,
        lighting=True,
        show_scalar_bar=False,
    )
    # 独立线框层，单独控制线宽和颜色
    pl.add_mesh(
        mesh,
        style='wireframe',
        color=EDGE_COLOR,
        line_width=EDGE_WIDTH,
        opacity=1.0,
    )


def render_single(mesh, scalars, title: str, out_path: str,
                  camera_position='iso'):
    """单视角渲染，保存单张PNG。"""
    pl = pv.Plotter(off_screen=True, window_size=[900, 900])
    pl.background_color = 'white'
    add_mesh_with_wireframe(pl, mesh, scalars)
    pl.add_title(title, font_size=11, color='black')
    pl.camera_position = camera_position
    pl.camera.elevation = 20
    pl.screenshot(out_path)
    pl.close()
    print(f"  已保存：{out_path}")


def render_single_multiview(mesh, scalars, title: str, out_path: str):
    """
    多视角渲染：iso / top / front 三个角度横向拼排，保存单张PNG。
    """
    views = list(CAMERA_VIEWS.items())
    ncols = len(views)
    pl = pv.Plotter(
        off_screen=True,
        shape=(1, ncols),
        window_size=[ncols * 600, 600],
    )
    pl.background_color = 'white'
    smin, smax = scalars.min(), scalars.max()

    for col, (view_name, cam) in enumerate(views):
        pl.subplot(0, col)
        add_mesh_with_wireframe(pl, mesh, scalars, smin, smax)
        pl.add_title(f"{title}\n[{view_name}]", font_size=9, color='black')
        pl.camera_position = cam
        if cam == 'iso':
            pl.camera.elevation = 20

    pl.screenshot(out_path)
    pl.close()
    print(f"  已保存：{out_path}")


def render_panel(meshes_scalars: list, titles: list,
                 out_path: str, ncols: int = 4,
                 multiview: bool = False):
    """
    多面板渲染。
    multiview=True 时每个面板内嵌三视角（iso/top/front），
    此时实际列数 = ncols * 3，行数不变。
    """
    n = len(meshes_scalars)
    nrows = int(np.ceil(n / ncols))

    if multiview:
        actual_cols = ncols * len(CAMERA_VIEWS)
        panel_w, panel_h = 400, 400
    else:
        actual_cols = ncols
        panel_w, panel_h = 700, 700

    pl = pv.Plotter(
        off_screen=True,
        shape=(nrows, actual_cols),
        window_size=[actual_cols * panel_w, nrows * panel_h],
    )
    pl.background_color = 'white'

    views = list(CAMERA_VIEWS.items())

    for i, ((mesh, scalars), title) in enumerate(zip(meshes_scalars, titles)):
        base_row, base_col = divmod(i, ncols)
        smin, smax = scalars.min(), scalars.max()

        if multiview:
            for v_idx, (view_name, cam) in enumerate(views):
                col = base_col * len(views) + v_idx
                pl.subplot(base_row, col)
                add_mesh_with_wireframe(pl, mesh, scalars, smin, smax)
                col_title = f"{title}\n[{view_name}]" if v_idx == 0 else f"[{view_name}]"
                pl.add_title(col_title, font_size=8, color='black')
                pl.camera_position = cam
                if cam == 'iso':
                    pl.camera.elevation = 20
        else:
            pl.subplot(base_row, base_col)
            add_mesh_with_wireframe(pl, mesh, scalars, smin, smax)
            pl.add_title(title, font_size=10, color='black')
            pl.camera_position = 'iso'
            pl.camera.elevation = 20

    # 空白子图
    for j in range(n, nrows * ncols):
        base_row, base_col = divmod(j, ncols)
        if multiview:
            for v_idx in range(len(views)):
                col = base_col * len(views) + v_idx
                pl.subplot(base_row, col)
                pl.background_color = 'white'
        else:
            pl.subplot(base_row, base_col)
            pl.background_color = 'white'

    pl.screenshot(out_path)
    pl.close()
    print(f"  已保存：{out_path}")


# ── 主重建函数 ────────────────────────────────────────────────────────────────

def full_reconstruction(df: pd.DataFrame, specimen_ids: list,
                        label: str = 'morph'):
    """全谱重建：指定标本逐一重建，多视角。"""
    print(f"\n=== 全谱重建（{label}）===")
    for sid in specimen_ids:
        rows = df[df['ID'] == sid]
        if rows.empty:
            print(f"  [跳过] 找不到 ID：{sid}")
            continue
        row = rows.iloc[0]
        cilm = get_coeff_array(row)
        grid = cilm_to_grid(cilm)
        mesh, scalars = grid_to_mesh(grid)
        out = os.path.join(OUT_DIR, f"recon_full_{sid}_{label}.png")
        render_single_multiview(mesh, scalars,
                                title=f"{sid} ({label})", out_path=out)


def exp_individual_reconstruction(df: pd.DataFrame, label: str = 'morph'):
    """EXP 全部标本逐件重建，保存到 EXP_individual/ 子目录。"""
    exp_df = df[df['ID'].str.startswith('EXP')].copy()
    print(f"\n=== EXP 逐件重建（{label}，n={len(exp_df)}）===")
    out_subdir = os.path.join(OUT_DIR, "EXP_individual")

    for _, row in exp_df.iterrows():
        sid = row['ID']
        try:
            cilm = get_coeff_array(row)
            grid = cilm_to_grid(cilm)
            mesh, scalars = grid_to_mesh(grid)
            out = os.path.join(out_subdir, f"{sid}_{label}.png")
            render_single_multiview(mesh, scalars,
                                    title=f"{sid}\n({label})", out_path=out)
        except Exception as e:
            print(f"  [跳过] {sid}: {e}")


def per_degree_reconstruction(df: pd.DataFrame, specimen_ids: list,
                               label: str = 'morph',
                               degrees: list = None):
    if degrees is None:
        degrees = DEGREES_TO_SHOW
    print(f"\n=== 分阶重建（{label}）===")
    for sid in specimen_ids:
        rows = df[df['ID'] == sid]
        if rows.empty:
            print(f"  [跳过] 找不到 ID：{sid}")
            continue
        row = rows.iloc[0]
        cilm = get_coeff_array(row)

        items, titles = [], []
        grid_full = cilm_to_grid(cilm)
        m, s = grid_to_mesh(grid_full)
        items.append((m, s))
        titles.append("Full")

        for l in degrees:
            cilm_l = cilm_single_degree(cilm, l)
            grid_l = cilm_to_grid(cilm_l, lmax_out=l)
            m, s = grid_to_mesh(grid_l)
            items.append((m, s))
            titles.append(f"l={l}")

        out = os.path.join(OUT_DIR, f"recon_degrees_{sid}_{label}.png")
        render_panel(items, titles, out_path=out,
                     ncols=len(degrees) + 1, multiview=True)


def axis_extremes_panel(df: pd.DataFrame,
                        extremes: dict = None,
                        label: str = 'morph'):
    if extremes is None:
        extremes = AXIS_EXTREMES
    print(f"\n=== 轴端对比图（{label}）===")

    items, titles = [], []
    for panel_title, sid in extremes.items():
        rows = df[df['ID'] == sid]
        if rows.empty:
            print(f"  [跳过] 找不到 ID：{sid}")
            items.append((pv.Sphere(), np.array([0.0])))
            titles.append(f"{panel_title}\n[NOT FOUND]")
            continue
        row = rows.iloc[0]
        cilm = get_coeff_array(row)
        grid = cilm_to_grid(cilm)
        m, s = grid_to_mesh(grid)
        items.append((m, s))
        titles.append(panel_title)

    out = os.path.join(OUT_DIR, f"recon_axis_extremes_{label}.png")
    render_panel(items, titles, out_path=out, ncols=2, multiview=True)


def im_reconstruction_panel(df: pd.DataFrame, label: str = 'morph'):
    im_df = df[df['ID'].str.startswith('IM_')].copy()
    if im_df.empty:
        print(f"  [跳过] 没有找到 IM_ 标本（{label}）")
        return

    print(f"\n=== IM 理想模型重建（{label}，n={len(im_df)}）===")
    items, titles = [], []
    for _, row in im_df.iterrows():
        sid = row['ID']
        try:
            cilm = get_coeff_array(row)
            grid = cilm_to_grid(cilm)
            m, s = grid_to_mesh(grid)
            items.append((m, s))
            titles.append(sid.replace('IM_', '').replace('_', '\n'))
            print(f"  {sid} ✓")
        except Exception as e:
            print(f"  [跳过] {sid}: {e}")

    out = os.path.join(OUT_DIR, f"recon_IM_{label}.png")
    ncols = min(len(items), 4)
    render_panel(items, titles, out_path=out, ncols=ncols, multiview=True)


def sdg_reconstruction_panel(df: pd.DataFrame, label: str = 'morph'):
    sdg_df = df[df['ID'].str.startswith('SDG')].copy()
    if sdg_df.empty:
        print(f"  [跳过] 没有找到 SDG 标本（{label}）")
        return

    print(f"\n=== SDG 考古标本重建（{label}，n={len(sdg_df)}）===")
    items, titles = [], []
    for _, row in sdg_df.iterrows():
        sid = row['ID']
        try:
            cilm = get_coeff_array(row)
            grid = cilm_to_grid(cilm)
            m, s = grid_to_mesh(grid)
            items.append((m, s))
            titles.append(sid)
            print(f"  {sid} ✓")
        except Exception as e:
            print(f"  [跳过] {sid}: {e}")

    out = os.path.join(OUT_DIR, f"recon_SDG_{label}.png")
    ncols = min(len(items), 4)
    render_panel(items, titles, out_path=out, ncols=ncols, multiview=True)


def typology_mean_reconstruction(df: pd.DataFrame,
                                 label: str = 'morph',
                                 typology_col: str = 'Typology',
                                 min_n: int = 3):
    print(f"\n=== EXP 分类型平均重建（{label}）===")
    df = df.copy()
    df = df[df['ID'].str.startswith('EXP')]
    df[typology_col] = df[typology_col].str.strip()
    lev_mask = df[typology_col].str.contains('evallois', case=False, na=False)
    df.loc[lev_mask, typology_col] = 'Levallois'
    df = df[df[typology_col] != 'Biface']

    typologies = df[typology_col].value_counts()
    typologies = typologies[typologies >= min_n].index.tolist()
    typologies = [t for t in TYPOLOGY_ORDER if t in typologies] + \
                 [t for t in typologies if t not in TYPOLOGY_ORDER]

    items, titles = [], []
    for typ in typologies:
        sub = df[df[typology_col] == typ]
        print(f"  {typ}: n={len(sub)}")
        cilm_list = []
        for _, row in sub.iterrows():
            try:
                cilm_list.append(get_coeff_array(row))
            except Exception as e:
                print(f"    [跳过] {row['ID']}: {e}")
        if not cilm_list:
            continue
        cilm_mean = np.mean(cilm_list, axis=0)
        grid = cilm_to_grid(cilm_mean)
        m, s = grid_to_mesh(grid)
        items.append((m, s))
        titles.append(f"{typ}\n(n={len(cilm_list)})")

    out = os.path.join(OUT_DIR, f"recon_typology_mean_{label}.png")
    ncols = min(len(items), 5)
    render_panel(items, titles, out_path=out, ncols=ncols, multiview=True)


def typology_mean_with_im(df: pd.DataFrame,
                          label: str = 'morph',
                          typology_col: str = 'Typology',
                          min_n: int = 3):
    print(f"\n=== EXP 分类型均值 + IM 对比（{label}）===")
    df = df.copy()

    exp_df = df[df['ID'].str.startswith('EXP')].copy()
    exp_df[typology_col] = exp_df[typology_col].str.strip()
    lev_mask = exp_df[typology_col].str.contains('evallois', case=False, na=False)
    exp_df.loc[lev_mask, typology_col] = 'Levallois'
    exp_df = exp_df[exp_df[typology_col] != 'Biface']

    typologies = exp_df[typology_col].value_counts()
    typologies = typologies[typologies >= min_n].index.tolist()
    typologies = [t for t in TYPOLOGY_ORDER if t in typologies] + \
                 [t for t in typologies if t not in TYPOLOGY_ORDER]

    exp_items, exp_titles = [], []
    for typ in typologies:
        sub = exp_df[exp_df[typology_col] == typ]
        cilm_list = []
        for _, row in sub.iterrows():
            try:
                cilm_list.append(get_coeff_array(row))
            except Exception as e:
                print(f"    [跳过] {row['ID']}: {e}")
        if not cilm_list:
            continue
        cilm_mean = np.mean(cilm_list, axis=0)
        grid = cilm_to_grid(cilm_mean)
        m, s = grid_to_mesh(grid)
        exp_items.append((m, s))
        exp_titles.append(f"EXP\n{typ}\n(n={len(cilm_list)})")
        print(f"  EXP {typ}: n={len(cilm_list)}")

    im_df = df[df['ID'].str.startswith('IM_')].copy()
    im_items, im_titles = [], []
    for _, row in im_df.iterrows():
        sid = row['ID']
        try:
            cilm = get_coeff_array(row)
            grid = cilm_to_grid(cilm)
            m, s = grid_to_mesh(grid)
            im_items.append((m, s))
            im_titles.append('IM\n' + sid.replace('IM_', '').replace('_', '\n'))
            print(f"  IM {sid} ✓")
        except Exception as e:
            print(f"  [跳过] {sid}: {e}")

    ncols = max(len(exp_items), len(im_items), 1)
    all_items  = exp_items  + im_items
    all_titles = exp_titles + im_titles

    out = os.path.join(OUT_DIR, f"recon_typology_vs_IM_{label}.png")
    render_panel(all_items, all_titles, out_path=out,
                 ncols=ncols, multiview=True)


# ── 轴连续变化轨迹 ────────────────────────────────────────────────────────────

def gaussian_weights(positions: np.ndarray,
                     target: float,
                     bandwidth: float) -> np.ndarray:
    w = np.exp(-0.5 * ((positions - target) / bandwidth) ** 2)
    w_sum = w.sum()
    if w_sum < 1e-10:
        return np.zeros_like(w)
    return w / w_sum


def axis_trajectory(df_coeffs: pd.DataFrame,
                    df_coords: pd.DataFrame,
                    axis_col: str,
                    label: str,
                    axis_name: str,
                    n_steps: int = TRAJECTORY_N_STEPS,
                    bandwidth: float = TRAJECTORY_BANDWIDTH,
                    pct: float = TRAJECTORY_PERCENTILE):
    print(f"\n=== 轴轨迹：{axis_name} × {label}（列：{axis_col}）===")

    merged = df_coords[['ID', axis_col]].merge(
        df_coeffs, on='ID', how='inner'
    )
    if merged.empty:
        print(f"  [跳过] 合并后无数据")
        return None

    positions  = merged[axis_col].values
    lo = np.percentile(positions, pct)
    hi = np.percentile(positions, 100 - pct)
    sample_points = np.linspace(lo, hi, n_steps)

    items, titles = [], []
    for pt in sample_points:
        w = gaussian_weights(positions, pt, bandwidth)
        if w.max() < 1e-6:
            print(f"  [跳过] 位置 {pt:.2f}：权重过小")
            continue

        coeff_cols   = sorted(
            [c for c in merged.columns if c.startswith('coeff_')],
            key=lambda x: int(x.split('_')[1])
        )
        coeff_matrix = merged[coeff_cols].values.astype(float)
        cilm_flat    = (w[:, np.newaxis] * coeff_matrix).sum(axis=0)
        cilm         = cilm_flat.reshape(2, LMAX + 1, LMAX + 1)

        grid = cilm_to_grid(cilm)
        m, s = grid_to_mesh(grid)
        items.append((m, s))
        titles.append(f"{axis_col}\n= {pt:.2f}")
        print(f"  位置 {pt:.2f}：有效权重标本数 ≈ {int(1/(w**2).sum())}")

    out = os.path.join(OUT_DIR,
                       f"recon_trajectory_{axis_name}_{label}.png")
    render_panel(items, titles, out_path=out,
                 ncols=n_steps, multiview=True)
    return items, titles


def all_axis_trajectories(df_morph: pd.DataFrame,
                          df_scar: pd.DataFrame,
                          df_coords: pd.DataFrame,
                          n_steps: int = TRAJECTORY_N_STEPS,
                          bandwidth: float = TRAJECTORY_BANDWIDTH):
    configs = [
        ('Morph_Axis1', df_morph, 'morph', 'Axis1'),
        ('Scar_Axis1',  df_scar,  'scar',  'Axis1'),
        ('Morph_Axis2', df_morph, 'morph', 'Axis2'),
        ('Scar_Axis2',  df_scar,  'scar',  'Axis2'),
    ]

    all_rows = []
    for axis_col, df_coeffs, label, axis_name in configs:
        result = axis_trajectory(
            df_coeffs=df_coeffs,
            df_coords=df_coords,
            axis_col=axis_col,
            label=label,
            axis_name=axis_name,
            n_steps=n_steps,
            bandwidth=bandwidth,
        )
        if result is None:
            continue
        items, titles = result
        all_rows.append((items, titles, f"{axis_name}\n{label}"))

    if not all_rows:
        print("  [跳过] 没有有效轨迹数据")
        return

    # ── 汇总图：4行 × (n_steps × 3视角) ──────────────────────────────────
    print("\n=== 生成四轨迹汇总对比图（多视角）===")
    views  = list(CAMERA_VIEWS.items())
    nv     = len(views)
    ncols  = n_steps * nv
    nrows  = len(all_rows)

    pl = pv.Plotter(
        off_screen=True,
        shape=(nrows, ncols),
        window_size=[ncols * 350, nrows * 350],
    )
    pl.background_color = 'white'

    for r, (items, titles, row_label) in enumerate(all_rows):
        for c_step, ((mesh, scalars), title) in enumerate(zip(items, titles)):
            smin, smax = scalars.min(), scalars.max()
            for v_idx, (view_name, cam) in enumerate(views):
                col = c_step * nv + v_idx
                pl.subplot(r, col)
                add_mesh_with_wireframe(pl, mesh, scalars, smin, smax)
                t = f"{row_label}\n{title}\n[{view_name}]" \
                    if (c_step == 0 and v_idx == 0) else f"[{view_name}]"
                pl.add_title(t, font_size=7, color='black')
                pl.camera_position = cam
                if cam == 'iso':
                    pl.camera.elevation = 20

        for c_step in range(len(items), n_steps):
            for v_idx in range(nv):
                col = c_step * nv + v_idx
                pl.subplot(r, col)
                pl.background_color = 'white'

    out = os.path.join(OUT_DIR, "recon_trajectory_all_axes.png")
    pl.screenshot(out)
    pl.close()
    print(f"  已保存：{out}")


# ── 主程序 ────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("读取数据...")
    df_morph  = load_csv(MORPH_CSV)
    df_scar   = load_csv(SCAR_CSV)
    df_coords = load_csv(COORD_CSV)
    df_coords = df_coords[df_coords['ID'].str.startswith('EXP')]

    coeff_cols = [c for c in df_morph.columns if c.startswith('coeff_')]
    n_coeffs   = len(coeff_cols)
    expected   = 2 * (LMAX + 1) ** 2
    print(f"检测到 {n_coeffs} 个系数列，期望 {expected}（lmax={LMAX}）")
    if n_coeffs != expected:
        actual_lmax = int(np.sqrt(n_coeffs / 2)) - 1
        print(f"[警告] 实际 lmax 可能是 {actual_lmax}，"
              f"请修改脚本顶部 LMAX = {actual_lmax}")

    df_morph_typed = df_morph.merge(
        df_scar[['ID', 'Typology']].drop_duplicates(),
        on='ID', how='left'
    )

    extreme_ids = list(AXIS_EXTREMES.values())

    # ── 形态端 ──────────────────────────────────────────────────────────────
    print("\n" + "="*60)
    print("形态谱重建")
    print("="*60)
    full_reconstruction(df_morph, extreme_ids, label='morph')
    exp_individual_reconstruction(df_morph, label='morph')
    per_degree_reconstruction(df_morph, extreme_ids,
                              label='morph', degrees=DEGREES_TO_SHOW)
    axis_extremes_panel(df_morph, label='morph')
    typology_mean_reconstruction(df_morph_typed, label='morph')
    im_reconstruction_panel(df_morph, label='morph')
    sdg_reconstruction_panel(df_morph, label='morph')
    typology_mean_with_im(df_morph_typed, label='morph')

    # ── 疤痕端 ──────────────────────────────────────────────────────────────
    print("\n" + "="*60)
    print("疤痕谱重建")
    print("="*60)
    full_reconstruction(df_scar, extreme_ids, label='scar')
    exp_individual_reconstruction(df_scar, label='scar')
    per_degree_reconstruction(df_scar, extreme_ids,
                              label='scar', degrees=DEGREES_TO_SHOW)
    axis_extremes_panel(df_scar, label='scar')
    typology_mean_reconstruction(df_scar, label='scar')
    im_reconstruction_panel(df_scar, label='scar')
    sdg_reconstruction_panel(df_scar, label='scar')
    typology_mean_with_im(df_scar, label='scar')

    # ── 轴连续变化轨迹 ──────────────────────────────────────────────────────
    print("\n" + "="*60)
    print("轴连续变化轨迹")
    print("="*60)
    df_morph_exp = df_morph[df_morph['ID'].str.startswith('EXP')]
    df_scar_exp  = df_scar[df_scar['ID'].str.startswith('EXP')]

    all_axis_trajectories(
        df_morph=df_morph_exp,
        df_scar=df_scar_exp,
        df_coords=df_coords,
        n_steps=TRAJECTORY_N_STEPS,
        bandwidth=TRAJECTORY_BANDWIDTH,
    )

    print("\n全部完成。输出目录：", OUT_DIR)
