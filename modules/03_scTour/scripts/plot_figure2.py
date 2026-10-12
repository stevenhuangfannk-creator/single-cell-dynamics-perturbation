"""读取保存的真实推断结果，独立重绘 Figure 2 a–i（没有数据的面板不画）。"""
import argparse
import json
import logging
from pathlib import Path

import anndata as ad
from anndata.experimental import read_elem
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import sctour as sct

PALETTES = {
    'A': {'nIPC': '#AC3CFF', 'Neuroblast_1': '#FEBB72', 'Neuroblast_2': '#FF7D00',
          'Granule-immature': '#9CC9E2', 'Granule-mature': '#2584C0'},
    'B': {'nIPC': '#AC3CFF', 'Neuroblast': '#FF7D00', 'Immature-Pyr': '#E01F25',
          'Immature-GC': '#9CC9E2', 'GC-juv': '#24BFD5', 'GC-adult': '#2584C0'},
}
TITLES = {'a': '单一颗粒细胞谱系：原始空间与 scTour 动力学', 'b': '多分化谱系的 PCA UMAP',
          'c': '相同 PCA UMAP 上的样本批次', 'd': 'PCA UMAP 上的发育拟时序',
          'e': '拟时序排序及批次和细胞类型', 'f': 'scTour 潜空间及批次 inset',
          'g': '潜空间中的官方向量场流线', 'h': '训练细胞比例的排序稳健性',
          'i': '20% 训练模型对完整谱系的重建'}
DEFAULT_DPI = 600
SUBSAMPLING_PERCENTAGES = {1, 5, 10, 20, 40, 60, 80, 95}

def style(cfg):
    logging.getLogger('fontTools.subset').setLevel(logging.WARNING)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': cfg['font_size'], 'axes.titlesize': 7,
                         'legend.fontsize': cfg['legend_size'], 'svg.fonttype': 'none',
                         'pdf.fonttype': 42, 'figure.facecolor': 'white', 'axes.facecolor': 'white'})

def clean(ax, title=None):
    ax.set_axis_off()
    ax.set_aspect('equal', adjustable='datalim')
    if title:
        ax.set_title(title, pad=3)

def load_plot_data(path):
    # Plotting needs observations/embeddings, not duplicate full count matrices in eight models.
    with h5py.File(path, 'r') as source:
        obs = read_elem(source['obs'])
        obsm = {key: read_elem(source['obsm'][key]) for key in
                ['X_TNODE', 'X_VF', 'X_umap_pca', 'X_umap_latent']}
        reverse = bool(read_elem(source['uns']['manual_reverse']))
    data = ad.AnnData(obs=obs, obsm=obsm)
    data.uns['manual_reverse'] = reverse
    configuration = json.loads(Path(path).with_name('config_used.json').read_text(encoding='utf-8'))
    data.uns['training_percent'] = configuration['percent']
    return data

def categorical(ax, data, key, embedding, palette, cfg, title=None):
    unknown = set(data.obs[key].astype(str)) - set(palette)
    assert not unknown, f'Figure 2 palette does not cover labels: {sorted(unknown)}'
    x = data.obsm[embedding]
    for label, col in palette.items():
        mask = data.obs[key].astype(str).to_numpy() == label
        ax.scatter(x[mask, 0], x[mask, 1], c=col, s=cfg['point_size'], alpha=cfg['point_alpha'],
                   linewidths=0, rasterized=True)
    clean(ax, title)

def pseudotime(ax, data, embedding, cfg, title=None, stream=False):
    if stream:
        # Recompute cosine similarities in the complete 5-D dynamics space and project using neighbors.
        sct.vf.plot_vector_field(data, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE',
                                 E_key=embedding, t_key='ptime', reverse=bool(data.uns['manual_reverse']),
                                 stream=True, stream_color='k', stream_density=cfg['stream_density'],
                                 stream_linewidth=cfg['stream_linewidth'], stream_arrowsize=cfg['stream_arrowsize'],
                                 color='ptime', color_map='viridis', vmin=0, vmax=1, colorbar_loc=None,
                                 size=cfg['point_size'], alpha=cfg['point_alpha'], frameon=False,
                                 show=False, ax=ax, title='')
        # Official streamplot interpolates masked grid speeds into some NaN linewidths.
        # PDF rejects these styling values; suppress only undefined-width marks, preserving all finite geometry/VF.
        for artist in ax.collections:
            widths = np.asarray(artist.get_linewidths())
            if not np.isfinite(widths).all():
                artist.set_linewidths(np.where(np.isfinite(widths), widths, 0))
        for artist in ax.patches:
            if not np.isfinite(artist.get_linewidth()):
                artist.set_linewidth(0)
    else:
        x = data.obsm[embedding]
        ax.scatter(x[:, 0], x[:, 1], c=data.obs.ptime, cmap='viridis', vmin=0, vmax=1,
                   s=cfg['point_size'], alpha=cfg['point_alpha'], linewidths=0, rasterized=True)
    clean(ax, title)

def batch_palette(data):
    batches = sorted(data.obs.batch.unique())
    cmap = plt.get_cmap('tab20')
    return {b: colors.to_hex(cmap(i % 20)) for i, b in enumerate(batches)}

def inset(ax, data, key, palette, cfg):
    # Inset uses the same latent coordinates; no separate layout or batch correction.
    sub = ax.inset_axes([.70, .74, .27, .25])
    categorical(sub, data, key, 'X_umap_latent', palette, {**cfg, 'point_size': .6})

def ordering(ax, data, cfg):
    unknown = set(data.obs.celltype.astype(str)) - set(PALETTES['B'])
    assert not unknown, f'Figure 2 palette does not cover labels: {sorted(unknown)}'
    ix = np.argsort(data.obs.ptime.to_numpy(), kind='stable')
    pt = data.obs.ptime.to_numpy()[ix]
    cmap = plt.get_cmap('viridis')
    batch = batch_palette(data)
    celltypes = data.obs.celltype.astype(str).to_numpy()[ix]
    batches = data.obs.batch.astype(str).to_numpy()[ix]
    # 1 raster column per observed cell; every cell is displayed in stable rank order.
    tracks = [cmap(pt)[:, :3], np.array([colors.to_rgb(batch[b]) for b in batches])]
    cats = list(PALETTES['B'])
    for c in cats:
        rgb = np.ones((len(ix), 3))
        rgb[celltypes == c] = colors.to_rgb(PALETTES['B'][c])
        tracks.append(rgb)
    ax.imshow(np.stack(tracks), aspect='auto', interpolation='nearest', extent=(0, len(ix), len(tracks), 0))
    ax.set_yticks(np.arange(len(tracks)) + .5, ['Pseudotime', 'Batch'] + cats, fontsize=6)
    ax.set_xticks([0, len(ix)], ['0', '1'])
    ax.set_xlabel('Pseudotime rank (all cells)', fontsize=6)
    ax.spines[['top', 'right']].set_visible(False)

def robustness(ax, runs, cfg):
    rows = []
    fractions = []
    for frac, d in sorted(runs.items()):
        ix = np.argsort(d.ptime.to_numpy(), kind='stable')
        rows.append(np.array([colors.to_rgb(PALETTES['B'][v]) for v in d.celltype.astype(str).to_numpy()[ix]]))
        fractions.append(f'{frac:g}%')
    ax.imshow(np.stack(rows), aspect='auto', interpolation='nearest')
    ax.set_yticks(np.arange(len(rows)), fractions, fontsize=6)
    ax.set_xticks([0, len(rows[0]) - 1], ['0', '1'])
    ax.set_xlabel('Pseudotime rank', fontsize=6)
    ax.set_ylabel('Training cells', fontsize=6)
    ax.spines[['top', 'right']].set_visible(False)

def legend(fig, palette, y):
    handles = [Line2D([], [], marker='o', color=c, linestyle='', markersize=3, label=k)
               for k, c in palette.items()]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, y),
               ncol=len(handles), frameon=False, handletextpad=.2, columnspacing=.75, fontsize=6)

def save(fig, path, cfg):
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(path.with_suffix('.' + ext), dpi=cfg.get('output_dpi', DEFAULT_DPI), facecolor='white')
    plt.close(fig)

def plot(a_path, b_path, output, config, panel='all', subsample_dir=None, prefix='scTour'):
    cfg = json.loads(Path(config).read_text(encoding='utf-8'))
    style(cfg)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    A = load_plot_data(a_path) if a_path else None
    B = load_plot_data(b_path) if b_path else None
    runs = {}
    if subsample_dir:
        for folder in Path(subsample_dir).glob('B_*'):
            p = folder / 'inference.tsv.gz'
            summary = folder / 'run_summary.json'
            if p.exists() and summary.exists():
                assert json.loads(summary.read_text(encoding='utf-8'))['epochs'] == 400
                runs[float(folder.name.split('_')[1])] = pd.read_csv(p, sep='\t', usecols=['cell_id', 'celltype', 'ptime'], index_col=0)
    if B is not None and not runs:
        # A main 20% model must never be relabelled as a 95% reference.
        runs[float(B.uns['training_percent']) * 100] = B.obs
    i_path = Path(subsample_dir) / 'B_20' / 'inferred.h5ad' if subsample_dir else None
    if i_path and i_path.exists():
        I = B if b_path and Path(b_path).resolve() == i_path.resolve() else load_plot_data(i_path)
    else:
        I = None
    available = ([] if A is None else ['a']) + ([] if B is None else list('bcdefg'))
    if SUBSAMPLING_PERCENTAGES.issubset(runs):
        available.append('h')
    if I is not None:
        available.append('i')
    requested = available if panel == 'all' else [panel]
    if not requested:
        raise ValueError('Provide at least one actual inference object; no figures or manifest were changed')
    if any(p not in available for p in requested):
        raise ValueError(f'No actual inference for requested panels: {set(requested)-set(available)}')
    def draw(p, axes):
        ax = axes[0] if isinstance(axes, list) else axes
        if p == 'a':
            categorical(axes[0], A, 'celltype', 'X_umap_pca', PALETTES['A'], cfg, 'PCA UMAP')
            pseudotime(axes[1], A, 'X_umap_pca', cfg, 'Pseudotime')
            pseudotime(axes[2], A, 'X_umap_pca', cfg, 'Vector field', stream=True)
            categorical(axes[3], A, 'celltype', 'X_umap_latent', PALETTES['A'], cfg, 'Latent space')
            pseudotime(axes[4], A, 'X_umap_latent', cfg, 'Pseudotime')
            pseudotime(axes[5], A, 'X_umap_latent', cfg, 'Vector field', stream=True)
        elif p == 'b':
            categorical(ax, B, 'celltype', 'X_umap_pca', PALETTES['B'], cfg)
        elif p == 'c':
            categorical(ax, B, 'batch', 'X_umap_pca', batch_palette(B), cfg)
        elif p == 'd':
            pseudotime(ax, B, 'X_umap_pca', cfg)
        elif p == 'e':
            ordering(ax, B, cfg)
        elif p == 'f':
            categorical(ax, B, 'celltype', 'X_umap_latent', PALETTES['B'], cfg)
            inset(ax, B, 'batch', batch_palette(B), cfg)
        elif p == 'g':
            pseudotime(ax, B, 'X_umap_latent', cfg, stream=True)
        elif p == 'h':
            robustness(ax, runs, cfg)
        elif p == 'i':
            pseudotime(ax, I, 'X_umap_latent', cfg, stream=True)
            inset(ax, I, 'celltype', PALETTES['B'], cfg)
        ax.text(-.05, 1.03, p, transform=ax.transAxes, fontsize=9, fontweight='bold')
    entries = []
    for p in requested:
        if p == 'a':
            fig, axs = plt.subplots(1, 6, figsize=(12, 2.6))
            draw(p, list(axs))
            fig.subplots_adjust(bottom=.17, wspace=.12)
            legend(fig, PALETTES['A'], .01)
        else:
            fig, ax = plt.subplots(figsize=(3.3, 3.1))
            fig.subplots_adjust(left=.23 if p == 'e' else .1, right=.97, top=.93, bottom=.29 if p == 'i' else .17)
            draw(p, ax)
            if p in 'bfhi':
                # Individual plots include a compact figure legend; composite uses shared legends.
                handles = [Line2D([], [], marker='o', color=c, linestyle='', markersize=3, label=k)
                           for k, c in PALETTES['B'].items()]
                fig.legend(handles=handles, loc='lower center', ncol=3, fontsize=5, frameon=False)
        if p in 'adgi':
            cax = fig.add_axes(cfg['colorbar_position'].get(p, cfg['colorbar_position']['default']))
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(0, 1), cmap='viridis'), cax=cax,
                              orientation='horizontal', ticks=[0, 1])
            cb.outline.set_visible(False)
            cax.tick_params(labelsize=5)
            cax.set_xlabel('Pseudotime', fontsize=6, labelpad=1)
        stem = out / f'{prefix}_Figure2{p}'
        save(fig, stem, cfg)
        base = 'modules/03_scTour/figures/' + stem.name
        validity = '真实原始 counts、官方模型与邻域投影；训练比例与方向审查见报告，非原论文数值逐点一致。'
        if p == 'a':
            stage_medians = A.obs.groupby('celltype', observed=True).ptime.median().reindex(PALETTES['A'])
            if not np.all(np.diff(stage_medians.to_numpy()) > 0):
                validity = '真实A模型未恢复单调的颗粒细胞成熟顺序，异常完整保留；不能将该面板解释为发育趋势复现成功。'
        entries.append({'id': f'scTour-2{p}', 'title_cn': TITLES[p], 'meaning': TITLES[p],
                        'dataset': 'GSE95315' if p == 'a' else 'GSE104323',
                        'reference': 'https://doi.org/10.1186/s13059-023-02988-9',
                        'analysis_script': 'modules/03_scTour/scripts/run_sctour.py',
                        'plot_script': 'modules/03_scTour/scripts/plot_figure2.py',
                        'png': base+'.png', 'pdf': base+'.pdf', 'svg': base+'.svg',
                        'status': 'PARTIAL', 'validity': validity})
    if prefix != 'scTour':
        # Diagnostic redraw uses preserved filenames and does not replace the main figure index.
        return
    if panel == 'all' and A is not None and B is not None and set(available) == set('abcdefghi'):
        fig = plt.figure(figsize=(10, 8.2))
        gs = fig.add_gridspec(3, 4, height_ratios=[1.05, 1.45, 1.45], hspace=.21, wspace=.22,
                             left=.04, right=.99, top=.96, bottom=.1)
        top = gs[0, :].subgridspec(1, 6, wspace=.13)
        draw('a', [fig.add_subplot(top[0, j]) for j in range(6)])
        mapping = {'b': (1, 0), 'c': (2, 0), 'd': (1, 1), 'e': (2, 1),
                   'f': (1, 2), 'g': (2, 2), 'h': (1, 3), 'i': (2, 3)}
        for p, position in mapping.items():
            ax = fig.add_subplot(gs[position])
            draw(p, ax)
        legend(fig, PALETTES['A'], .722)
        legend(fig, PALETTES['B'], .025)
        cax = fig.add_axes(cfg['colorbar_position']['composite'])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(0, 1), cmap='viridis'), cax=cax,
                          orientation='horizontal', ticks=[0, 1])
        cb.outline.set_visible(False)
        cax.set_xlabel('Pseudotime', fontsize=6)
        save(fig, out / 'scTour_Figure2_Reproduction', cfg)
        base = 'modules/03_scTour/figures/scTour_Figure2_Reproduction'
        entries.append({'id':'scTour-Figure2', 'title_cn':'scTour Figure 2 a–i 完整真实重算组合图',
                        'meaning':'联合比较原始空间、拟时序、潜空间、向量场与训练抽样稳定性。',
                        'dataset':'GSE95315 + GSE104323', 'reference':'https://doi.org/10.1186/s13059-023-02988-9',
                        'analysis_script':'modules/03_scTour/scripts/run_subsampling.py',
                        'plot_script':'modules/03_scTour/scripts/plot_figure2.py',
                        'redraw_command': '& ./work/sctour/.venv/Scripts/python.exe ./modules/03_scTour/scripts/plot_figure2.py --a-path work/sctour/full_A20/inferred.h5ad --b-path work/sctour/subsamples/B_20/inferred.h5ad --subsample-dir work/sctour/subsamples --config modules/03_scTour/config/plot.json --output modules/03_scTour/figures',
                        'png':base+'.png', 'pdf':base+'.pdf', 'svg':base+'.svg', 'status':'PARTIAL',
                        'validity':'所有面板均来自实际400epochs模型；生物学偏差与原论文未确认配置逐项保留，不以视觉一致宣布严格复现。'})
    manifest = out / 'figure_manifest.json'
    if manifest.exists():
        previous = json.loads(manifest.read_text(encoding='utf-8'))
        keep = [e for e in previous if '/diagnostics/' in e['png'] or panel != 'all']
        entries = [e for e in keep if e['id'] not in {n['id'] for n in entries}] + entries
    for entry in entries:
        if entry['id'] in {f'scTour-2{p}' for p in 'abcdefghi'}:
            p = entry['id'][-1]
            inputs = ('--a-path work/sctour/full_A20/inferred.h5ad' if p == 'a' else
                      '--subsample-dir work/sctour/subsamples' if p in 'hi' else
                      '--b-path work/sctour/subsamples/B_20/inferred.h5ad')
            entry['redraw_command'] = ('& ./work/sctour/.venv/Scripts/python.exe ./modules/03_scTour/scripts/plot_figure2.py '
                                       f'{inputs} --panel {p} --config modules/03_scTour/config/plot.json --output modules/03_scTour/figures')
        elif '/diagnostics/' in entry['png']:
            stem = Path(entry['png']).stem
            p = stem[-1]
            prefix_used = stem.removesuffix(f'_Figure2{p}')
            inputs = ('--a-path work/sctour/full_A/inferred.h5ad' if p == 'a' else
                      '--b-path work/sctour/subsamples/B_95/inferred.h5ad')
            entry['redraw_command'] = ('& ./work/sctour/.venv/Scripts/python.exe ./modules/03_scTour/scripts/plot_figure2.py '
                                       f'{inputs} --panel {p} --prefix {prefix_used} --config modules/03_scTour/config/plot.json '
                                       '--output modules/03_scTour/figures/diagnostics')
    manifest.write_text(json.dumps(entries, ensure_ascii=False, indent=2), encoding='utf-8')

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--a-path')
    p.add_argument('--b-path')
    p.add_argument('--output', required=True)
    p.add_argument('--config', required=True)
    p.add_argument('--panel', default='all', choices=['all'] + list('abcdefghi'))
    p.add_argument('--subsample-dir')
    p.add_argument('--prefix', default='scTour')
    plot(**vars(p.parse_args()))
